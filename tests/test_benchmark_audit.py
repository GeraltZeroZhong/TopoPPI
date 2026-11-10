import csv
import gzip
import hashlib
import json
import os
import shutil
import sys
import tempfile
import threading
import unittest
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import trimesh

from topoppi.atlas.uv import set_uv_layout
from topoppi.benchmarking import BenchmarkRunner
from topoppi.benchmarking.coordinate_audit import (
    AUDIT_PROTOCOL,
    AUDIT_SCHEMA_VERSION,
    validate_coordinate_audit,
)
from topoppi.benchmarking.manifest_metadata import INFERENCE_DEPENDENCY_BASIS
from topoppi.benchmarking.runner import _afdb_complex_confidence, _worker_cpu_time
from topoppi.config import BenchmarkConfig, OptCutsConfig
from topoppi.errors import ConfigurationError
from topoppi.io.io_loader import PDBLoader
from topoppi.mesh.provenance import provenance_summary

FIXTURE_DIR = Path(__file__).parent / "fixtures"
TINY_PDB = FIXTURE_DIR / "tiny_complex.pdb"
OPTCUTS_BIN = Path(__file__).parents[1] / "tools" / "OptCuts" / "OptCuts_bin"
OPTCUTS_SHA256 = "d7990fc4f1ca46e0ba06b70801b64701dfdeb795f7efee6f7b9f197aa3b426eb"


class BenchmarkAuditTests(unittest.TestCase):
    def test_residue_objective_uses_declared_prolif_partner_degrees(self):
        loader = PDBLoader(TINY_PDB)
        coords_a, atoms_a = loader.get_chain_data("A")
        coords_b, atoms_b = loader.get_chain_data("B")
        runner = BenchmarkRunner(BenchmarkConfig(str(FIXTURE_DIR), "unused"))
        input_sha256 = hashlib.sha256(TINY_PDB.read_bytes()).hexdigest()

        labels, interaction_weights, objective_weights, source, _definition = runner._residue_objective(
            job_metadata={
                "input_sha256": input_sha256,
                "prolif_file": str(FIXTURE_DIR / "prolif_interactions.json"),
            },
            chain_a="A",
            chain_b="B",
            atoms_a=atoms_a,
            coords_a=coords_a,
            atoms_b=atoms_b,
            coords_b=coords_b,
        )
        _, geometric_weights, _geometric_objective, geometric_source, _ = runner._residue_objective(
            job_metadata={"input_sha256": input_sha256},
            chain_a="A",
            chain_b="B",
            atoms_a=atoms_a,
            coords_a=coords_a,
            atoms_b=atoms_b,
            coords_b=coords_b,
        )

        self.assertEqual(source, "prolif")
        self.assertEqual(geometric_source, "geometric_fallback")
        self.assertEqual(interaction_weights, {"A:GLY:1": 1.0})
        self.assertNotEqual(interaction_weights, geometric_weights)
        self.assertEqual(objective_weights["A:GLY:1"], 2.0)
        self.assertEqual(objective_weights["A:ALA:2"], 1.0)
        self.assertEqual(set(labels), {"A:GLY:1", "A:ALA:2"})

    def test_formal_residue_objective_requires_declared_prolif(self):
        runner = BenchmarkRunner(BenchmarkConfig(str(FIXTURE_DIR), "unused"))
        runner.config = replace(runner.config, formal_mode=True)

        with self.assertRaisesRegex(ValueError, "require prolif_file and prolif_sha256"):
            runner._interaction_job_metadata(
                {},
                structure_path=str(TINY_PDB),
                chain_a="A",
                chain_b="B",
                input_sha256=hashlib.sha256(TINY_PDB.read_bytes()).hexdigest(),
                atoms_a=[],
                atoms_b=[],
            )

    def test_interaction_evidence_excludes_unresolved_partner_residues(self):
        loader = PDBLoader(TINY_PDB)
        _coords_a, atoms_a = loader.get_chain_data("A")
        _coords_b, atoms_b = loader.get_chain_data("B")
        input_sha256 = hashlib.sha256(TINY_PDB.read_bytes()).hexdigest()
        runner = BenchmarkRunner(BenchmarkConfig(str(FIXTURE_DIR), "unused"))
        with tempfile.TemporaryDirectory() as tmp:
            interactions = Path(tmp) / "interactions.json"
            interactions.write_text(
                json.dumps(
                    {
                        "engine": "prolif",
                        "chain_a": "A",
                        "chain_b": "B",
                        "source_sha256": input_sha256,
                        "interactions": [
                            {"res_a_seq": "1", "res_b_seq": "1", "interaction": "VdWContact"},
                            {"res_a_seq": "1", "res_b_seq": "999", "interaction": "HBAcceptor"},
                        ],
                    }
                ),
                encoding="utf-8",
            )

            resolved, _provenance = runner._load_declared_interactions(
                interactions,
                expected_chain_a="A",
                expected_chain_b="B",
                expected_input_sha256=input_sha256,
                atoms_a=atoms_a,
                atoms_b=atoms_b,
            )

        self.assertEqual(len(resolved), 1)
        self.assertIn("1--1--VdWContact", next(iter(resolved)))

    def test_formal_preflight_rejects_zero_resolvable_interaction_pairs(self):
        loader = PDBLoader(TINY_PDB)
        _coords_a, atoms_a = loader.get_chain_data("A")
        _coords_b, atoms_b = loader.get_chain_data("B")
        input_sha256 = hashlib.sha256(TINY_PDB.read_bytes()).hexdigest()
        runner = BenchmarkRunner(BenchmarkConfig(str(FIXTURE_DIR), "unused"))
        runner.config = replace(runner.config, formal_mode=True)
        with tempfile.TemporaryDirectory() as tmp:
            interactions = Path(tmp) / "interactions.json"
            interactions.write_text(
                json.dumps(
                    {
                        "engine": "prolif",
                        "chain_a": "A",
                        "chain_b": "B",
                        "source_sha256": input_sha256,
                        "interactions": [{"res_a_seq": "1", "res_b_seq": "999", "interaction": "HBAcceptor"}],
                    }
                ),
                encoding="utf-8",
            )
            interaction_sha256 = hashlib.sha256(interactions.read_bytes()).hexdigest()

            with self.assertRaisesRegex(ValueError, "no residue pair"):
                runner._interaction_job_metadata(
                    {
                        "prolif_file": str(interactions),
                        "prolif_sha256": interaction_sha256,
                    },
                    structure_path=str(TINY_PDB),
                    chain_a="A",
                    chain_b="B",
                    input_sha256=input_sha256,
                    atoms_a=atoms_a,
                    atoms_b=atoms_b,
                )

    def test_manifest_template_example_matches_its_header(self):
        template = Path(__file__).parents[1] / "docs" / "benchmark_manifest_template.csv"
        with template.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertNotIn(None, row)
        self.assertEqual(row["hotspot_residues_a"], "A:ALA:1;A:PHE:2")
        self.assertEqual(row["prolif_file"], "example.prolif.json")
        self.assertEqual(row["prolif_sha256"], "<64-character-sha256>")
        self.assertEqual(row["record_id"], "record-0001")

    def test_worker_cpu_time_includes_reaped_child_usage(self):
        usage = SimpleNamespace(ru_utime=3.25, ru_stime=0.5)
        with (
            patch("topoppi.benchmarking.runner.time.process_time", return_value=2.0),
            patch("topoppi.benchmarking.runner.resource") as mocked_resource,
        ):
            mocked_resource.RUSAGE_CHILDREN = 7
            mocked_resource.getrusage.return_value = usage
            measured = _worker_cpu_time()

        self.assertEqual(measured, 5.75)

    def test_supervisor_timeout_is_returned_as_a_right_censored_measurement(self):
        class HangingProcess:
            pid = 424242
            returncode = -15

            @staticmethod
            def poll():
                return None

            @staticmethod
            def wait():
                return -15

        with tempfile.TemporaryDirectory() as tmp:
            config = BenchmarkConfig(
                str(FIXTURE_DIR),
                str(Path(tmp) / "out"),
                repetitions=1,
                worker_timeout_sec=1.0,
            )
            runner = BenchmarkRunner(config)
            result_path = Path(config.output_root) / config.worker_log_folder / "example.pdb.measured.000.result.json"
            result_path.parent.mkdir(parents=True)
            result_path.write_text('{"status":"ok","result":{"status":"ok"}}', encoding="utf-8")
            with (
                patch.object(runner, "_worker_cpu_affinity", return_value=[]),
                patch.object(runner, "_terminate_process_tree"),
                patch.object(runner, "_process_tree_rss_mb", return_value=12.0),
                patch("topoppi.benchmarking.runner.subprocess.Popen", return_value=HangingProcess()),
                patch("topoppi.benchmarking.runner.time.perf_counter", side_effect=[0.0, 2.0, 2.5]),
            ):
                outcome = runner._execute_worker(
                    {"pdb": "example.pdb", "chain_a": "A", "chain_b": "B"},
                    run_index=0,
                    is_warmup=False,
                )

        self.assertEqual(outcome["payload"]["status"], "failed")
        self.assertIn("timed out", outcome["payload"]["error"])
        self.assertTrue(outcome["measurement"]["right_censored"])
        self.assertFalse(outcome["measurement"]["worker_completed"])
        self.assertEqual(outcome["measurement"]["termination_reason"], "timeout")
        self.assertEqual(outcome["measurement"]["censoring_threshold_sec"], 1.0)
        self.assertEqual(outcome["measurement"]["censoring_event_elapsed_sec"], 2.0)
        self.assertEqual(outcome["measurement"]["runtime_observation_sec"], 2.0)
        self.assertEqual(outcome["measurement"]["wall_sec"], 2.5)

    def test_operational_method_timeout_is_returned_as_a_right_censored_measurement(self):
        class CompletedProcess:
            pid = 424243

            @staticmethod
            def poll():
                return 0

            @staticmethod
            def wait():
                return 0

        with tempfile.TemporaryDirectory() as tmp:
            config = BenchmarkConfig(
                str(FIXTURE_DIR),
                str(Path(tmp) / "out"),
                repetitions=1,
                execution_profile="operational_optcuts",
                optcuts_variants=("residue_aware_optcuts",),
                include_topology_ablation=False,
                optcuts=replace(
                    OptCutsConfig(),
                    residue_fragmentation_weight=1.0,
                    timeout_sec=300.0,
                ),
            )
            runner = BenchmarkRunner(config)
            result_path = Path(config.output_root) / config.worker_log_folder / "example.pdb.measured.000.result.json"
            result_path.parent.mkdir(parents=True)
            worker_payload = json.dumps(
                {
                    "status": "ok",
                    "result": {
                        "pdb": "example.pdb",
                        "status": "failed",
                        "error": "method budget exhausted",
                        "execution_profile": "operational_optcuts",
                        "operational_method": "residue_aware_optcuts",
                        "method_execution": {
                            "residue_aware_optcuts": {
                                "method_arm_time_budget_sec": 300.0,
                                "failures": [
                                    {
                                        "failure_type": "timeout",
                                        "reason": "OptCuts timed out after 300.0s.",
                                    }
                                ],
                            }
                        },
                        "timing": {"end_to_end": {"wall_sec": 300.5}},
                    },
                }
            )

            def completed_process(*_args, **_kwargs):
                result_path.write_text(worker_payload, encoding="utf-8")
                return CompletedProcess()

            with (
                patch.object(runner, "_worker_cpu_affinity", return_value=[]),
                patch.object(runner, "_process_tree_rss_mb", return_value=12.0),
                patch("topoppi.benchmarking.runner.subprocess.Popen", side_effect=completed_process),
                patch("topoppi.benchmarking.runner.time.perf_counter", side_effect=[0.0, 300.7]),
            ):
                outcome = runner._execute_worker(
                    {"pdb": "example.pdb", "chain_a": "A", "chain_b": "B"},
                    run_index=0,
                    is_warmup=False,
                )

        measurement = outcome["measurement"]
        self.assertTrue(measurement["worker_completed"])
        self.assertTrue(measurement["right_censored"])
        self.assertEqual(measurement["termination_reason"], "optcuts_method_timeout")
        self.assertEqual(measurement["censoring_threshold_sec"], 300.0)
        self.assertEqual(measurement["censoring_event_elapsed_sec"], 300.5)
        self.assertEqual(measurement["runtime_observation_sec"], 300.5)
        self.assertEqual(measurement["wall_sec"], 300.7)

    def test_operational_profile_is_single_method_and_performance_only(self):
        valid = BenchmarkConfig(
            str(FIXTURE_DIR),
            "unused",
            execution_profile="operational_optcuts",
            optcuts_variants=("optcuts_automatic",),
            include_topology_ablation=False,
        )
        valid.validate()

        BenchmarkConfig(
            str(FIXTURE_DIR),
            "unused",
            execution_profile="operational_optcuts",
            optcuts_variants=("residue_aware_optcuts",),
            include_topology_ablation=False,
            optcuts=replace(OptCutsConfig(), residue_fragmentation_weight=1.0),
        ).validate()

        with self.assertRaisesRegex(ConfigurationError, "exactly one automatic"):
            BenchmarkConfig(
                str(FIXTURE_DIR),
                "unused",
                execution_profile="operational_optcuts",
                include_topology_ablation=False,
            ).validate()
        with self.assertRaisesRegex(ConfigurationError, "performance-only"):
            BenchmarkConfig(
                str(FIXTURE_DIR),
                "unused",
                benchmark_purpose="quality",
                execution_profile="operational_optcuts",
                optcuts_variants=("optcuts_automatic",),
                include_topology_ablation=False,
            ).validate()

    def test_optcuts_method_arm_budget_is_shared_across_patches(self):
        meshes = []
        for index in range(2):
            mesh = trimesh.Trimesh(
                vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
                faces=np.array([[0, 1, 2]]),
                process=False,
            )
            mesh.metadata["patch_id"] = f"patch_{index:04d}"
            meshes.append(mesh)

        class FailingOptimizer:
            timeout_values = []

            def optimize_patches(self, _patches, **kwargs):
                self.timeout_values.append(kwargs["timeout_sec"])
                raise RuntimeError("fixture execution failure")

        with tempfile.TemporaryDirectory() as tmp:
            runner = BenchmarkRunner(
                BenchmarkConfig(
                    tmp,
                    str(Path(tmp) / "output"),
                    optcuts=replace(OptCutsConfig(), timeout_sec=10.0),
                )
            )
            optimizer = FailingOptimizer()
            with patch(
                "topoppi.benchmarking.runner.time.perf_counter",
                side_effect=[0.0, 0.0, 11.0, 11.0],
            ):
                output, diagnostics = runner._run_optcuts(
                    meshes,
                    initialization="automatic",
                    optimizer=optimizer,
                )

        self.assertEqual(output, [])
        self.assertEqual(optimizer.timeout_values, [10.0])
        self.assertEqual(diagnostics["attempted"], 2)
        self.assertEqual(diagnostics["invoked"], 1)
        self.assertEqual(diagnostics["not_invoked"], 1)
        self.assertEqual(
            [failure["failure_type"] for failure in diagnostics["failures"]],
            ["execution_failure", "arm_budget_exhausted"],
        )
        self.assertTrue(diagnostics["method_arm_budget_exhausted"])

    def test_formal_quality_mode_uses_one_measurement_per_structure(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir = Path(tmp) / "input"
            input_dir.mkdir()
            manifest = Path(tmp) / "manifest.csv"
            manifest.write_text("pdb,chain_a,chain_b\n", encoding="utf-8")
            config = BenchmarkConfig(
                str(input_dir),
                str(Path(tmp) / "out"),
                chain_selection_mode="manifest",
                manifest_path=str(manifest),
                formal_mode=True,
                benchmark_purpose="quality",
                repetitions=1,
                warmup_runs=0,
                max_workers=2,
                threads_per_worker=1,
                optcuts=OptCutsConfig(expected_binary_sha256="a" * 64),
            )

            config.validate()

    def test_formal_mode_accepts_a_declared_multicore_worker(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir = Path(tmp) / "input"
            input_dir.mkdir()
            manifest = Path(tmp) / "manifest.csv"
            manifest.write_text("pdb,chain_a,chain_b\n", encoding="utf-8")
            config = BenchmarkConfig(
                str(input_dir),
                str(Path(tmp) / "out"),
                chain_selection_mode="manifest",
                manifest_path=str(manifest),
                formal_mode=True,
                benchmark_purpose="quality",
                repetitions=1,
                warmup_runs=0,
                max_workers=1,
                threads_per_worker=4,
                optcuts=OptCutsConfig(expected_binary_sha256="a" * 64),
            )

            config.validate()

    def test_worker_affinity_uses_disjoint_physical_first_slots(self):
        config = BenchmarkConfig(
            str(FIXTURE_DIR),
            "unused",
            max_workers=2,
            threads_per_worker=2,
        )
        runner = BenchmarkRunner(config)
        runner._available_worker_cpus = [0, 2, 4, 6]
        barrier = threading.Barrier(2)

        def affinity_for_worker(index):
            barrier.wait()
            return runner._worker_cpu_affinity({"record_id": str(index)})

        with ThreadPoolExecutor(max_workers=2) as executor:
            affinities = list(executor.map(affinity_for_worker, range(2)))

        self.assertEqual({tuple(item) for item in affinities}, {(0, 2), (4, 6)})
        self.assertTrue(set(affinities[0]).isdisjoint(affinities[1]))

    def test_worker_count_is_capped_by_affinity_capacity(self):
        runner = BenchmarkRunner(
            BenchmarkConfig(
                str(FIXTURE_DIR),
                "unused",
                max_workers=8,
                threads_per_worker=2,
            )
        )
        runner._available_worker_cpus = [0, 2, 4, 6]

        self.assertEqual(runner._resolve_worker_count(20), 2)

    def test_unspecified_worker_count_uses_affinity_capacity(self):
        runner = BenchmarkRunner(
            BenchmarkConfig(
                str(FIXTURE_DIR),
                "unused",
                max_workers=None,
                threads_per_worker=2,
            )
        )
        runner._available_worker_cpus = [0, 2, 4, 6]

        self.assertEqual(runner._resolve_worker_count(20), 2)
        self.assertEqual(runner._resolve_worker_count(1), 1)

    @unittest.skipUnless(OPTCUTS_BIN.is_file(), "OptCuts binary is not included in source distributions")
    def test_preflight_is_read_only_and_resolves_the_binary(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = BenchmarkConfig(
                input_folder=str(FIXTURE_DIR),
                output_root=os.path.join(tmp, "out"),
                min_chain_residues=1,
                repetitions=1,
                optcuts=replace(BenchmarkConfig(str(FIXTURE_DIR), tmp).optcuts, optcuts_bin=str(OPTCUTS_BIN)),
            )
            preflight = BenchmarkRunner(config).preflight()

            self.assertTrue(preflight["ready"], msg=preflight["blockers"])
            self.assertGreaterEqual(preflight["accepted_job_count"], 1)
            self.assertEqual(
                preflight["planned_worker_process_count"],
                preflight["accepted_job_count"],
            )
            self.assertFalse(os.path.exists(config.output_root))

    def test_nonempty_unmatched_output_is_blocked_for_every_execution_profile(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_folder = Path(tmp) / "inputs"
            input_folder.mkdir()
            shutil.copy2(TINY_PDB, input_folder / TINY_PDB.name)

            for profile in ("comparative", "operational_optcuts"):
                with self.subTest(profile=profile):
                    output_root = Path(tmp) / f"output-{profile}"
                    output_root.mkdir()
                    marker = output_root / "keep.txt"
                    marker.write_text("user data\n", encoding="utf-8")
                    optcuts = replace(
                        BenchmarkConfig(str(input_folder), str(output_root)).optcuts,
                        optcuts_bin=str(OPTCUTS_BIN),
                        residue_fragmentation_weight=0.0,
                    )
                    config = BenchmarkConfig(
                        input_folder=str(input_folder),
                        output_root=str(output_root),
                        execution_profile=profile,
                        optcuts_variants=("optcuts_automatic",),
                        include_topology_ablation=profile == "comparative",
                        min_chain_residues=1,
                        repetitions=1,
                        optcuts=optcuts,
                    )
                    runner = BenchmarkRunner(config)

                    preflight = runner.preflight()

                    self.assertFalse(preflight["ready"])
                    self.assertEqual(preflight["output_state"]["state"], "nonempty_unmatched")
                    self.assertTrue(any("output_root" in blocker for blocker in preflight["blockers"]))
                    with self.assertRaisesRegex(ValueError, "output_root"):
                        runner.run()
                    self.assertEqual(marker.read_text(encoding="utf-8"), "user data\n")
                    self.assertFalse((output_root / config.report_filename).exists())

    def test_repeated_preflight_rechecks_output_without_reprocessing_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_folder = Path(tmp) / "inputs"
            input_folder.mkdir()
            shutil.copy2(TINY_PDB, input_folder / TINY_PDB.name)
            output_root = Path(tmp) / "output"
            output_root.mkdir()
            marker = output_root / "old-result.txt"
            marker.write_text("old\n", encoding="utf-8")
            optcuts = replace(
                BenchmarkConfig(str(input_folder), str(output_root)).optcuts,
                optcuts_bin=sys.executable,
                residue_fragmentation_weight=0.0,
            )
            runner = BenchmarkRunner(
                BenchmarkConfig(
                    input_folder=str(input_folder),
                    output_root=str(output_root),
                    execution_profile="operational_optcuts",
                    optcuts_variants=("optcuts_automatic",),
                    include_topology_ablation=False,
                    min_chain_residues=1,
                    repetitions=1,
                    optcuts=optcuts,
                )
            )

            with patch.object(
                runner,
                "_prepare_benchmark_jobs",
                wraps=runner._prepare_benchmark_jobs,
            ) as prepare_jobs:
                first = runner.preflight()
                marker.unlink()
                second = runner.preflight()

        self.assertFalse(first["ready"])
        self.assertEqual(first["output_state"]["state"], "nonempty_unmatched")
        self.assertTrue(second["ready"], msg=second["blockers"])
        self.assertEqual(second["output_state"]["state"], "empty")
        prepare_jobs.assert_called_once()

    def test_formal_coordinate_audit_is_bound_to_the_exact_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.csv"
            manifest.write_text("pdb,chain_a,chain_b\n", encoding="utf-8")
            manifest_sha256 = hashlib.sha256(manifest.read_bytes()).hexdigest()
            audit = Path(tmp) / "coordinate_audit.json"
            audit.write_text(
                json.dumps(
                    {
                        "schema_version": AUDIT_SCHEMA_VERSION,
                        "audit_protocol": AUDIT_PROTOCOL,
                        "status": "passed",
                        "coordinate_failure_count": 0,
                        "coordinate_record_count": 0,
                        "manifest_records": {"experimental": 0},
                        "manifest_sha256": {"experimental": manifest_sha256},
                        "coordinate_results": [],
                    }
                ),
                encoding="utf-8",
            )
            audit_sha256 = hashlib.sha256(audit.read_bytes()).hexdigest()
            config = BenchmarkConfig(
                input_folder=tmp,
                output_root=os.path.join(tmp, "out"),
                chain_selection_mode="manifest",
                manifest_path=str(manifest),
                formal_mode=True,
                benchmark_purpose="quality",
                repetitions=1,
                warmup_runs=0,
                coordinate_audit_path=str(audit),
                expected_coordinate_audit_sha256=audit_sha256,
                optcuts=replace(
                    BenchmarkConfig(tmp, os.path.join(tmp, "x")).optcuts,
                    expected_binary_sha256=OPTCUTS_SHA256,
                ),
            )
            runner = BenchmarkRunner(config)

            self.assertEqual(runner._coordinate_audit_preflight()["status"], "validated")
            manifest.write_text("pdb,chain_a,chain_b\nchanged.pdb,A,B\n", encoding="utf-8")
            mismatch = runner._coordinate_audit_preflight()

        self.assertEqual(mismatch["status"], "manifest_mismatch")

    def test_coordinate_audit_requires_exact_per_record_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.csv"
            fields = (
                "record_id,pdb,input_sha256,chain_a,chain_b,structure_type,"
                "sequence_a_sha256,sequence_b_sha256,chain_a_residue_count,"
                "chain_b_residue_count\n"
            )
            manifest.write_text(
                fields + f"record-1,one.pdb,{'a' * 64},A,B,experimental,{'b' * 64},{'c' * 64},10,12\n",
                encoding="utf-8",
            )
            manifest_sha256 = hashlib.sha256(manifest.read_bytes()).hexdigest()
            evidence = {
                "manifest": "experimental",
                "record_id": "record-1",
                "status": "passed",
                "input_sha256": "a" * 64,
                "chain_a": "A",
                "chain_b": "B",
                "structure_type": "experimental",
                "sequence_a_sha256": "b" * 64,
                "sequence_b_sha256": "c" * 64,
                "chain_a_residue_count": 10,
                "chain_b_residue_count": 12,
                "heavy_atom_count": 150,
            }
            payload = {
                "schema_version": AUDIT_SCHEMA_VERSION,
                "audit_protocol": AUDIT_PROTOCOL,
                "status": "passed",
                "coordinate_failure_count": 0,
                "coordinate_record_count": 1,
                "manifest_records": {"experimental": 1},
                "manifest_sha256": {"experimental": manifest_sha256},
                "coordinate_results": [evidence],
            }
            audit = Path(tmp) / "coordinate_audit.json"
            audit.write_text(json.dumps(payload), encoding="utf-8")
            audit_sha256 = hashlib.sha256(audit.read_bytes()).hexdigest()

            self.assertEqual(
                validate_coordinate_audit(audit, audit_sha256, manifest)["status"],
                "validated",
            )
            payload["coordinate_results"][0]["record_id"] = "forged-record"
            audit.write_text(json.dumps(payload), encoding="utf-8")
            forged_sha256 = hashlib.sha256(audit.read_bytes()).hexdigest()

            self.assertEqual(
                validate_coordinate_audit(audit, forged_sha256, manifest)["status"],
                "manifest_mismatch",
            )

    def test_formal_preflight_rejects_a_different_frozen_git_revision(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.csv"
            manifest.write_text("pdb,chain_a,chain_b\n", encoding="utf-8")
            config = BenchmarkConfig(
                input_folder=tmp,
                output_root=os.path.join(tmp, "out"),
                chain_selection_mode="manifest",
                manifest_path=str(manifest),
                formal_mode=True,
                expected_git_commit="0" * 40,
                benchmark_purpose="quality",
                repetitions=1,
                warmup_runs=0,
                optcuts=replace(
                    BenchmarkConfig(tmp, os.path.join(tmp, "x")).optcuts,
                    optcuts_bin=str(OPTCUTS_BIN),
                    expected_binary_sha256=OPTCUTS_SHA256,
                ),
            )

            preflight = BenchmarkRunner(config).preflight()

        self.assertFalse(preflight["ready"])
        self.assertTrue(any("expected_git_commit" in blocker for blocker in preflight["blockers"]))

    def test_formal_mode_rejects_missing_warmup(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.csv"
            manifest.write_text("pdb,chain_a,chain_b,input_sha256,cluster_id\n", encoding="utf-8")
            config = BenchmarkConfig(
                input_folder=tmp,
                output_root=os.path.join(tmp, "out"),
                chain_selection_mode="manifest",
                manifest_path=str(manifest),
                formal_mode=True,
                repetitions=3,
                warmup_runs=0,
                max_workers=1,
                optcuts=replace(
                    BenchmarkConfig(tmp, os.path.join(tmp, "x")).optcuts,
                    expected_binary_sha256=OPTCUTS_SHA256,
                ),
            )
            with self.assertRaises(ConfigurationError):
                BenchmarkRunner(config)

    def test_formal_preflight_detects_input_checksum_mismatch_without_running(self):
        with tempfile.TemporaryDirectory() as tmp:
            shutil.copy(TINY_PDB, Path(tmp) / TINY_PDB.name)
            manifest = Path(tmp) / "manifest.csv"
            manifest.write_text(
                f"pdb,chain_a,chain_b,input_sha256,cluster_id\n{TINY_PDB.name},A,B,{'0' * 64},cluster-1\n",
                encoding="utf-8",
            )
            config = BenchmarkConfig(
                input_folder=tmp,
                output_root=os.path.join(tmp, "out"),
                chain_selection_mode="manifest",
                manifest_path=str(manifest),
                formal_mode=True,
                repetitions=3,
                warmup_runs=1,
                max_workers=1,
                min_chain_residues=1,
                optcuts=replace(
                    BenchmarkConfig(tmp, os.path.join(tmp, "x")).optcuts,
                    optcuts_bin=str(OPTCUTS_BIN),
                    expected_binary_sha256=OPTCUTS_SHA256,
                ),
            )
            preflight = BenchmarkRunner(config).preflight()

            self.assertFalse(preflight["ready"])
            self.assertEqual(preflight["preprocessing"]["integrity_error_count"], 1)
            self.assertFalse(os.path.exists(config.output_root))
            with self.assertRaisesRegex(ValueError, "no jobs were started"):
                BenchmarkRunner(config).run()
            self.assertFalse(os.path.exists(config.output_root))

    def test_formal_preflight_treats_unlisted_unparseable_input_as_fatal(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "unlisted.pdb").write_text("not a structure\n", encoding="utf-8")
            manifest = Path(tmp) / "manifest.csv"
            manifest.write_text(
                "pdb,chain_a,chain_b,input_sha256,cluster_id,dataset_source,"
                "source_accession,license_or_terms,structure_type\n",
                encoding="utf-8",
            )
            config = BenchmarkConfig(
                input_folder=tmp,
                output_root=os.path.join(tmp, "out"),
                chain_selection_mode="manifest",
                manifest_path=str(manifest),
                formal_mode=True,
                repetitions=3,
                warmup_runs=1,
                max_workers=1,
                min_chain_residues=1,
                optcuts=replace(
                    BenchmarkConfig(tmp, os.path.join(tmp, "x")).optcuts,
                    optcuts_bin=str(OPTCUTS_BIN),
                    expected_binary_sha256=OPTCUTS_SHA256,
                ),
            )

            preflight = BenchmarkRunner(config).preflight()

        self.assertFalse(preflight["ready"])
        self.assertEqual(preflight["preprocessing"]["integrity_error_count"], 1)
        self.assertIn("absent from the explicit manifest", preflight["preprocessing"]["skipped"][0]["reason"])

    def test_preflight_reports_duplicate_manifest_records_as_integrity_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.csv"
            manifest.write_text(
                "pdb,chain_a,chain_b\nrepeat.pdb,A,B\nrepeat.pdb,A,B\n",
                encoding="utf-8",
            )
            config = BenchmarkConfig(
                input_folder=tmp,
                output_root=os.path.join(tmp, "out"),
                chain_selection_mode="manifest",
                manifest_path=str(manifest),
                optcuts=replace(
                    BenchmarkConfig(tmp, os.path.join(tmp, "x")).optcuts,
                    optcuts_bin=str(OPTCUTS_BIN),
                ),
            )

            preflight = BenchmarkRunner(config).preflight()

        self.assertFalse(preflight["ready"])
        self.assertEqual(preflight["preprocessing"]["integrity_error_count"], 1)
        self.assertIn("Duplicate benchmark manifest record", preflight["preprocessing"]["skipped"][0]["reason"])

    def test_formal_manifest_rejects_role_swapped_family_splitting(self):
        runner = BenchmarkRunner(BenchmarkConfig(str(FIXTURE_DIR), "unused"))
        first = {
            "record_id": "first",
            "family_id": "ordered-a-b",
            "sequence_cluster_a": "seq-a",
            "sequence_cluster_b": "seq-b",
            "cluster_id": "component",
            "analysis_split": "test",
        }
        swapped = {
            "record_id": "second",
            "family_id": "ordered-b-a",
            "sequence_cluster_a": "seq-b",
            "sequence_cluster_b": "seq-a",
            "cluster_id": "component",
            "analysis_split": "test",
        }

        with self.assertRaisesRegex(ValueError, "unordered partner pair"):
            runner._validate_manifest_cohort({"first": first, "second": swapped})

    def test_formal_manifest_rejects_sequence_cluster_split_leakage(self):
        runner = BenchmarkRunner(BenchmarkConfig(str(FIXTURE_DIR), "unused"))
        development = {
            "record_id": "development",
            "family_id": "family-a",
            "sequence_cluster_a": "shared",
            "sequence_cluster_b": "dev-only",
            "cluster_id": "component-a",
            "analysis_split": "development",
        }
        test = {
            "record_id": "test",
            "family_id": "family-b",
            "sequence_cluster_a": "shared",
            "sequence_cluster_b": "test-only",
            "cluster_id": "component-b",
            "analysis_split": "test",
        }

        with self.assertRaisesRegex(ValueError, "sequence cluster"):
            runner._validate_manifest_cohort({"development": development, "test": test})

    def test_formal_manifest_rejects_sequence_cluster_component_fragmentation(self):
        runner = BenchmarkRunner(BenchmarkConfig(str(FIXTURE_DIR), "unused"))
        runner.config = replace(runner.config, formal_mode=True)
        first = {
            "record_id": "first",
            "family_id": "family-a",
            "sequence_cluster_a": "shared",
            "sequence_cluster_b": "partner-a",
            "cluster_id": "homology-a",
            "analysis_split": "test",
            "analysis_split_component_id": "component-a",
            "analysis_split_basis": "frozen-components",
        }
        second = {
            "record_id": "second",
            "family_id": "family-b",
            "sequence_cluster_a": "shared",
            "sequence_cluster_b": "partner-b",
            "cluster_id": "homology-b",
            "analysis_split": "test",
            "analysis_split_component_id": "component-b",
            "analysis_split_basis": "frozen-components",
        }

        with self.assertRaisesRegex(ValueError, "sequence-cluster component"):
            runner._validate_manifest_cohort({"first": first, "second": second})

    def test_formal_manifest_record_requires_an_explicit_analysis_split(self):
        runner = BenchmarkRunner(BenchmarkConfig(str(FIXTURE_DIR), "unused"))
        record = {
            "record_id": "record-1",
            "input_sha256": "a" * 64,
            "cluster_id": "component",
            "family_id": "family",
            "sequence_cluster_a": "sequence-a",
            "sequence_cluster_b": "sequence-b",
            "dataset_source": "source",
            "source_accession": "accession",
            "license_or_terms": "terms",
            "structure_type": "experimental",
        }
        runner.config = replace(runner.config, formal_mode=True)

        with self.assertRaisesRegex(ValueError, "missing analysis_split"):
            runner._validate_manifest_record(record, actual_sha256="a" * 64)

    def test_formal_manifest_requires_a_unique_nonempty_record_id(self):
        runner = BenchmarkRunner(BenchmarkConfig(str(FIXTURE_DIR), "unused"))
        runner.config = replace(runner.config, formal_mode=True)
        shared = {
            "family_id": "family",
            "sequence_cluster_a": "sequence-a",
            "sequence_cluster_b": "sequence-b",
            "cluster_id": "component",
            "analysis_split": "test",
            "analysis_split_component_id": "split-component",
            "analysis_split_basis": "frozen-components",
        }

        with self.assertRaisesRegex(ValueError, "without record_id"):
            runner._validate_manifest_cohort({"first": dict(shared)})
        with self.assertRaisesRegex(ValueError, "Duplicate formal manifest record_id"):
            runner._validate_manifest_cohort(
                {
                    "first": {**shared, "record_id": "duplicate"},
                    "second": {**shared, "record_id": "duplicate"},
                }
            )

    def test_formal_manifest_rejects_an_unknown_structure_type(self):
        runner = BenchmarkRunner(BenchmarkConfig(str(FIXTURE_DIR), "unused"))
        record = {
            "record_id": "record-1",
            "input_sha256": "a" * 64,
            "cluster_id": "component",
            "family_id": "family",
            "sequence_cluster_a": "sequence-a",
            "sequence_cluster_b": "sequence-b",
            "analysis_split": "test",
            "analysis_split_component_id": "split-component",
            "analysis_split_basis": "frozen-components",
            "dataset_source": "source",
            "source_accession": "accession",
            "license_or_terms": "terms",
            "structure_type": "afbd",
        }
        runner.config = replace(runner.config, formal_mode=True)

        with self.assertRaisesRegex(ValueError, "unsupported structure_type"):
            runner._validate_manifest_record(record, actual_sha256="a" * 64)

    def test_formal_interaction_evidence_requires_semantic_bindings(self):
        runner = BenchmarkRunner(BenchmarkConfig(str(FIXTURE_DIR), "unused"))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "interactions.json"
            path.write_text(
                json.dumps({"interactions": [{"res_a_seq": "1", "res_b_seq": "2"}]}),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "must declare chain_a"):
                runner._load_declared_interactions(
                    path,
                    expected_chain_a="A",
                    expected_chain_b="B",
                    expected_input_sha256="a" * 64,
                    require_bindings=True,
                )

    def test_interaction_evidence_requires_object_root_and_canonical_residue_tokens(self):
        runner = BenchmarkRunner(BenchmarkConfig(str(FIXTURE_DIR), "unused"))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "interactions.json"
            path.write_text(
                json.dumps([{"res_a_seq": "ALA001A", "res_b_seq": "GLY002"}]),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "root must be an object"):
                runner._load_declared_interactions(
                    path,
                    expected_chain_a="A",
                    expected_chain_b="B",
                )

            path.write_text(
                json.dumps(
                    {
                        "interactions": [
                            {
                                "res_a_seq": "ALA001A",
                                "res_b_seq": "GLY002",
                                "interaction": "Hydrophobic",
                            }
                        ]
                    }
                ),
                encoding="utf-8",
            )
            interactions, _provenance = runner._load_declared_interactions(
                path,
                expected_chain_a="A",
                expected_chain_b="B",
            )

        self.assertEqual(interactions, {"1A--2--Hydrophobic--record0": "1A"})

    def test_worker_reuses_preflight_interaction_checksum(self):
        runner = BenchmarkRunner(BenchmarkConfig(str(FIXTURE_DIR), "unused"))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "interactions.json"
            path.write_text(json.dumps({"interactions": []}), encoding="utf-8")
            with patch("topoppi.benchmarking.runner.sha256_file") as sha256:
                _interactions, provenance = runner._load_declared_interactions(
                    path,
                    expected_chain_a="A",
                    expected_chain_b="B",
                    known_file_sha256="a" * 64,
                )

        sha256.assert_not_called()
        self.assertEqual(provenance["sha256"], "a" * 64)

    def test_formal_predicted_structure_requires_valid_plddt_bfactors(self):
        with tempfile.TemporaryDirectory() as tmp:
            structure = Path(tmp) / TINY_PDB.name
            structure.write_text(
                TINY_PDB.read_text(encoding="utf-8").replace(" 20.00", "120.00"),
                encoding="utf-8",
            )
            input_sha256 = hashlib.sha256(structure.read_bytes()).hexdigest()
            manifest = Path(tmp) / "manifest.csv"
            manifest.write_text(
                "record_id,pdb,chain_a,chain_b,input_sha256,cluster_id,family_id,sequence_cluster_a,"
                "sequence_cluster_b,inference_sequence_cluster_a,inference_sequence_cluster_b,"
                "inference_family_id,inference_dependency_basis,analysis_split,analysis_split_component_id,"
                "analysis_split_basis,dataset_source,source_accession,"
                "license_or_terms,structure_type,"
                "confidence_metric,confidence_source\n"
                f"record-1,{structure.name},A,B,{input_sha256},cluster-1,family-1,seq-a,seq-b,pdep-a,pdep-b,"
                f"pifam_pdep-a_pdep-b,{INFERENCE_DEPENDENCY_BASIS},test,split-component,"
                "frozen-components,AFDB,"
                "AF-test,CC-BY-4.0,"
                "predicted,plddt_bfactor,AlphaFold\n",
                encoding="utf-8",
            )
            config = BenchmarkConfig(
                input_folder=tmp,
                output_root=os.path.join(tmp, "out"),
                chain_selection_mode="manifest",
                manifest_path=str(manifest),
                formal_mode=True,
                repetitions=3,
                warmup_runs=1,
                max_workers=1,
                min_chain_residues=1,
                optcuts=replace(
                    BenchmarkConfig(tmp, os.path.join(tmp, "x")).optcuts,
                    optcuts_bin=str(OPTCUTS_BIN),
                    expected_binary_sha256=OPTCUTS_SHA256,
                ),
            )

            preflight = BenchmarkRunner(config).preflight()

        self.assertFalse(preflight["ready"])
        self.assertEqual(preflight["preprocessing"]["integrity_error_count"], 1)
        self.assertIn("finite 0-100 pLDDT", preflight["preprocessing"]["skipped"][0]["reason"])

    def test_formal_predicted_structure_validates_partner_chain_plddt(self):
        with tempfile.TemporaryDirectory() as tmp:
            structure = Path(tmp) / TINY_PDB.name
            lines = []
            for line in TINY_PDB.read_text(encoding="utf-8").splitlines(keepends=True):
                if line.startswith("ATOM") and len(line) > 21 and line[21] == "B":
                    line = f"{line[:60]}120.00{line[66:]}"
                lines.append(line)
            structure.write_text("".join(lines), encoding="utf-8")
            input_sha256 = hashlib.sha256(structure.read_bytes()).hexdigest()
            manifest = Path(tmp) / "manifest.csv"
            manifest.write_text(
                "record_id,pdb,chain_a,chain_b,input_sha256,cluster_id,family_id,sequence_cluster_a,"
                "sequence_cluster_b,inference_sequence_cluster_a,inference_sequence_cluster_b,"
                "inference_family_id,inference_dependency_basis,analysis_split,analysis_split_component_id,"
                "analysis_split_basis,dataset_source,source_accession,"
                "license_or_terms,structure_type,confidence_metric,confidence_source\n"
                f"record-1,{structure.name},A,B,{input_sha256},cluster-1,family-1,seq-a,seq-b,pdep-a,pdep-b,"
                f"pifam_pdep-a_pdep-b,{INFERENCE_DEPENDENCY_BASIS},test,split-component,"
                "frozen-components,AFDB,AF-test,CC-BY-4.0,"
                "predicted,plddt_bfactor,AlphaFold\n",
                encoding="utf-8",
            )
            config = BenchmarkConfig(
                input_folder=tmp,
                output_root=os.path.join(tmp, "out"),
                chain_selection_mode="manifest",
                manifest_path=str(manifest),
                formal_mode=True,
                benchmark_purpose="quality",
                repetitions=1,
                warmup_runs=0,
                max_workers=1,
                min_chain_residues=1,
                optcuts=replace(
                    BenchmarkConfig(tmp, os.path.join(tmp, "x")).optcuts,
                    optcuts_bin=str(OPTCUTS_BIN),
                    expected_binary_sha256=OPTCUTS_SHA256,
                ),
            )

            preflight = BenchmarkRunner(config).preflight()

        self.assertEqual(preflight["preprocessing"]["integrity_error_count"], 1)
        self.assertIn("partner-chain", preflight["preprocessing"]["skipped"][0]["reason"])

    def test_formal_predicted_manifest_propagates_frozen_dependency_and_method_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            structure = Path(tmp) / TINY_PDB.name
            shutil.copy(TINY_PDB, structure)
            input_sha256 = hashlib.sha256(structure.read_bytes()).hexdigest()
            interactions = Path(tmp) / "interactions.json"
            interactions.write_text(
                json.dumps(
                    {
                        "engine": "prolif",
                        "chain_a": "A",
                        "chain_b": "B",
                        "source_sha256": input_sha256,
                        "interactions": [
                            {
                                "res_a_seq": "1",
                                "res_b_seq": "1",
                                "interaction": "Hydrophobic",
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            interactions_sha256 = hashlib.sha256(interactions.read_bytes()).hexdigest()
            manifest = Path(tmp) / "manifest.csv"
            manifest.write_text(
                "record_id,pdb,chain_a,chain_b,input_sha256,cluster_id,family_id,sequence_cluster_a,"
                "sequence_cluster_b,inference_sequence_cluster_a,inference_sequence_cluster_b,"
                "inference_family_id,inference_dependency_basis,analysis_split,analysis_split_component_id,"
                "analysis_split_basis,dataset_source,source_accession,"
                "license_or_terms,structure_type,confidence_metric,confidence_source,"
                "experimental_methods_json,experimental_method_group,experimental_method_contains_nmr,"
                "pdbbind_index_resolution_angstrom,rcsb_resolution_combined_angstrom_json,"
                "rcsb_experiment_metadata_source,prolif_file,prolif_sha256\n"
                f"record-1,{structure.name},A,B,{input_sha256},cluster-1,family-1,seq-a,seq-b,pdep-a,pdep-b,"
                f"pifam_pdep-a_pdep-b,{INFERENCE_DEPENDENCY_BASIS},test,split-component,"
                "frozen-components,AFDB,AF-test,CC-BY-4.0,"
                "predicted,plddt_bfactor,AlphaFold,"
                f'"[""SOLUTION NMR""]",solution_nmr,True,NMR,"[2.0]",rcsb_graphql,'
                f"{interactions.name},{interactions_sha256}\n",
                encoding="utf-8",
            )
            config = BenchmarkConfig(
                input_folder=tmp,
                output_root=os.path.join(tmp, "out"),
                chain_selection_mode="manifest",
                manifest_path=str(manifest),
                formal_mode=True,
                repetitions=3,
                warmup_runs=1,
                max_workers=1,
                min_chain_residues=1,
                optcuts=replace(
                    BenchmarkConfig(tmp, os.path.join(tmp, "x")).optcuts,
                    optcuts_bin=str(OPTCUTS_BIN),
                    expected_binary_sha256=OPTCUTS_SHA256,
                ),
            )

            jobs, preprocessing = BenchmarkRunner(config)._prepare_benchmark_jobs([structure.name])

        self.assertEqual(preprocessing["accepted_files"], 1)
        self.assertEqual(jobs[0]["inference_sequence_cluster_a"], "pdep-a")
        self.assertEqual(jobs[0]["inference_family_id"], "pifam_pdep-a_pdep-b")
        self.assertEqual(jobs[0]["inference_dependency_basis"], INFERENCE_DEPENDENCY_BASIS)
        self.assertEqual(jobs[0]["experimental_methods_json"], '["SOLUTION NMR"]')
        self.assertEqual(jobs[0]["experimental_method_group"], "solution_nmr")
        self.assertEqual(jobs[0]["experimental_method_contains_nmr"], "True")
        self.assertEqual(jobs[0]["pdbbind_index_resolution_angstrom"], "NMR")
        self.assertEqual(jobs[0]["rcsb_resolution_combined_angstrom_json"], "[2.0]")
        self.assertEqual(jobs[0]["rcsb_experiment_metadata_source"], "rcsb_graphql")

    def test_prediction_dependency_clusters_cannot_cross_analysis_splits(self):
        runner = BenchmarkRunner(BenchmarkConfig(str(FIXTURE_DIR), "unused"))
        development = {
            "record_id": "development",
            "family_id": "family-dev",
            "sequence_cluster_a": "sequence-dev-a",
            "sequence_cluster_b": "sequence-dev-b",
            "cluster_id": "component-dev",
            "analysis_split": "development",
            "inference_sequence_cluster_a": "shared-prediction",
            "inference_sequence_cluster_b": "prediction-dev-b",
            "inference_family_id": "pifam_prediction-dev-b_shared-prediction",
            "inference_dependency_basis": INFERENCE_DEPENDENCY_BASIS,
        }
        test = {
            "record_id": "test",
            "family_id": "family-test",
            "sequence_cluster_a": "sequence-test-a",
            "sequence_cluster_b": "sequence-test-b",
            "cluster_id": "component-test",
            "analysis_split": "test",
            "inference_sequence_cluster_a": "shared-prediction",
            "inference_sequence_cluster_b": "prediction-test-b",
            "inference_family_id": "pifam_prediction-test-b_shared-prediction",
            "inference_dependency_basis": INFERENCE_DEPENDENCY_BASIS,
        }

        with self.assertRaisesRegex(ValueError, "prediction-dependency sequence cluster"):
            runner._validate_manifest_cohort({"development": development, "test": test})

    def test_failed_worker_result_preserves_inference_dependencies(self):
        runner = BenchmarkRunner(BenchmarkConfig(str(FIXTURE_DIR), "unused", repetitions=1))
        runner._execute_worker = lambda *_args, **_kwargs: {
            "payload": {"status": "failed", "error": "expected failure"},
            "measurement": {"wall_sec": 0.1, "peak_rss_mb": 1.0},
        }
        job = {
            "pdb": "example.pdb",
            "chain_a": "A",
            "chain_b": "B",
            "inference_sequence_cluster_a": "pdep-a",
            "inference_sequence_cluster_b": "pdep-b",
            "inference_family_id": "pifam_pdep-a_pdep-b",
            "inference_dependency_basis": INFERENCE_DEPENDENCY_BASIS,
        }

        result = runner._run_isolated_job(job)

        self.assertEqual(result["inference_family_id"], "pifam_pdep-a_pdep-b")
        self.assertEqual(result["inference_dependency_basis"], INFERENCE_DEPENDENCY_BASIS)

    def test_manifest_preserves_a_zero_confidence_threshold(self):
        with tempfile.TemporaryDirectory() as tmp:
            structure = Path(tmp) / TINY_PDB.name
            shutil.copy(TINY_PDB, structure)
            manifest = Path(tmp) / "manifest.csv"
            manifest.write_text(
                f"pdb,chain_a,chain_b,confidence_threshold\n{structure.name},A,B,0\n",
                encoding="utf-8",
            )
            config = BenchmarkConfig(
                input_folder=tmp,
                output_root=os.path.join(tmp, "out"),
                chain_selection_mode="manifest",
                manifest_path=str(manifest),
                min_chain_residues=1,
                optcuts=replace(
                    BenchmarkConfig(tmp, os.path.join(tmp, "x")).optcuts,
                    optcuts_bin=str(OPTCUTS_BIN),
                ),
            )
            jobs, preprocessing = BenchmarkRunner(config)._prepare_benchmark_jobs([structure.name])

        self.assertEqual(preprocessing["accepted_files"], 1)
        self.assertEqual(jobs[0]["confidence_threshold"], "0")

    def test_complex_confidence_is_not_attached_to_a_monomer_replacement(self):
        metadata = {
            "structure_type": "afdb_monomer_replacement",
            "afdb_model_id": "AF-P1-F1+AF-P2-F1",
            "afdb_iptm": 0.9,
            "afdb_ipsae": 0.8,
        }

        self.assertEqual(_afdb_complex_confidence(metadata), {})
        self.assertEqual(
            _afdb_complex_confidence({**metadata, "structure_type": "afdb"}),
            {
                "model_id": "AF-P1-F1+AF-P2-F1",
                "iptm": 0.9,
                "ipsae": 0.8,
                "pdockq": None,
                "pdockq2": None,
                "lis": None,
            },
        )

    def test_formal_paired_record_rejects_an_unknown_geometry_stratum(self):
        with tempfile.TemporaryDirectory() as tmp:
            structure = Path(tmp) / TINY_PDB.name
            shutil.copy(TINY_PDB, structure)
            input_sha256 = hashlib.sha256(structure.read_bytes()).hexdigest()
            manifest = Path(tmp) / "manifest.csv"
            manifest.write_text(
                "record_id,pdb,chain_a,chain_b,input_sha256,cluster_id,family_id,sequence_cluster_a,"
                "sequence_cluster_b,analysis_split,analysis_split_component_id,analysis_split_basis,"
                "dataset_source,source_accession,license_or_terms,structure_type,"
                "paired_record_id,paired_geometry_stratum,paired_contact_recall_fnat,"
                "paired_contact_precision,paired_contact_jaccard,"
                "paired_interface_ligand_ca_rmsd_angstrom,"
                "paired_cross_chain_clash_atom_fraction\n"
                f"record-1,{structure.name},A,B,{input_sha256},cluster-1,family-1,seq-a,seq-b,test,"
                "split-component,frozen-components,"
                "source,accession,terms,experimental,pair-1,unknown,0.8,0.8,0.7,1.0,0.01\n",
                encoding="utf-8",
            )
            config = BenchmarkConfig(
                input_folder=tmp,
                output_root=os.path.join(tmp, "out"),
                chain_selection_mode="manifest",
                manifest_path=str(manifest),
                formal_mode=True,
                repetitions=3,
                warmup_runs=1,
                max_workers=1,
                min_chain_residues=1,
                optcuts=replace(
                    BenchmarkConfig(tmp, os.path.join(tmp, "x")).optcuts,
                    optcuts_bin=str(OPTCUTS_BIN),
                    expected_binary_sha256=OPTCUTS_SHA256,
                ),
            )

            preflight = BenchmarkRunner(config).preflight()

        self.assertFalse(preflight["ready"])
        self.assertIn("invalid paired_geometry_stratum", preflight["preprocessing"]["skipped"][0]["reason"])

    def test_optcuts_failure_is_recorded_without_lscm_fallback(self):
        mesh = trimesh.Trimesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2]]),
            process=False,
        )
        mesh.metadata["patch_id"] = "patch_0000"
        set_uv_layout(mesh, np.asarray(mesh.vertices[:, :2]))
        with tempfile.TemporaryDirectory() as tmp:
            config = BenchmarkConfig(
                input_folder=tmp,
                output_root=os.path.join(tmp, "out"),
                repetitions=1,
                optcuts=replace(
                    BenchmarkConfig(tmp, os.path.join(tmp, "x")).optcuts,
                    optcuts_bin=os.path.join(tmp, "missing-optcuts"),
                ),
            )
            output, diagnostic = BenchmarkRunner(config)._run_optcuts([mesh], initialization="provided")

        self.assertEqual(output, [])
        self.assertEqual(diagnostic["success"], 0)
        self.assertFalse(diagnostic["fallback_used"])
        self.assertIn("not found", diagnostic["failures"][0]["reason"])

    def test_incomplete_residue_aware_pair_does_not_emit_partial_domain_efficacy(self):
        record = BenchmarkRunner._residue_aware_ablation_record(
            comparison_complete=False,
            residue_fragmentation_weight=0.5,
            comparisons={"automatic": {"objective_weighted_fragmentation_treatment": 0.1}},
        )

        self.assertEqual(record["status"], "incomplete_comparison")
        self.assertFalse(record["efficacy_values_available"])
        self.assertEqual(record["comparisons"], {})

    def test_written_report_json_uses_null_not_nonstandard_nan_tokens(self):
        with tempfile.TemporaryDirectory() as tmp:
            runner = BenchmarkRunner(BenchmarkConfig(input_folder=tmp, output_root=os.path.join(tmp, "out")))
            os.makedirs(runner.config.output_root)
            runner._write_outputs(
                {
                    "preprocessing": {"accepted": [], "skipped": []},
                    "summary": {"value": float("nan")},
                },
                [],
            )
            text = Path(runner.config.output_root, runner.config.report_filename).read_text(encoding="utf-8")
            parsed = json.loads(text, parse_constant=lambda token: self.fail(f"non-standard JSON token: {token}"))

        self.assertIsNone(parsed["summary"]["value"])

    def test_per_face_evidence_uses_the_full_prepared_reference_domain(self):
        reference = trimesh.Trimesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2]]),
            process=False,
        )
        reference.metadata["patch_id"] = "patch_0000"
        reference.metadata["source_face_ids"] = np.array([17])
        output = reference.copy()
        set_uv_layout(output, np.asarray(output.vertices[:, :2]), key="uv_optcuts")
        runner = BenchmarkRunner(BenchmarkConfig(str(FIXTURE_DIR), "unused"))

        records = runner._per_face_sample_records(
            [reference],
            {"optcuts_automatic": [output]},
        )

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["source_face_id"], 17)
        self.assertIn("optcuts_automatic_distortion", records[0])

    def test_externalized_residue_evidence_retains_exact_pair_domain(self):
        residue = {"residue": "A:GLY:1", "fragmentation": 1.0}
        with tempfile.TemporaryDirectory() as tmp:
            output_root = Path(tmp) / "out"
            output_root.mkdir()
            runner = BenchmarkRunner(BenchmarkConfig(tmp, str(output_root)))
            payload = {
                "status": "ok",
                "result": {
                    "pdb": "pair.pdb",
                    "comparison_domain": {"signature": "standard-domain"},
                    "residue_aware_comparison_domain": {"signature": "topoppi-domain"},
                    "residue_footprint_fragmentation": {"methods": {"optcuts_automatic": {"residues": [residue]}}},
                    "residue_aware_pair_quality": {
                        "complete": True,
                        "domain_signature": "topoppi-domain",
                        "methods": {
                            "optcuts_automatic": {"residue_footprint_fragmentation": {"residues": [residue]}},
                            "residue_aware_optcuts": {"residue_footprint_fragmentation": {"residues": [residue]}},
                        },
                    },
                },
            }
            result_path = output_root / "pair.result.json"
            externalized = runner._externalize_worker_payload(payload, result_path, preserve_details=True)
            result = externalized["result"]
            with gzip.open(
                output_root / result["detail_artifact"]["path"],
                "rt",
                encoding="utf-8",
            ) as handle:
                detail = json.load(handle)

        domains = {record["evidence_domain"] for record in detail["per_residue_records"]}
        self.assertEqual(domains, {"top_level_method_domain", "residue_aware_exact_pair"})
        exact_methods = {
            record["method"]
            for record in detail["per_residue_records"]
            if record["evidence_domain"] == "residue_aware_exact_pair"
        }
        self.assertEqual(exact_methods, {"optcuts_automatic", "residue_aware_optcuts"})

    def test_topology_rejected_components_remain_in_biological_retention_denominator(self):
        loader = PDBLoader(str(TINY_PDB))
        coords_a, atoms_a = loader.get_chain_data("A")
        coords_b, atoms_b = loader.get_chain_data("B")
        with tempfile.TemporaryDirectory() as tmp:
            runner = BenchmarkRunner(BenchmarkConfig(input_folder=tmp, output_root=os.path.join(tmp, "out")))
            records = runner._patch_retention_records(
                [],
                atoms_a,
                coords_a,
                atoms_b,
                coords_b,
                extracted_patches=[],
                topology_components=[
                    {
                        "component_index": 0,
                        "patch_id": "component_0000",
                        "status": "dropped",
                        "reason": "below_min_patch_vertices",
                        "before_sanitation": {
                            "face_count": 1,
                            "vertex_count": 3,
                            "area": 2.0,
                            "source_face_ids": [4],
                            "source_vertex_ids": [1, 2, 3],
                            "source_atom_indices": [0],
                        },
                        "after_sanitation": None,
                    }
                ],
                preparation={},
                job_metadata={"chain_a": "A", "chain_b": "B"},
            )

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["retention_status"], "rejected")
        self.assertEqual(records[0]["face_retention_ratio"], 0.0)
        self.assertEqual(records[0]["residue_retention_ratio"], 0.0)
        self.assertEqual(records[0]["rejection_stage"], "topology_extraction")

    def test_vertex_duplication_is_not_reported_as_source_vertex_retention_above_one(self):
        loader = PDBLoader(str(TINY_PDB))
        coords_a, atoms_a = loader.get_chain_data("A")
        coords_b, atoms_b = loader.get_chain_data("B")
        vertices = np.asarray(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, -1.0, 0.0],
            ],
            dtype=np.float64,
        )
        extracted = trimesh.Trimesh(
            vertices=vertices,
            faces=np.asarray([[0, 1, 2], [0, 3, 4]], dtype=np.int64),
            process=False,
        )
        extracted.metadata.update(
            {
                "patch_id": "patch_0000",
                "source_vertex_ids": np.arange(5, dtype=np.int64),
                "source_face_ids": np.arange(2, dtype=np.int64),
                "source_atom_indices": np.arange(5, dtype=np.int64),
            }
        )
        extracted.metadata["topology_component_before"] = provenance_summary(extracted)
        prepared = trimesh.Trimesh(
            vertices=np.vstack((vertices, vertices[0])),
            faces=np.asarray([[0, 1, 2], [5, 3, 4]], dtype=np.int64),
            process=False,
        )
        prepared.metadata.update(
            {
                "patch_id": "patch_0000",
                "source_vertex_ids": np.asarray([0, 1, 2, 3, 4, 0]),
                "source_face_ids": np.arange(2, dtype=np.int64),
                "source_atom_indices": np.asarray([0, 1, 2, 3, 4, 0]),
            }
        )

        with tempfile.TemporaryDirectory() as tmp:
            runner = BenchmarkRunner(BenchmarkConfig(input_folder=tmp, output_root=os.path.join(tmp, "out")))
            records = runner._patch_retention_records(
                [prepared],
                atoms_a,
                coords_a,
                atoms_b,
                coords_b,
                extracted_patches=[extracted],
                topology_components=[],
                preparation={},
                job_metadata={"chain_a": "A", "chain_b": "B"},
            )

        self.assertEqual(records[0]["source_vertex_retention_ratio"], 1.0)
        self.assertEqual(records[0]["parameterization_source_vertex_retention_ratio"], 1.0)
        self.assertEqual(records[0]["materialized_vertex_count_ratio"], 1.2)

    def test_resume_provenance_is_rebuilt_from_structure_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_root = Path(tmp) / "out"
            output_root.mkdir()
            runner = BenchmarkRunner(BenchmarkConfig(input_folder=tmp, output_root=str(output_root)))
            provenance_path = output_root / runner.config.provenance_filename
            fields = ["pdb", "patch_id", "entity", "final_index", "source_id", "source_atom_index"]
            with gzip.open(provenance_path, "wt", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerow(
                    {
                        "pdb": "old.pdb",
                        "patch_id": "patch_0000",
                        "entity": "face",
                        "final_index": 0,
                        "source_id": 7,
                        "source_atom_index": "",
                    }
                )
            runner._write_provenance_csv(
                [
                    {
                        "pdb": "old.pdb",
                        "provenance_records": [
                            {"patch_id": "duplicate", "entity": "face", "final_index": 0, "source_id": 8}
                        ],
                    },
                    {
                        "pdb": "new.pdb",
                        "provenance_records": [
                            {"patch_id": "patch_0000", "entity": "face", "final_index": 0, "source_id": 9}
                        ],
                    },
                ]
            )
            with gzip.open(provenance_path, "rt", newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))

        self.assertEqual([row["pdb"] for row in rows], ["old.pdb", "new.pdb"])
        self.assertEqual(rows[0]["patch_id"], "duplicate")
        self.assertFalse((output_root / f"{runner.config.provenance_filename}.tmp").exists())

    def test_resume_preserves_a_terminal_failed_attempt(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_root = Path(tmp) / "out"
            output_root.mkdir()
            runner = BenchmarkRunner(BenchmarkConfig(input_folder=tmp, output_root=str(output_root)))
            job = {
                "pdb": "failed.pdb",
                "input_sha256": "input-hash",
                "chain_a": "A",
                "chain_b": "B",
                "structure_type": "experimental",
            }
            checkpoint = {
                "config_fingerprint": runner._checkpoint_fingerprint,
                "files": [
                    {
                        "pdb": "failed.pdb",
                        "input_sha256": "input-hash",
                        "interaction_sha256": None,
                        "chain_selection": {"chain_a": "A", "chain_b": "B"},
                        "structure_type": "experimental",
                        "status": "failed",
                        "error": "terminal attempt failed",
                    }
                ],
            }
            (output_root / runner.config.checkpoint_filename).write_text(
                json.dumps(checkpoint),
                encoding="utf-8",
            )

            completed, remaining = runner._load_resume_state([job])

        self.assertEqual(len(completed), 1)
        self.assertEqual(completed[0]["status"], "failed")
        self.assertEqual(remaining, [])

    def test_resume_rejects_duplicate_checkpoint_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_root = Path(tmp) / "out"
            output_root.mkdir()
            runner = BenchmarkRunner(BenchmarkConfig(input_folder=tmp, output_root=str(output_root)))
            record = {
                "pdb": "same.pdb",
                "input_sha256": "hash",
                "chain_selection": {"chain_a": "A", "chain_b": "B"},
                "structure_type": "experimental",
                "status": "failed",
            }
            (output_root / runner.config.checkpoint_filename).write_text(
                json.dumps(
                    {
                        "config_fingerprint": runner._checkpoint_fingerprint,
                        "files": [record, dict(record)],
                    }
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "duplicate structure records"):
                runner._load_resume_state([{"pdb": "same.pdb", "input_sha256": "hash", "chain_a": "A", "chain_b": "B"}])

    def test_resume_rejects_a_same_name_result_when_input_hash_changed(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_root = Path(tmp) / "out"
            output_root.mkdir()
            runner = BenchmarkRunner(BenchmarkConfig(input_folder=tmp, output_root=str(output_root)))
            checkpoint = {
                "config_fingerprint": runner._checkpoint_fingerprint,
                "files": [
                    {
                        "pdb": "same.pdb",
                        "input_sha256": "old-hash",
                        "chain_selection": {"chain_a": "A", "chain_b": "B"},
                        "status": "ok",
                    }
                ],
            }
            (output_root / runner.config.checkpoint_filename).write_text(
                json.dumps(checkpoint),
                encoding="utf-8",
            )
            jobs = [
                {
                    "pdb": "same.pdb",
                    "input_sha256": "new-hash",
                    "chain_a": "A",
                    "chain_b": "B",
                }
            ]

            completed, remaining = runner._load_resume_state(jobs)

        self.assertEqual(completed, [])
        self.assertEqual(remaining, jobs)


if __name__ == "__main__":
    unittest.main()
