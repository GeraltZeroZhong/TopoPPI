"""Auditable benchmark runner with frozen domains and isolated measurements."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import nullcontext
from dataclasses import asdict, replace
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy.spatial import cKDTree

try:
    from tqdm.auto import tqdm
except ImportError:  # optional dependency
    tqdm = None

try:
    import psutil
except ImportError:  # optional dependency
    psutil = None

try:
    import resource
except ImportError:  # unavailable on Windows
    resource = None

from topoppi import __version__
from topoppi.atlas.footprints import (
    atom_residue_label,
    contact_partner_degrees,
    geometric_contact_partner_map,
    residue_aware_residue_weights,
    residue_fragmentation_report,
    source_atom_residue_labels,
)
from topoppi.atlas.metrics import UVAtlasMetrics
from topoppi.atlas.packing import apply_packed_uv, pack_mesh_charts
from topoppi.atlas.uv import as_corner_uv, face_domain_hash, set_uv_layout
from topoppi.benchmark_methods import (
    PARAMETERIZATION_METHODS,
    RESIDUE_AWARE_BASELINE,
    RESIDUE_AWARE_OPTCUTS_METHODS,
    STANDARD_OPTCUTS_METHODS,
)
from topoppi.benchmarking.coordinate_audit import validate_coordinate_audit
from topoppi.benchmarking.manifest_metadata import (
    FORMAL_STRUCTURE_TYPES,
    INFERENCE_DEPENDENCY_BASIS,
    INFERENCE_DEPENDENCY_FIELDS,
    PREDICTED_STRUCTURE_TYPES,
    inference_family_id,
)
from topoppi.benchmarking.metrics_utils import (
    atlas_trainability_metrics,
    avg_energy,
    avg_seam_length,
    improvement_rate,
    quality_block,
    rasterize_feature_maps,
)
from topoppi.benchmarking.reporting import aggregate_results, write_csv
from topoppi.config import BenchmarkConfig
from topoppi.file_utils import git_worktree_state, sha256_file
from topoppi.interactions.interaction_engine import (
    load_prolif_partner_map,
    residue_sequence_token,
)
from topoppi.io.io_loader import PDBLoader
from topoppi.io.pdb_records import residue_plddt_values
from topoppi.json_utils import dump_json_atomic, json_safe
from topoppi.mesh.parameterization import Parameterizer
from topoppi.mesh.provenance import provenance_summary
from topoppi.mesh.surface import SurfaceGenerator
from topoppi.mesh.topology import TopologyManager
from topoppi.optimization.optcuts import (
    OptCutsUVOptimizer,
    resolve_optcuts_binary,
    supports_residue_footprint_energy,
)

STRUCTURE_SUFFIXES = (".pdb", ".cif", ".mmcif")
RESULT_IDENTITY_FIELDS = (
    "manifest_record_id",
    "cluster_id",
    "family_id",
    "sequence_cluster_a",
    "sequence_cluster_b",
    *INFERENCE_DEPENDENCY_FIELDS,
    "analysis_split",
    "analysis_split_component_id",
    "analysis_split_basis",
    "chain_a_residue_count",
    "chain_b_residue_count",
    "candidate_chain_pair_count",
    "selected_atom_contact_fraction",
    "selected_residue_contact_fraction",
    "structure_method",
    "resolution_angstrom",
    "experimental_methods_json",
    "experimental_method_group",
    "experimental_method_contains_nmr",
    "pdbbind_index_resolution_angstrom",
    "rcsb_resolution_combined_angstrom_json",
    "rcsb_experiment_metadata_source",
    "confidence_metric",
    "confidence_stratum",
    "afdb_ipsae_stratum",
    "paired_record_id",
    "paired_experimental_record_id",
    "paired_geometry_stratum",
    "paired_contact_cutoff_angstrom",
    "paired_predicted_contact_count_total",
    "paired_alignment_a_optimal_correspondence_count",
    "paired_alignment_b_optimal_correspondence_count",
    "paired_alignment_a_selected_pair_consensus_fraction",
    "paired_alignment_b_selected_pair_consensus_fraction",
)


def _result_identity_metadata(metadata: Dict[str, object]) -> Dict[str, object]:
    result = {
        field: (
            None
            if metadata.get(field) is None
            or (isinstance(metadata.get(field), str) and not str(metadata.get(field)).strip())
            else metadata.get(field)
        )
        for field in RESULT_IDENTITY_FIELDS
    }
    result["structure_type"] = metadata.get("structure_type") or "experimental"
    confidence_preflight = metadata.get("confidence_preflight")
    if isinstance(confidence_preflight, dict):
        result["confidence_preflight"] = confidence_preflight
    return result


def _uses_optcuts_uv(method: str) -> bool:
    return method in STANDARD_OPTCUTS_METHODS or method in RESIDUE_AWARE_OPTCUTS_METHODS


def _operational_method_censoring(worker_payload: Dict[str, object]) -> Dict[str, object] | None:
    """Describe an internal method-budget censoring event in an operational run."""

    result = worker_payload.get("result")
    if not isinstance(result, dict):
        return None
    if str(result.get("execution_profile") or "").strip().lower() != "operational_optcuts":
        return None
    if str(result.get("status") or "").strip().lower() == "ok":
        return None
    method = str(result.get("operational_method") or "").strip()
    executions = result.get("method_execution")
    diagnostics = executions.get(method) if isinstance(executions, dict) else None
    if not isinstance(diagnostics, dict):
        return None
    failures = diagnostics.get("failures")
    failure_types = {
        str(failure.get("failure_type") or "").strip().lower()
        for failure in failures or []
        if isinstance(failure, dict)
    }
    if not failure_types.intersection({"timeout", "arm_budget_exhausted"}):
        return None
    reason = "optcuts_method_timeout" if "timeout" in failure_types else "optcuts_method_budget_exhausted"
    timing = result.get("timing")
    end_to_end = timing.get("end_to_end") if isinstance(timing, dict) else None
    return {
        "termination_reason": reason,
        "censoring_threshold_sec": diagnostics.get("method_arm_time_budget_sec"),
        "censoring_event_elapsed_sec": (end_to_end.get("wall_sec") if isinstance(end_to_end, dict) else None),
    }


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _first_present(*values: object, default: object = None) -> object:
    for value in values:
        if value is None or (isinstance(value, str) and not value.strip()):
            continue
        return value
    return default


def _afdb_complex_confidence(metadata: Dict[str, object]) -> Dict[str, object]:
    """Expose complex-level scores only for an actual predicted complex."""

    if str(metadata.get("structure_type") or "").strip().lower() != "afdb":
        return {}
    return {
        "model_id": _first_present(metadata.get("afdb_model_id")),
        "iptm": _first_present(metadata.get("afdb_iptm")),
        "ipsae": _first_present(metadata.get("afdb_ipsae")),
        "pdockq": _first_present(metadata.get("afdb_pdockq")),
        "pdockq2": _first_present(metadata.get("afdb_pdockq2")),
        "lis": _first_present(metadata.get("afdb_lis")),
    }


def _optional_finite_float(value: object) -> float | None:
    """Normalize an optional manifest scalar without inventing a value."""

    if value is None or (isinstance(value, str) and not value.strip()):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def _optional_nonnegative_int(value: object) -> int | None:
    if value is None or (isinstance(value, str) and not value.strip()):
        return None
    try:
        number = int(value)
    except (TypeError, ValueError):
        return None
    return number if number >= 0 else None


def _worker_cpu_time() -> float:
    """Return CPU seconds for this worker and child processes it has reaped."""

    total = float(time.process_time())
    if resource is not None:
        children = resource.getrusage(resource.RUSAGE_CHILDREN)
        total += float(children.ru_utime + children.ru_stime)
    return total


class BenchmarkRunner:
    def __init__(
        self,
        config: BenchmarkConfig,
        log_fn: Optional[Callable[[str], None]] = None,
        progress_fn: Optional[Callable[[int, int, str], None]] = None,
        cancel_event=None,
        worker_mode: bool = False,
    ):
        self.config = config
        self.log = log_fn or (lambda _message: None)
        self.progress = progress_fn or (lambda *_args: None)
        self.cancel_event = cancel_event
        self._proc = psutil.Process(os.getpid()) if psutil else None
        self.config.validate()
        self._checkpoint_path = os.path.join(self.config.output_root, self.config.checkpoint_filename)
        self._resolved_optcuts_binary = None if worker_mode else resolve_optcuts_binary(self.config.optcuts)
        self._resolved_optcuts_sha256 = (
            sha256_file(self._resolved_optcuts_binary) if self._resolved_optcuts_binary else None
        )
        # Isolated workers never read or write the parent checkpoint.  Avoid
        # hashing the binary and invoking Git inside every timed subprocess.
        self._checkpoint_fingerprint = "worker-not-applicable" if worker_mode else self._config_fingerprint()
        self._preflight_context: Dict[str, object] | None = None
        self._preflight_jobs: List[Dict[str, object]] | None = None
        self._preflight_preprocessing: Dict[str, object] | None = None
        self._available_worker_cpus = self._ordered_available_cpus()
        self._affinity_lock = threading.Lock()
        self._worker_affinity_by_thread: Dict[int, List[int]] = {}
        self._next_affinity_slot = 0
        self._validated_detail_artifacts: set[str] = set()

    def _configured_methods(self) -> tuple[str, ...]:
        """Return only methods actually executed by the selected profile."""

        if self.config.execution_profile.strip().lower() == "operational_optcuts":
            return self.config.resolved_optcuts_variants()
        return (*PARAMETERIZATION_METHODS, *self.config.resolved_optcuts_variants())

    def _check_cancelled(self) -> None:
        if self.cancel_event is not None and self.cancel_event.is_set():
            raise RuntimeError("Benchmark cancelled by user.")

    def run(self) -> Dict[str, object]:
        self._check_cancelled()
        preflight = self.preflight()
        prepared_jobs = list(self._preflight_jobs)
        preprocessing_log = dict(self._preflight_preprocessing)
        if not preflight["ready"]:
            if not preflight["structure_file_count"]:
                raise ValueError("No PDB/mmCIF files found for benchmark.")
            if not self.config.formal_mode and not prepared_jobs and bool(preflight["output_state"]["acceptable"]):
                os.makedirs(self.config.output_root, exist_ok=True)
                output = self._build_output(
                    preprocessing_log=preprocessing_log,
                    files=[],
                    worker_count=0,
                )
                self._write_outputs(output, [])
                raise ValueError("No valid structures after preprocessing. See benchmark_report.json for reasons.")
            blockers = "; ".join(str(item) for item in preflight["blockers"])
            label = "Formal benchmark" if self.config.formal_mode else "Benchmark"
            raise ValueError(f"{label} preflight failed; no jobs were started: {blockers}")

        os.makedirs(self.config.output_root, exist_ok=True)
        self._check_cancelled()
        completed_results, prepared_jobs = self._load_resume_state(prepared_jobs)
        if self.config.resume and not os.path.isfile(self._checkpoint_path):
            if not self._save_checkpoint(completed_results):
                raise OSError("Failed to initialize the benchmark checkpoint.")

        worker_count = self._resolve_worker_count(len(prepared_jobs))
        self.log(
            f"[Benchmark] Supervising {worker_count} isolated process(es); "
            f"measured repetitions={self.config.repetitions}, warmups={self.config.warmup_runs}."
        )
        total_jobs = len(completed_results) + len(prepared_jobs)
        self._safe_progress(len(completed_results), total_jobs, "Benchmark started")
        new_results = (
            self._run_files_concurrently(
                prepared_jobs,
                worker_count,
                completed_results=completed_results,
                total_jobs=total_jobs,
            )
            if prepared_jobs
            else []
        )
        self._check_cancelled()
        all_results = completed_results + new_results
        output = self._build_output(preprocessing_log=preprocessing_log, files=all_results, worker_count=worker_count)
        self._write_outputs(output, all_results)
        return output

    def preflight(self) -> Dict[str, object]:
        """Validate inputs and frozen artifacts without starting benchmark jobs."""

        if self._preflight_context is None:
            structure_files = sorted(
                filename
                for filename in os.listdir(self.config.input_folder)
                if filename.lower().endswith(STRUCTURE_SUFFIXES)
            )
            try:
                jobs, preprocessing = self._prepare_benchmark_jobs(structure_files)
            except (OSError, ValueError, csv.Error, json.JSONDecodeError) as exc:
                self._check_cancelled()
                jobs = []
                preprocessing = {
                    "total_files": int(len(structure_files)),
                    "accepted_files": 0,
                    "skipped_files": 1,
                    "integrity_error_count": 1,
                    "accepted": [],
                    "skipped": [
                        {
                            "pdb": os.path.basename(self.config.manifest_path) or "benchmark_manifest",
                            "reason": f"Manifest/preprocessing validation failed: {exc}",
                            "fatal_integrity_error": True,
                        }
                    ],
                    "rules": ["Preflight exceptions are reported as fatal integrity errors."],
                }
            resolved = self._resolved_optcuts_binary
            binary_sha256 = self._resolved_optcuts_sha256
            expected_sha256 = self.config.optcuts.expected_binary_sha256.strip().lower() or None
            binary_matches = bool(
                resolved and (expected_sha256 is None or str(binary_sha256).lower() == expected_sha256)
            )
            residue_aware_capable = (
                bool(resolved and supports_residue_footprint_energy(resolved))
                if self.config.optcuts.residue_fragmentation_weight > 0.0
                else None
            )
            self._preflight_context = {
                "structure_files": structure_files,
                "jobs": jobs,
                "preprocessing": preprocessing,
                "resolved": resolved,
                "binary_sha256": binary_sha256,
                "expected_sha256": expected_sha256,
                "binary_matches": binary_matches,
                "residue_aware_capable": residue_aware_capable,
                "environment": self._environment_metadata(),
                "coordinate_audit": self._coordinate_audit_preflight(),
            }
            self._preflight_jobs = jobs
            self._preflight_preprocessing = preprocessing

        context = self._preflight_context
        structure_files = context["structure_files"]
        jobs = context["jobs"]
        preprocessing = context["preprocessing"]
        resolved = context["resolved"]
        binary_sha256 = context["binary_sha256"]
        expected_sha256 = context["expected_sha256"]
        binary_matches = context["binary_matches"]
        residue_aware_capable = context["residue_aware_capable"]
        environment = context["environment"]
        coordinate_audit = context["coordinate_audit"]
        output_state = self._output_state()
        remaining_jobs = jobs
        resumed_structure_count = 0
        if output_state.get("state") == "matching_resume_checkpoint":
            try:
                resumed, remaining_jobs = self._load_resume_state(jobs)
                resumed_structure_count = len(resumed)
            except ValueError as exc:
                output_state = {
                    "acceptable": False,
                    "state": "invalid_resume_checkpoint",
                    "reason": str(exc),
                }
        blockers = self._preflight_blockers(
            structure_files,
            jobs,
            preprocessing,
            resolved=resolved,
            binary_matches=binary_matches,
            residue_aware_capable=residue_aware_capable,
            environment=environment,
            output_state=output_state,
            coordinate_audit=coordinate_audit,
        )
        result = {
            "ready": not blockers,
            "formal_mode": bool(self.config.formal_mode),
            "benchmark_purpose": self.config.benchmark_purpose.strip().lower(),
            "execution_profile": self.config.execution_profile.strip().lower(),
            "optcuts_variants": list(self.config.resolved_optcuts_variants()),
            "structure_file_count": int(len(structure_files)),
            "accepted_job_count": int(len(jobs)),
            "resumed_structure_count": int(resumed_structure_count),
            "remaining_job_count": int(len(remaining_jobs)),
            "planned_worker_process_count": int(
                len(remaining_jobs) * (self.config.warmup_runs + self.config.repetitions)
            ),
            "preprocessing": preprocessing,
            "optcuts": {
                "resolved_path": resolved,
                "actual_sha256": binary_sha256,
                "expected_sha256": expected_sha256,
                "checksum_matches": binary_matches,
                "residue_aware_capable": residue_aware_capable,
            },
            "config_fingerprint": self._checkpoint_fingerprint,
            "environment": environment,
            "output_state": output_state,
            "coordinate_audit": coordinate_audit,
            "blockers": blockers,
            "note": "Preflight is read-only and does not execute surface generation or OptCuts.",
        }
        return result

    def _preflight_blockers(
        self,
        structure_files,
        jobs,
        preprocessing,
        *,
        resolved,
        binary_matches,
        residue_aware_capable,
        environment,
        output_state,
        coordinate_audit,
    ) -> List[str]:
        blockers = []
        if not structure_files:
            blockers.append("No PDB/mmCIF inputs were found.")
        if not jobs:
            blockers.append("No input passed preprocessing.")
        if int(preprocessing.get("integrity_error_count", 0)):
            blockers.append("Dataset manifest has fatal integrity errors.")
        if not resolved:
            blockers.append("OptCuts binary could not be resolved.")
        elif not binary_matches:
            blockers.append("OptCuts binary checksum does not match the frozen configuration.")
        if self.config.optcuts.residue_fragmentation_weight > 0.0 and not residue_aware_capable:
            blockers.append("OptCuts binary lacks residue-footprint energy support.")
        if self.config.formal_mode and psutil is None:
            blockers.append("psutil is required for formal peak-RSS measurement.")
        if self.config.formal_mode and not environment.get("git_commit"):
            blockers.append("Formal benchmark requires a resolvable committed Git revision.")
        if self.config.formal_mode and environment.get("git_worktree_dirty") is not False:
            blockers.append("Formal benchmark requires a clean Git worktree with known status.")
        expected_commit = str(self.config.expected_git_commit).strip().lower()
        actual_commit = str(environment.get("git_commit") or "").strip().lower()
        if self.config.formal_mode and expected_commit and actual_commit != expected_commit:
            blockers.append(
                "Git revision does not match benchmark.expected_git_commit "
                f"(expected {expected_commit}, got {actual_commit or 'unavailable'})."
            )
        if not output_state["acceptable"]:
            blockers.append(str(output_state["reason"]))
        if self.config.formal_mode and coordinate_audit.get("status") != "validated":
            blockers.append(
                "Formal benchmark requires a passing coordinate audit bound to the current manifest: "
                + str(coordinate_audit.get("reason") or coordinate_audit.get("status") or "unknown error")
            )
        return blockers

    def _coordinate_audit_preflight(self) -> Dict[str, object]:
        """Validate the frozen coordinate audit without recomputing it."""

        if not self.config.formal_mode:
            return {"status": "not_required"}
        return validate_coordinate_audit(
            self.config.coordinate_audit_path,
            self.config.expected_coordinate_audit_sha256,
            self.config.manifest_path,
        )

    def _output_state(self) -> Dict[str, object]:
        root = Path(self.config.output_root)
        if not root.exists():
            return {"acceptable": True, "state": "absent"}
        if not root.is_dir():
            return {
                "acceptable": False,
                "state": "not_a_directory",
                "reason": "Benchmark output_root exists but is not a directory.",
            }
        entries = list(root.iterdir())
        if not entries:
            return {"acceptable": True, "state": "empty"}
        checkpoint = root / self.config.checkpoint_filename
        if self.config.resume and checkpoint.is_file():
            try:
                with checkpoint.open("r", encoding="utf-8") as handle:
                    payload = json.load(handle)
                if str(payload.get("config_fingerprint") or "") == self._checkpoint_fingerprint:
                    return {
                        "acceptable": True,
                        "state": "matching_resume_checkpoint",
                        "checkpoint": str(checkpoint),
                    }
            except (OSError, json.JSONDecodeError):
                pass
        return {
            "acceptable": False,
            "state": "nonempty_unmatched",
            "reason": "Benchmark output_root must be empty or contain a matching resume checkpoint.",
        }

    def _build_output(
        self,
        preprocessing_log: Dict[str, object],
        files: List[Dict[str, object]],
        worker_count: int,
    ) -> Dict[str, object]:
        return {
            "schema_version": "2.0",
            "created_at": _utc_now(),
            "topoppi_version": __version__,
            "config": asdict(self.config),
            "runtime": {
                "worker_count": int(worker_count),
                "execution_model": "fresh_subprocess_per_structure_repetition",
                "formal_mode": bool(self.config.formal_mode),
                "config_fingerprint": self._checkpoint_fingerprint,
                "python": sys.version.split()[0],
                "platform": platform.platform(),
                "cpu_count": os.cpu_count(),
                "threads_per_worker": int(self.config.threads_per_worker),
                "optcuts_method_arm_time_budget_sec": float(self.config.optcuts.timeout_sec),
                "optcuts_method_arm_time_budget_scope": (
                    "shared across all external solver invocations for one method and structure"
                ),
                "environment": self._environment_metadata(),
                "coordinate_audit": self._coordinate_audit_preflight(),
            },
            "metric_protocol": {
                "domain": (
                    "not_applicable_operational_runtime_only"
                    if self.config.execution_profile.strip().lower() == "operational_optcuts"
                    else "exact_common_source_faces_across_all_methods"
                ),
                "methods": list(self._configured_methods()),
                "uv_representation": "per_face_corner",
                "shape_and_area_scale_alignment": "global_total_area_similarity_per_patch",
                "symmetric_dirichlet_scale_alignment": (
                    "analytic_global_minimum_per_patch_after_similarity-normalized_jacobian_construction"
                ),
                "face_weighting": "original_3d_area",
                "reflection": "one_global_orientation_correction_per_patch",
                "atlas_overlap": ("polygonal_triangle_union_in_normalized_double_precision_with_reported_tolerance"),
                "seam": "discontinuous_internal_3d_edges",
                "residue_footprint_fragmentation": ("mass_aware_fragmentation_within_original_face_dual_components"),
                "detailed_evidence": {
                    "per_face_sample": self.config.per_face_sample_filename,
                    "per_residue": self.config.per_residue_filename,
                    "provenance": self.config.provenance_filename,
                    "optcuts_executions": self.config.optcuts_execution_filename,
                },
            },
            "preprocessing": preprocessing_log,
            "files": files,
            "summary": aggregate_results(
                files,
                methods=self._configured_methods(),
                bootstrap_iterations=self.config.bootstrap_iterations,
                random_seed=self.config.random_seed,
            ),
        }

    def _write_outputs(self, output: Dict[str, object], results: List[Dict[str, object]]) -> None:
        self._write_provenance_csv(results)
        self._write_per_face_sample_csv(results)
        self._write_per_residue_csv(results)
        self._write_optcuts_execution_jsonl(results)
        for result in results:
            result.pop("provenance_records", None)
            result.pop("per_face_sample_records", None)
            result.pop("per_residue_records", None)
            result["provenance_artifact"] = self.config.provenance_filename
            result["per_face_sample_artifact"] = self.config.per_face_sample_filename
            result["per_residue_artifact"] = self.config.per_residue_filename
            result["optcuts_execution_artifact"] = self.config.optcuts_execution_filename
        report_path = os.path.join(self.config.output_root, self.config.report_filename)
        dump_json_atomic(output, report_path)
        if not self._save_checkpoint(results):
            raise OSError("Failed to write the final benchmark checkpoint.")
        write_csv(results, self.config.output_root, filename=self.config.summary_filename)
        self._write_manifest_csv(output["preprocessing"])
        self._write_failure_csv(output["preprocessing"], results)
        self._write_per_patch_csv(results)
        self._write_artifact_checksums()

    def _write_provenance_csv(self, results: List[Dict[str, object]]) -> None:
        path = os.path.join(self.config.output_root, self.config.provenance_filename)
        temporary_path = path + ".tmp"
        fields = [
            "pdb",
            "patch_id",
            "entity",
            "final_index",
            "source_id",
            "source_atom_index",
        ]
        try:
            with gzip.open(temporary_path, "wt", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
                writer.writeheader()
                for result in results:
                    provenance_records = self._detail_records(result, "provenance_records")
                    for record in provenance_records:
                        writer.writerow({"pdb": result.get("pdb"), **record})
            os.replace(temporary_path, path)
        except Exception:
            Path(temporary_path).unlink(missing_ok=True)
            raise

    def _write_manifest_csv(self, preprocessing: Dict[str, object]) -> None:
        path = os.path.join(self.config.output_root, self.config.manifest_filename)
        fields = [
            "pdb",
            "input_sha256",
            "chain_a",
            "chain_b",
            "selection_mode",
            "cluster_id",
            "family_id",
            "sequence_cluster_a",
            "sequence_cluster_b",
            *INFERENCE_DEPENDENCY_FIELDS,
            "analysis_split",
            "analysis_split_component_id",
            "analysis_split_basis",
            "chain_a_residue_count",
            "chain_b_residue_count",
            "candidate_chain_pair_count",
            "selected_atom_contact_fraction",
            "selected_residue_contact_fraction",
            "dataset_source",
            "source_accession",
            "license_or_terms",
            "structure_type",
            "structure_method",
            "resolution_angstrom",
            "experimental_methods_json",
            "experimental_method_group",
            "experimental_method_contains_nmr",
            "pdbbind_index_resolution_angstrom",
            "rcsb_resolution_combined_angstrom_json",
            "rcsb_experiment_metadata_source",
            "confidence_metric",
            "confidence_stratum",
            "afdb_ipsae_stratum",
            "confidence_source",
            "confidence_threshold",
            "afdb_model_id",
            "afdb_iptm",
            "afdb_ipsae",
            "afdb_pdockq",
            "afdb_pdockq2",
            "afdb_lis",
            "paired_record_id",
            "paired_experimental_record_id",
            "paired_geometry_stratum",
            "paired_contact_cutoff_angstrom",
            "paired_predicted_contact_count_total",
            "paired_contact_recall_fnat",
            "paired_contact_precision",
            "paired_contact_jaccard",
            "paired_experimental_contact_mapping_coverage",
            "paired_interface_residue_a_mapping_coverage",
            "paired_interface_residue_b_mapping_coverage",
            "paired_interface_ligand_ca_mapping_coverage",
            "paired_interface_ligand_ca_rmsd_angstrom",
            "paired_cross_chain_clash_atom_fraction",
            "paired_alignment_a_optimal_correspondence_count",
            "paired_alignment_b_optimal_correspondence_count",
            "paired_alignment_a_selected_pair_consensus_fraction",
            "paired_alignment_b_selected_pair_consensus_fraction",
            "confidence_preflight_status",
            "confidence_preflight_summary_unit",
            "confidence_preflight_atom_count",
            "confidence_preflight_residue_count",
            "confidence_preflight_minimum",
            "confidence_preflight_mean",
            "confidence_preflight_maximum",
            "surface_estimate_status",
            "surface_requested_voxel_count",
            "surface_effective_voxel_count",
            "surface_effective_resolution_angstrom",
            "surface_dense_field_bytes_lower_bound",
            "manifest_record_id",
            "hotspot_residues_a",
            "prolif_file",
            "prolif_sha256",
            "available_chains",
            "status",
            "reason",
            "fatal_integrity_error",
        ]
        accepted = [{**item, "status": "accepted", "reason": ""} for item in preprocessing.get("accepted", [])]
        skipped = [{**item, "status": "skipped"} for item in preprocessing.get("skipped", [])]
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            for row in accepted + skipped:
                row = dict(row)
                row["available_chains"] = json.dumps(row.get("available_chains", []), separators=(",", ":"))
                confidence = row.get("confidence_preflight", {})
                if isinstance(confidence, dict):
                    row["confidence_preflight_status"] = confidence.get("status", "")
                    row["confidence_preflight_summary_unit"] = confidence.get("summary_unit", "")
                    row["confidence_preflight_atom_count"] = confidence.get("atom_count", "")
                    row["confidence_preflight_residue_count"] = confidence.get("residue_count", "")
                    row["confidence_preflight_minimum"] = confidence.get("minimum", "")
                    row["confidence_preflight_mean"] = confidence.get("mean", "")
                    row["confidence_preflight_maximum"] = confidence.get("maximum", "")
                writer.writerow(row)

    def _write_failure_csv(self, preprocessing: Dict[str, object], results: List[Dict[str, object]]) -> None:
        path = os.path.join(self.config.output_root, self.config.failures_filename)
        fields = ["pdb", "stage", "method", "patch_id", "reason"]
        rows = [
            {
                "pdb": item.get("pdb"),
                "stage": "preprocessing",
                "method": "",
                "patch_id": "",
                "reason": item.get("reason"),
            }
            for item in preprocessing.get("skipped", [])
        ]
        for result in results:
            rows.extend(self._result_failure_rows(result))
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    @staticmethod
    def _result_failure_rows(result: Dict[str, object]) -> List[Dict[str, object]]:
        pdb = result.get("pdb")
        rows = []
        if error := result.get("error"):
            rows.append({"pdb": pdb, "stage": "structure", "method": "", "patch_id": "", "reason": error})
        for method, block in result.get("method_execution", {}).items():
            if isinstance(block, dict):
                rows.extend(
                    {
                        "pdb": pdb,
                        "stage": "method",
                        "method": method,
                        "patch_id": failure.get("patch_id", ""),
                        "reason": failure.get("reason", "unknown"),
                    }
                    for failure in block.get("failures", [])
                )
        preparation = result.get("topology_preparation", {})
        if isinstance(preparation, dict):
            rows.extend(
                {
                    "pdb": pdb,
                    "stage": "topology_preparation",
                    "method": "shared_domain_gate",
                    "patch_id": patch.get("patch_id", ""),
                    "reason": patch["failure_reason"],
                }
                for patch in preparation.get("patches", [])
                if patch.get("failure_reason")
            )
        extraction = result.get("topology_extraction", {})
        if isinstance(extraction, dict):
            rows.extend(
                {
                    "pdb": pdb,
                    "stage": "topology_extraction",
                    "method": "",
                    "patch_id": component.get("patch_id", ""),
                    "reason": component.get("reason", "component_dropped"),
                }
                for component in extraction.get("components", [])
                if component.get("status") != "accepted"
            )
        return rows

    def _write_per_patch_csv(self, results: List[Dict[str, object]]) -> None:
        path = os.path.join(self.config.output_root, self.config.per_patch_filename)
        fields = [
            "pdb",
            "patch_id",
            "source_face_hash",
            "retention_status",
            "rejection_stage",
            "failure_reason",
            "face_count_before",
            "face_count_after_topology_sanitation",
            "face_count_after",
            "materialized_vertex_count_before",
            "materialized_vertex_count_after_topology_sanitation",
            "materialized_vertex_count_after",
            "source_vertex_count_before",
            "source_vertex_count_after_topology_sanitation",
            "source_vertex_count_after",
            "area_before",
            "area_after_topology_sanitation",
            "area_after",
            "source_atom_count_before",
            "source_atom_count_after_topology_sanitation",
            "source_atom_count_after",
            "face_retention_ratio",
            "materialized_vertex_count_ratio",
            "source_vertex_retention_ratio",
            "area_retention_ratio",
            "source_atom_retention_ratio",
            "topology_face_retention_ratio",
            "topology_materialized_vertex_count_ratio",
            "topology_source_vertex_retention_ratio",
            "topology_area_retention_ratio",
            "topology_source_atom_retention_ratio",
            "parameterization_face_retention_ratio",
            "parameterization_materialized_vertex_count_ratio",
            "parameterization_source_vertex_retention_ratio",
            "parameterization_area_retention_ratio",
            "parameterization_source_atom_retention_ratio",
            "residue_retention_ratio",
            "topology_residue_retention_ratio",
            "parameterization_residue_retention_ratio",
            "geometric_contact_pair_retention_ratio",
            "topology_geometric_contact_pair_retention_ratio",
            "parameterization_geometric_contact_pair_retention_ratio",
            "declared_hotspot_retention_ratio",
            "topology_declared_hotspot_retention_ratio",
            "parameterization_declared_hotspot_retention_ratio",
            "declared_interaction_retention_ratio",
            "topology_declared_interaction_retention_ratio",
            "parameterization_declared_interaction_retention_ratio",
            "confidence_metric",
            "confidence_threshold",
            "confidence_atom_count_before",
            "confidence_atom_count_after_topology_sanitation",
            "confidence_atom_count_after",
            "confidence_atom_retention_ratio",
            "topology_confidence_atom_retention_ratio",
            "parameterization_confidence_atom_retention_ratio",
            "residue_count_before",
            "residue_count_after_topology_sanitation",
            "residue_count_after",
            "geometric_contact_pair_count_before",
            "geometric_contact_pair_count_after_topology_sanitation",
            "geometric_contact_pair_count_after",
            "declared_hotspot_count_on_patch_before",
            "declared_hotspot_count_after_topology_sanitation",
            "declared_hotspot_count_after",
            "declared_interaction_count_on_patch_before",
            "declared_interaction_count_after_topology_sanitation",
            "declared_interaction_count_after",
            "confidence_mean_before",
            "confidence_mean_after",
            "low_confidence_atom_fraction_before",
            "low_confidence_atom_fraction_after",
            "retention_denominator",
        ]
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            for result in results:
                for record in result.get("patch_records", []):
                    writer.writerow({"pdb": result.get("pdb"), **record})

    def _write_per_face_sample_csv(self, results: List[Dict[str, object]]) -> None:
        path = os.path.join(self.config.output_root, self.config.per_face_sample_filename)
        fields = [
            "pdb",
            "patch_id",
            "source_face_id",
            "face_area_3d",
            "sampling_rank",
            "sampling_rule",
        ]
        for method in self._configured_methods():
            fields.extend(
                [
                    f"{method}_distortion",
                    f"{method}_symmetric_dirichlet",
                    f"{method}_angle_distortion_rad",
                    f"{method}_area_distortion",
                    f"{method}_flipped_after_global_reflection",
                ]
            )
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            for result in results:
                for record in self._detail_records(result, "per_face_sample_records"):
                    writer.writerow({"pdb": result.get("pdb"), **record})

    def _write_per_residue_csv(self, results: List[Dict[str, object]]) -> None:
        path = os.path.join(self.config.output_root, self.config.per_residue_filename)
        fields = [
            "pdb",
            "evidence_domain",
            "domain_signature",
            "method",
            "residue",
            "face_count",
            "dual_edge_count",
            "baseline_component_count",
            "component_count_after_seams",
            "extra_component_count",
            "cycle_rank",
            "seam_crossing_edge_count",
            "nonseparating_seam_crossing_edge_count",
            "seam_crossing_length_3d",
            "footprint_area",
            "fragmentation",
            "interaction_weight",
            "objective_weight",
        ]
        with gzip.open(path, "wt", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            for result in results:
                for record in self._detail_records(result, "per_residue_records"):
                    writer.writerow({"pdb": result.get("pdb"), **record})

    def _write_optcuts_execution_jsonl(self, results: List[Dict[str, object]]) -> None:
        path = os.path.join(self.config.output_root, self.config.optcuts_execution_filename)
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            for result in results:
                detail = self._load_detail_payload(result)
                methods = detail.get("optcuts_execution_details", {})
                topology_ablation = detail.get("topology_ablation_execution")
                if not methods and topology_ablation is None:
                    methods = {
                        method: block.get("executions", [])
                        for method, block in result.get("method_execution", {}).items()
                        if isinstance(block, dict) and block.get("executions")
                    }
                record = {
                    "pdb": result.get("pdb"),
                    "methods": methods,
                    "topology_ablation": topology_ablation,
                }
                handle.write(json.dumps(json_safe(record), sort_keys=True, separators=(",", ":")))
                handle.write("\n")

    def _detail_records(self, result: Dict[str, object], key: str) -> List[Dict[str, object]]:
        direct = result.get(key)
        if isinstance(direct, list):
            return direct
        records = self._load_detail_payload(result).get(key, [])
        return records if isinstance(records, list) else []

    def _load_detail_payload(self, result: Dict[str, object]) -> Dict[str, object]:
        artifact = result.get("detail_artifact")
        if not isinstance(artifact, dict):
            return {}
        path = str(artifact.get("path") or "")
        expected = str(artifact.get("sha256") or "").lower()
        if path and not os.path.isabs(path):
            path = os.path.join(self.config.output_root, path)
        if not path or not os.path.isfile(path):
            raise OSError(f"Detailed benchmark evidence is missing for {result.get('pdb')}: {path}")
        if path not in self._validated_detail_artifacts:
            if sha256_file(path).lower() != expected:
                raise OSError(f"Detailed benchmark evidence checksum mismatch for {result.get('pdb')}.")
            self._validated_detail_artifacts.add(path)
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            payload = json.load(handle)
        if not isinstance(payload, dict) or payload.get("schema_version") != "1.0":
            raise ValueError(f"Detailed benchmark evidence has an invalid schema for {result.get('pdb')}.")
        return payload

    def _write_artifact_checksums(self) -> None:
        checksum_name = self.config.artifact_checksums_filename
        names = (
            self.config.report_filename,
            self.config.summary_filename,
            self.config.checkpoint_filename,
            self.config.manifest_filename,
            self.config.failures_filename,
            self.config.per_patch_filename,
            self.config.per_face_sample_filename,
            self.config.per_residue_filename,
            self.config.provenance_filename,
            self.config.optcuts_execution_filename,
        )
        artifacts = []
        missing = []
        for name in names:
            path = os.path.join(self.config.output_root, name)
            if os.path.isfile(path):
                artifacts.append(
                    {
                        "filename": name,
                        "bytes": int(os.path.getsize(path)),
                        "sha256": sha256_file(path),
                    }
                )
            else:
                missing.append(name)
        if missing:
            raise OSError("Benchmark evidence bundle is incomplete; missing: " + ", ".join(missing))
        path = os.path.join(self.config.output_root, checksum_name)
        dump_json_atomic(
            {
                "created_at": _utc_now(),
                "algorithm": "sha256",
                "config_fingerprint": self._checkpoint_fingerprint,
                "artifacts": artifacts,
                "note": "Worker logs are intentionally excluded because they are numerous and run-specific.",
            },
            path,
        )

    def _config_fingerprint(self) -> str:
        payload = asdict(self.config)
        payload.pop("output_root", None)
        payload["topoppi_version"] = __version__
        repo_root = Path(__file__).resolve().parents[3]
        payload["git_commit"], _git_dirty = git_worktree_state(repo_root)
        payload["resolved_optcuts_binary_sha256"] = self._resolved_optcuts_sha256
        if self.config.manifest_path and os.path.isfile(self.config.manifest_path):
            payload["manifest_sha256"] = sha256_file(self.config.manifest_path)
        data = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
        return hashlib.sha256(data).hexdigest()

    def _run_files_concurrently(
        self,
        jobs: List[Dict[str, object]],
        worker_count: int,
        completed_results: Optional[List[Dict[str, object]]] = None,
        total_jobs: Optional[int] = None,
    ) -> List[Dict[str, object]]:
        completed_results = completed_results or []
        total_jobs = int(total_jobs if total_jobs is not None else len(jobs))
        all_results: List[Optional[Dict[str, object]]] = [None] * len(jobs)
        completed = len(completed_results)
        progress_ctx = (
            tqdm(
                total=total_jobs,
                desc="Benchmark",
                unit="file",
                disable=not bool(self.config.show_tqdm),
                initial=completed,
            )
            if tqdm is not None
            else nullcontext(None)
        )
        with self._affinity_lock:
            self._worker_affinity_by_thread.clear()
            self._next_affinity_slot = 0
        completions_since_checkpoint = 0
        with ThreadPoolExecutor(max_workers=worker_count) as executor, progress_ctx as progress_bar:
            future_to_job = {}
            for index, job in enumerate(jobs):
                self._check_cancelled()
                self.log(
                    f"[Benchmark] Queued {job['pdb']} with chains {job['chain_a']}/{job['chain_b']} in a fresh process."
                )
                future_to_job[executor.submit(self._run_isolated_job, job)] = (index, job)

            for future in as_completed(future_to_job):
                self._check_cancelled()
                output_index, failed_job = future_to_job[future]
                pdb_name = str(failed_job["pdb"])
                completed += 1
                try:
                    all_results[output_index] = future.result()
                    self.log(f"[Benchmark] Finished {pdb_name}")
                except Exception as exc:
                    self._check_cancelled()
                    all_results[output_index] = {
                        "pdb": pdb_name,
                        "input_sha256": failed_job.get("input_sha256") or None,
                        "interaction_sha256": failed_job.get("prolif_sha256") or None,
                        "chain_selection": {
                            "chain_a": failed_job.get("chain_a"),
                            "chain_b": failed_job.get("chain_b"),
                            "mode": failed_job.get("selection_mode", "configured"),
                            "details": failed_job.get("selection_details", {}),
                        },
                        "patch_count": 0,
                        "status": "failed",
                        "error": f"Isolated benchmark execution failed: {exc}",
                        "benchmark_purpose": self.config.benchmark_purpose,
                        "execution_profile": self.config.execution_profile,
                        "topology_ablation_configured": bool(self.config.include_topology_ablation),
                        **_result_identity_metadata(failed_job),
                    }
                    self.log(f"[Benchmark] Failed {pdb_name}: {exc}")
                self._safe_progress(completed, total_jobs, f"Finished {pdb_name}")
                if progress_bar is not None:
                    progress_bar.update(1)
                    progress_bar.set_postfix_str(pdb_name)
                completions_since_checkpoint += 1
                checkpoint_due = bool(
                    completions_since_checkpoint >= int(self.config.checkpoint_interval_structures)
                    or completed == total_jobs
                )
                if checkpoint_due:
                    checkpoint = completed_results + [item for item in all_results if item is not None]
                    if self._save_checkpoint(checkpoint):
                        completions_since_checkpoint = 0
        return [item for item in all_results if item is not None]

    def _run_isolated_job(self, job: Dict[str, object]) -> Dict[str, object]:
        measured = []
        all_runs = []
        total_runs = int(self.config.warmup_runs + self.config.repetitions)
        for run_index in range(total_runs):
            self._check_cancelled()
            is_warmup = run_index < int(self.config.warmup_runs)
            outcome = self._execute_worker(job, run_index=run_index, is_warmup=is_warmup)
            all_runs.append(outcome["measurement"])
            if outcome["payload"].get("status") != "ok":
                measurement = outcome["measurement"]
                measured_attempt = not is_warmup
                runtime_observation = measurement.get(
                    "runtime_observation_sec",
                    measurement["wall_sec"],
                )
                return {
                    "pdb": job["pdb"],
                    "input_sha256": job.get("input_sha256") or None,
                    "interaction_sha256": job.get("prolif_sha256") or None,
                    "chain_selection": {
                        "chain_a": job["chain_a"],
                        "chain_b": job["chain_b"],
                        "mode": job.get("selection_mode", self.config.chain_selection_mode),
                        "details": job.get("selection_details", {}),
                    },
                    "patch_count": 0,
                    "status": "failed",
                    "error": outcome["payload"].get("error", "Worker failed"),
                    "benchmark_purpose": self.config.benchmark_purpose,
                    "execution_profile": self.config.execution_profile,
                    "topology_ablation_configured": bool(self.config.include_topology_ablation),
                    "timing": {
                        "isolated_repetitions": {
                            "count": int(measured_attempt),
                            "wall_sec_values": [runtime_observation] if measured_attempt else [],
                            "wall_sec_median": runtime_observation if measured_attempt else None,
                            "wall_sec_mean": runtime_observation if measured_attempt else None,
                            "wall_sec_std": 0.0 if measured_attempt else None,
                            "cpu_sec_values": [],
                            "cpu_sec_median": None,
                            "stopped_after_incomplete_or_failed_attempt": True,
                            "right_censored": bool(measurement.get("right_censored", False)),
                            "termination_reason": measurement.get("termination_reason"),
                        }
                    },
                    "memory": {
                        "measurement": "parent_sampled_worker_plus_descendants",
                        "peak_rss_mb_values": [measurement["peak_rss_mb"]],
                        "peak_rss_mb": measurement["peak_rss_mb"],
                        "peak_rss_mb_median": measurement["peak_rss_mb"],
                    },
                    "worker_measurements": all_runs,
                    **_result_identity_metadata(job),
                }
            run_result = outcome["payload"]["result"]
            if run_result.get("status") != "ok":
                run_result["worker_measurements"] = all_runs
                run_result.setdefault("timing", {})["isolated_repetitions"] = {
                    "count": 0 if is_warmup else 1,
                    "wall_sec_values": [] if is_warmup else [outcome["measurement"]["wall_sec"]],
                    "wall_sec_median": None if is_warmup else outcome["measurement"]["wall_sec"],
                    "wall_sec_mean": None if is_warmup else outcome["measurement"]["wall_sec"],
                    "wall_sec_std": None if is_warmup else 0.0,
                    "cpu_sec_values": [],
                    "cpu_sec_median": None,
                    "stopped_after_incomplete_or_failed_attempt": True,
                }
                run_result["memory"] = {
                    "measurement": "parent_sampled_worker_plus_descendants",
                    "peak_rss_mb_values": [outcome["measurement"]["peak_rss_mb"]],
                    "peak_rss_mb": outcome["measurement"]["peak_rss_mb"],
                    "peak_rss_mb_median": outcome["measurement"]["peak_rss_mb"],
                }
                return run_result
            if not is_warmup:
                measured.append({"result": run_result, "measurement": outcome["measurement"]})

        base = measured[0]["result"]
        signatures = [
            (
                item["result"].get("execution_domain", {}).get("signature")
                if self.config.execution_profile.strip().lower() == "operational_optcuts"
                else item["result"].get("comparison_domain", {}).get("signature")
            )
            for item in measured
        ]
        if len(set(signatures)) > 1:
            base["status"] = "failed"
            base["error"] = "Comparison domain changed across isolated repetitions."
        wall = np.asarray([item["measurement"]["wall_sec"] for item in measured], dtype=np.float64)
        memory = np.asarray([item["measurement"]["peak_rss_mb"] for item in measured], dtype=np.float64)
        cpu = np.asarray(
            [item["result"].get("timing", {}).get("end_to_end", {}).get("cpu_sec", float("nan")) for item in measured],
            dtype=np.float64,
        )
        base.setdefault("timing", {})["isolated_repetitions"] = {
            "count": int(len(measured)),
            "wall_sec_values": wall.tolist(),
            "wall_sec_median": float(np.median(wall)),
            "wall_sec_mean": float(np.mean(wall)),
            "wall_sec_std": float(np.std(wall, ddof=1)) if len(wall) > 1 else 0.0,
            "cpu_sec_values": cpu.tolist(),
            "cpu_sec_median": float(np.nanmedian(cpu)),
        }
        base["memory"] = {
            "measurement": "parent_sampled_worker_plus_descendants",
            "peak_rss_mb_values": memory.tolist(),
            "peak_rss_mb": float(np.max(memory)),
            "peak_rss_mb_median": float(np.median(memory)),
        }
        base["repeatability"] = self._repeatability_block(
            [item["result"] for item in measured],
            signatures=signatures,
        )
        base["worker_measurements"] = all_runs
        return base

    def _repeatability_block(self, results: List[Dict[str, object]], signatures: List[object]) -> Dict[str, object]:
        def _value(record: Dict[str, object], path: Tuple[str, ...]) -> float:
            current: object = record
            for key in path:
                if not isinstance(current, dict):
                    return float("nan")
                current = current.get(key)
            try:
                return float(current)
            except (TypeError, ValueError):
                return float("nan")

        methods = self._configured_methods()
        metrics = {}
        for method in methods:
            values = np.asarray(
                [_value(result, (method, "distortion", "mean")) for result in results],
                dtype=np.float64,
            )
            finite = values[np.isfinite(values)]
            metrics[method] = {
                "distortion_mean_values": values.tolist(),
                "finite_count": int(len(finite)),
                "mean": float(np.mean(finite)) if len(finite) else float("nan"),
                "std": float(np.std(finite, ddof=1)) if len(finite) > 1 else 0.0 if len(finite) else float("nan"),
                "range": float(np.ptp(finite)) if len(finite) else float("nan"),
            }
        return {
            "measured_repetition_count": int(len(results)),
            "comparison_domain_signatures": list(signatures),
            "domain_stable": len(set(signatures)) <= 1,
            "method_distortion": metrics,
        }

    def _execute_worker(self, job: Dict[str, object], *, run_index: int, is_warmup: bool) -> Dict[str, object]:
        log_root = Path(self.config.output_root) / self.config.worker_log_folder
        log_root.mkdir(parents=True, exist_ok=True)
        safe_name = Path(str(job["pdb"])).name.replace(" ", "_")
        role = "warmup" if is_warmup else "measured"
        stem = f"{safe_name}.{role}.{run_index:03d}"
        job_path = log_root / f"{stem}.job.json"
        result_path = log_root / f"{stem}.result.json"
        stdout_path = log_root / f"{stem}.stdout.log"
        stderr_path = log_root / f"{stem}.stderr.log"
        payload = {
            "config": asdict(self.config),
            "pdb_path": os.path.join(self.config.input_folder, str(job["pdb"])),
            "chain_a": job["chain_a"],
            "chain_b": job["chain_b"],
            "job_metadata": {
                **job,
                "benchmark_run_index": int(run_index),
                "benchmark_is_warmup": bool(is_warmup),
                "emit_provenance": bool(not is_warmup and run_index == int(self.config.warmup_runs)),
            },
        }
        dump_json_atomic(payload, job_path)
        # A resumed failed attempt reuses the same deterministic filename. An
        # old result must never be mistaken for output from the new process.
        result_path.unlink(missing_ok=True)

        command = [sys.executable, "-m", "topoppi.benchmarking.worker", str(job_path), str(result_path)]
        environment = os.environ.copy()
        source_root = str(Path(__file__).resolve().parents[2])
        existing_path = environment.get("PYTHONPATH", "")
        environment["PYTHONPATH"] = source_root + (os.pathsep + existing_path if existing_path else "")
        thread_count = str(int(self.config.threads_per_worker))
        for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
            environment[variable] = thread_count
        environment["PYTHONHASHSEED"] = str(int(self.config.random_seed))
        environment["MPLBACKEND"] = "Agg"
        environment["TOPOPPI_BENCHMARK_RUN_INDEX"] = str(int(run_index))
        cpu_affinity = self._worker_cpu_affinity(job)
        environment["TOPOPPI_CPU_AFFINITY"] = ",".join(str(cpu) for cpu in cpu_affinity)
        started = time.perf_counter()
        peak_rss = 0.0
        termination_reason = None
        termination_error = None
        censoring_threshold_sec = None
        censoring_event_elapsed_sec = None
        method_censoring = None
        with (
            open(stdout_path, "w", encoding="utf-8") as stdout_handle,
            open(
                stderr_path,
                "w",
                encoding="utf-8",
            ) as stderr_handle,
        ):
            process = subprocess.Popen(command, stdout=stdout_handle, stderr=stderr_handle, env=environment)
            applied_cpu_affinity = None
            if psutil is not None and cpu_affinity:
                try:
                    worker_process = psutil.Process(process.pid)
                    worker_process.cpu_affinity(cpu_affinity)
                    applied_cpu_affinity = worker_process.cpu_affinity()
                except (psutil.Error, OSError, ValueError) as exc:
                    self._terminate_process_tree(process)
                    raise RuntimeError(f"Could not apply worker CPU affinity: {exc}") from exc
            while process.poll() is None:
                if self.cancel_event is not None and self.cancel_event.is_set():
                    self._terminate_process_tree(process)
                    raise RuntimeError("Benchmark cancelled by user.")
                elapsed = time.perf_counter() - started
                peak_rss = max(peak_rss, self._process_tree_rss_mb(process.pid))
                if elapsed > float(self.config.worker_timeout_sec):
                    termination_reason = "timeout"
                    censoring_threshold_sec = float(self.config.worker_timeout_sec)
                    censoring_event_elapsed_sec = float(elapsed)
                    termination_error = f"Worker timed out after {self.config.worker_timeout_sec:.1f}s."
                    self._terminate_process_tree(process)
                    break
                memory_limit = self.config.worker_memory_limit_mb
                if memory_limit is not None and peak_rss > float(memory_limit):
                    termination_reason = "memory_limit"
                    censoring_event_elapsed_sec = float(elapsed)
                    termination_error = (
                        f"Worker exceeded the {float(memory_limit):.1f} MB RSS limit (observed {peak_rss:.1f} MB)."
                    )
                    self._terminate_process_tree(process)
                    break
                time.sleep(float(self.config.worker_poll_interval_sec))
            return_code = process.wait()
            peak_rss = max(peak_rss, self._process_tree_rss_mb(process.pid))

        if termination_reason is not None:
            worker_payload = {"status": "failed", "error": termination_error}
        elif result_path.exists():
            with open(result_path, "r", encoding="utf-8") as handle:
                worker_payload = json.load(handle)
            worker_payload = self._externalize_worker_payload(
                worker_payload,
                result_path,
                preserve_details=bool(not is_warmup and run_index == int(self.config.warmup_runs)),
            )
            method_censoring = _operational_method_censoring(worker_payload)
        else:
            termination_reason = "worker_exit_without_result"
            worker_payload = {"status": "failed", "error": f"Worker exited {return_code} without a result file."}
        observed_wall_sec = float(time.perf_counter() - started)
        reported_termination_reason = termination_reason
        if termination_reason is None and method_censoring is not None:
            reported_termination_reason = str(method_censoring["termination_reason"])
            censoring_threshold_sec = method_censoring.get("censoring_threshold_sec")
            censoring_event_elapsed_sec = method_censoring.get("censoring_event_elapsed_sec")
        runtime_observation_sec = (
            float(censoring_event_elapsed_sec) if censoring_event_elapsed_sec is not None else observed_wall_sec
        )
        measurement = {
            "run_index": int(run_index),
            "warmup": bool(is_warmup),
            "wall_sec": observed_wall_sec,
            "runtime_observation_sec": runtime_observation_sec,
            "peak_rss_mb": float(peak_rss),
            "return_code": int(return_code),
            "worker_completed": termination_reason is None,
            "right_censored": bool(termination_reason in {"timeout", "memory_limit"} or method_censoring is not None),
            "termination_reason": reported_termination_reason,
            "censoring_threshold_sec": censoring_threshold_sec,
            "censoring_event_elapsed_sec": censoring_event_elapsed_sec,
            "job_file": str(job_path),
            "result_file": str(result_path),
            "stdout_log": str(stdout_path),
            "stderr_log": str(stderr_path),
            "command": [sys.executable, "-m", "topoppi.benchmarking.worker", "<job.json>", "<result.json>"],
            "cpu_affinity_requested": cpu_affinity,
            "cpu_affinity_applied": applied_cpu_affinity,
        }
        return {"payload": worker_payload, "measurement": measurement}

    def _externalize_worker_payload(
        self,
        worker_payload: Dict[str, object],
        result_path: Path,
        *,
        preserve_details: bool,
    ) -> Dict[str, object]:
        """Move large row-level evidence out of the structure summary payload."""

        result = worker_payload.get("result")
        if not isinstance(result, dict):
            return worker_payload

        per_residue_records = []

        def append_residue_records(
            methods: object,
            *,
            evidence_domain: str,
            domain_signature: object,
        ) -> None:
            if not isinstance(methods, dict):
                return
            for method in sorted(methods):
                report = methods[method]
                if not isinstance(report, dict):
                    continue
                footprint = report.get("residue_footprint_fragmentation", report)
                if not isinstance(footprint, dict):
                    continue
                for record in footprint.get("residues", []):
                    if isinstance(record, dict):
                        per_residue_records.append(
                            {
                                "evidence_domain": evidence_domain,
                                "domain_signature": domain_signature,
                                "method": method,
                                **record,
                            }
                        )

        fragmentation = result.get("residue_footprint_fragmentation")
        fragmentation_methods = fragmentation.get("methods", {}) if isinstance(fragmentation, dict) else {}
        standard_signature = (result.get("comparison_domain") or {}).get("signature")
        initialization_signature = (result.get("initialization_comparison_domain") or {}).get("signature")
        residue_aware_signature = (result.get("residue_aware_comparison_domain") or {}).get("signature")
        if isinstance(fragmentation_methods, dict):
            for method, report in fragmentation_methods.items():
                signature = (
                    initialization_signature
                    if method == "optcuts_lscm_initialized"
                    else residue_aware_signature
                    if method in RESIDUE_AWARE_OPTCUTS_METHODS
                    else standard_signature
                )
                append_residue_records(
                    {method: report},
                    evidence_domain="top_level_method_domain",
                    domain_signature=signature,
                )

        pair_quality = result.get("residue_aware_pair_quality")
        if isinstance(pair_quality, dict) and bool(pair_quality.get("complete", False)):
            append_residue_records(
                pair_quality.get("methods", {}),
                evidence_domain="residue_aware_exact_pair",
                domain_signature=pair_quality.get("domain_signature"),
            )

        independent_arms = result.get("independent_optcuts_arm_quality")
        if isinstance(independent_arms, dict):
            for method, arm in independent_arms.items():
                if not isinstance(arm, dict) or not isinstance(arm.get("quality"), dict):
                    continue
                append_residue_records(
                    {method: arm["quality"]},
                    evidence_domain="independent_complete_arm",
                    domain_signature=arm.get("domain_signature"),
                )

        method_execution = result.get("method_execution")
        execution_details = {}
        if isinstance(method_execution, dict):
            for method, block in method_execution.items():
                if isinstance(block, dict) and isinstance(block.get("executions"), list):
                    execution_details[str(method)] = block.pop("executions")

        topology_ablation_execution = None
        topology_ablation = result.get("topology_preprocessing_ablation")
        if isinstance(topology_ablation, dict):
            topology_ablation_execution = topology_ablation.pop("execution", None)

        detail = {
            "schema_version": "1.0",
            "pdb": result.get("pdb"),
            "per_face_sample_records": result.pop("per_face_sample_records", []),
            "per_residue_records": per_residue_records,
            "provenance_records": result.pop("provenance_records", []),
            "optcuts_execution_details": execution_details,
            "topology_ablation_execution": topology_ablation_execution,
        }
        self._remove_nested_residue_rows(result)

        if preserve_details:
            stem = result_path.name.removesuffix(".result.json")
            detail_path = result_path.with_name(f"{stem}.detail.json.gz")
            self._dump_gzip_json_atomic(detail, detail_path)
            result["detail_artifact"] = {
                "schema_version": "1.0",
                "path": os.path.relpath(detail_path, self.config.output_root),
                "sha256": sha256_file(detail_path),
                "per_face_sample_record_count": int(len(detail["per_face_sample_records"])),
                "per_residue_record_count": int(len(detail["per_residue_records"])),
                "provenance_record_count": int(len(detail["provenance_records"])),
                "optcuts_method_count": int(len(execution_details)),
            }

        dump_json_atomic(worker_payload, result_path)
        return worker_payload

    @classmethod
    def _remove_nested_residue_rows(cls, value: object) -> None:
        if isinstance(value, dict):
            value.pop("residues", None)
            for item in value.values():
                cls._remove_nested_residue_rows(item)
        elif isinstance(value, list):
            for item in value:
                cls._remove_nested_residue_rows(item)

    @staticmethod
    def _dump_gzip_json_atomic(payload: object, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
        try:
            with gzip.open(temporary, "wt", encoding="utf-8") as handle:
                json.dump(
                    json_safe(payload),
                    handle,
                    allow_nan=False,
                    sort_keys=True,
                    separators=(",", ":"),
                )
            os.replace(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)

    @staticmethod
    def _ordered_available_cpus() -> List[int]:
        if hasattr(os, "sched_getaffinity"):
            available = sorted(os.sched_getaffinity(0))
        elif psutil is not None:
            process = psutil.Process()
            if not hasattr(process, "cpu_affinity"):
                return []
            available = sorted(process.cpu_affinity())
        else:
            return []

        core_groups: Dict[Tuple[int, int], List[int]] = {}
        for cpu in available:
            topology_root = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
            try:
                package = int((topology_root / "physical_package_id").read_text(encoding="ascii").strip())
                core = int((topology_root / "core_id").read_text(encoding="ascii").strip())
            except (OSError, ValueError):
                return available
            core_groups.setdefault((package, core), []).append(cpu)

        ordered: List[int] = []
        groups = [sorted(core_groups[key]) for key in sorted(core_groups)]
        for sibling_index in range(max((len(group) for group in groups), default=0)):
            ordered.extend(group[sibling_index] for group in groups if sibling_index < len(group))
        return ordered

    def _worker_cpu_affinity(self, _job: Dict[str, object]) -> List[int]:
        available = self._available_worker_cpus
        count = min(int(self.config.threads_per_worker), len(available))
        if count <= 0:
            return []
        thread_id = threading.get_ident()
        with self._affinity_lock:
            assigned = self._worker_affinity_by_thread.get(thread_id)
            if assigned is not None:
                return list(assigned)
            start = self._next_affinity_slot * count
            if start + count > len(available):
                raise RuntimeError("Benchmark worker allocation exceeds the CPUs available to the supervisor.")
            assigned = list(available[start : start + count])
            self._worker_affinity_by_thread[thread_id] = assigned
            self._next_affinity_slot += 1
            return list(assigned)

    @staticmethod
    def _terminate_process_tree(process: subprocess.Popen) -> None:
        if psutil is not None:
            try:
                parent = psutil.Process(process.pid)
                children = parent.children(recursive=True)
                for child in children:
                    child.terminate()
                parent.terminate()
                _gone, alive = psutil.wait_procs([parent, *children], timeout=3.0)
                for item in alive:
                    item.kill()
                return
            except (psutil.Error, OSError):
                pass
        process.terminate()
        try:
            process.wait(timeout=3.0)
        except subprocess.TimeoutExpired:
            process.kill()

    @staticmethod
    def _process_tree_rss_mb(pid: int) -> float:
        if psutil is None:
            return 0.0
        try:
            parent = psutil.Process(pid)
            processes = [parent, *parent.children(recursive=True)]
            return float(sum(item.memory_info().rss for item in processes if item.is_running())) / (1024.0 * 1024.0)
        except (psutil.Error, OSError):
            return 0.0

    def _resolve_worker_count(self, file_count: int) -> int:
        if file_count == 0:
            return 0
        if self._available_worker_cpus:
            capacity = max(1, len(self._available_worker_cpus) // int(self.config.threads_per_worker))
        else:
            capacity = 1
        configured = capacity if self.config.max_workers is None else min(int(self.config.max_workers), capacity)
        return max(1, min(int(configured), int(file_count)))

    def _run_single(
        self,
        pdb_path: str,
        chain_a: str,
        chain_b: str,
        job_metadata: Dict[str, object],
    ) -> Dict[str, object]:
        self._check_cancelled()
        overall_wall = time.perf_counter()
        overall_cpu = _worker_cpu_time()
        self._log_thread(f"Start processing {os.path.basename(pdb_path)} ({chain_a}/{chain_b})")
        stage = {}

        stage_wall, stage_cpu = time.perf_counter(), _worker_cpu_time()
        loader = PDBLoader(pdb_path)
        coords_a, atoms_a = loader.get_chain_data(chain_a)
        coords_b, atoms_b = loader.get_chain_data(chain_b)
        stage["load_structure"] = self._stage_stats(stage_wall, stage_cpu)

        stage_wall, stage_cpu = time.perf_counter(), _worker_cpu_time()
        surface_generator = SurfaceGenerator(coords_a, config=self.config.surface)
        mesh_a = surface_generator.generate_mesh()
        if mesh_a is None or len(mesh_a.vertices) == 0:
            raise ValueError(
                f"Surface generation failed for {os.path.basename(pdb_path)}: {surface_generator.last_report}"
            )
        stage["surface_generation"] = self._stage_stats(stage_wall, stage_cpu)

        stage_wall, stage_cpu = time.perf_counter(), _worker_cpu_time()
        topology_manager = TopologyManager(mesh_a, coords_b, config=self.config.topology)
        patches = topology_manager.get_interface_patches()
        stage["interface_and_topology"] = self._stage_stats(stage_wall, stage_cpu)
        if not patches:
            preparation = {"attempted": 0, "success": 0, "failed": 0, "patches": []}
            patch_records = self._patch_retention_records(
                [],
                atoms_a,
                coords_a,
                atoms_b,
                coords_b,
                extracted_patches=[],
                topology_components=topology_manager.component_provenance,
                preparation=preparation,
                job_metadata=job_metadata,
            )
            return self._failed_single_result(
                pdb_path=pdb_path,
                chain_a=chain_a,
                chain_b=chain_b,
                job_metadata=job_metadata,
                error="No interface patch passed topology extraction",
                mesh_a=mesh_a,
                surface_report=surface_generator.last_report,
                topology_report=topology_manager.last_report,
                preparation=preparation,
                patch_records=patch_records,
                stage=stage,
                overall_wall=overall_wall,
                overall_cpu=overall_cpu,
                raw_patch_count=0,
            )

        parameterizer = Parameterizer(config=self.config.parameterization)
        stage_wall, stage_cpu = time.perf_counter(), _worker_cpu_time()
        prepared_patches, preparation = self._prepare_parameterization_domains(patches, parameterizer)
        stage["parameterization_domain_preparation"] = self._stage_stats(
            stage_wall,
            stage_cpu,
        )
        if not prepared_patches:
            patch_records = self._patch_retention_records(
                [],
                atoms_a,
                coords_a,
                atoms_b,
                coords_b,
                extracted_patches=patches,
                topology_components=topology_manager.component_provenance,
                preparation=preparation,
                job_metadata=job_metadata,
            )
            return self._failed_single_result(
                pdb_path=pdb_path,
                chain_a=chain_a,
                chain_b=chain_b,
                job_metadata=job_metadata,
                error="No patch passed the shared topology preparation gate",
                mesh_a=mesh_a,
                surface_report=surface_generator.last_report,
                topology_report=topology_manager.last_report,
                preparation=preparation,
                patch_records=patch_records,
                stage=stage,
                overall_wall=overall_wall,
                overall_cpu=overall_cpu,
                raw_patch_count=len(patches),
            )

        if self.config.execution_profile.strip().lower() == "operational_optcuts":
            return self._run_operational_optcuts(
                pdb_path=pdb_path,
                chain_a=chain_a,
                chain_b=chain_b,
                job_metadata=job_metadata,
                mesh_a=mesh_a,
                surface_report=surface_generator.last_report,
                topology_report=topology_manager.last_report,
                patches=patches,
                prepared_patches=prepared_patches,
                preparation=preparation,
                atoms_a=atoms_a,
                coords_a=coords_a,
                atoms_b=atoms_b,
                coords_b=coords_b,
                stage=stage,
                overall_wall=overall_wall,
                overall_cpu=overall_cpu,
            )

        stage_wall, stage_cpu = time.perf_counter(), _worker_cpu_time()
        tree_b = cKDTree(coords_b) if len(coords_b) else None
        contact_pairs_by_atom = self._geometric_contact_pairs_by_atom(
            range(len(atoms_a)),
            atoms_a,
            coords_a,
            atoms_b,
            tree_b,
        )
        (
            source_residue_labels,
            interaction_weights,
            objective_weights,
            interaction_source,
            contact_weight_definition,
        ) = self._residue_objective(
            job_metadata=job_metadata,
            chain_a=chain_a,
            chain_b=chain_b,
            atoms_a=atoms_a,
            coords_a=coords_a,
            atoms_b=atoms_b,
            coords_b=coords_b,
        )
        stage["contact_weight_preparation"] = self._stage_stats(stage_wall, stage_cpu)

        patch_results = {}
        method_execution = {}
        for method in PARAMETERIZATION_METHODS:
            self._check_cancelled()
            method_patches, wall_times, cpu_times, diagnostics = self._parameterize_patches(
                prepared_patches,
                method=method,
                parameterizer=parameterizer,
            )
            patch_results[method] = method_patches
            method_execution[method] = diagnostics
            stage[f"{method}_parameterization"] = self._from_timing_list(wall_times, cpu_times)

        active_optcuts = self.config.resolved_optcuts_variants()
        standard_optcuts_patches: Dict[str, List[object]] = {}
        residue_aware_method_patches: Dict[str, List[object]] = {}
        baseline_optcuts_config = replace(self.config.optcuts, residue_fragmentation_weight=0.0)
        optcuts_optimizer = OptCutsUVOptimizer(baseline_optcuts_config, cancel_event=self.cancel_event)
        optcuts_optimizer.preflight_binary()
        for method in STANDARD_OPTCUTS_METHODS:
            if method not in active_optcuts:
                continue
            initialization = "provided" if method.endswith("lscm_initialized") else "automatic"
            inputs = patch_results["lscm"] if initialization == "provided" else prepared_patches
            stage_wall, stage_cpu = time.perf_counter(), _worker_cpu_time()
            method_patches, diagnostics = self._run_optcuts(
                inputs,
                initialization=initialization,
                optimizer=optcuts_optimizer,
                source_residue_labels=source_residue_labels,
                residue_weights=objective_weights,
            )
            stage[method] = self._stage_stats(stage_wall, stage_cpu)
            method_execution[method] = diagnostics
            standard_optcuts_patches[method] = method_patches

        active_residue_aware = tuple(method for method in RESIDUE_AWARE_OPTCUTS_METHODS if method in active_optcuts)
        if active_residue_aware:
            residue_aware_optimizer = OptCutsUVOptimizer(self.config.optcuts, cancel_event=self.cancel_event)
            residue_aware_optimizer.preflight_binary()
            for method in active_residue_aware:
                initialization = "provided" if method.endswith("lscm_initialized") else "automatic"
                inputs = patch_results["lscm"] if initialization == "provided" else prepared_patches
                stage_wall, stage_cpu = time.perf_counter(), _worker_cpu_time()
                method_patches, diagnostics = self._run_optcuts(
                    inputs,
                    initialization=initialization,
                    optimizer=residue_aware_optimizer,
                    source_residue_labels=source_residue_labels,
                    residue_weights=objective_weights,
                )
                stage[method] = self._stage_stats(stage_wall, stage_cpu)
                method_execution[method] = diagnostics
                residue_aware_method_patches[method] = method_patches

        topology_ablation = {"status": "disabled"}
        topology_pair_quality: Dict[str, object] = {
            "status": "disabled",
            "complete": False,
            "methods": {},
        }
        if self.config.include_topology_ablation:
            stage_wall, stage_cpu = time.perf_counter(), _worker_cpu_time()
            raw_domain_optcuts, raw_domain_diag = self._run_optcuts(
                patches,
                initialization="automatic",
                optimizer=optcuts_optimizer,
                source_residue_labels=source_residue_labels,
                residue_weights=objective_weights,
            )
            stage["optcuts_without_parameterization_topology_gate"] = self._stage_stats(stage_wall, stage_cpu)
            raw_quality = quality_block(
                raw_domain_optcuts,
                patch_gap=self.config.optcuts.patch_gap,
                uv_key="uv_optcuts",
            )
            prepared_auto = standard_optcuts_patches["optcuts_automatic"]
            raw_inputs = {str(patch.metadata.get("patch_id", "unknown")): patch for patch in patches}
            prepared_inputs = {str(patch.metadata.get("patch_id", "unknown")): patch for patch in prepared_patches}
            raw_outputs = {str(patch.metadata.get("patch_id", "unknown")): patch for patch in raw_domain_optcuts}
            prepared_outputs = {str(patch.metadata.get("patch_id", "unknown")): patch for patch in prepared_auto}
            expected_ids = set(raw_inputs)
            exact_input_ids = {
                patch_id
                for patch_id in expected_ids & set(prepared_inputs)
                if face_domain_hash(raw_inputs[patch_id]) == face_domain_hash(prepared_inputs[patch_id])
            }
            exact_output_ids = {
                patch_id
                for patch_id in exact_input_ids & set(raw_outputs) & set(prepared_outputs)
                if face_domain_hash(raw_outputs[patch_id]) == face_domain_hash(raw_inputs[patch_id])
                and face_domain_hash(prepared_outputs[patch_id]) == face_domain_hash(raw_inputs[patch_id])
            }
            exact_ids = sorted(exact_output_ids)
            unique_ids = all(
                len(mapping) == len(collection)
                for mapping, collection in (
                    (raw_inputs, patches),
                    (prepared_inputs, prepared_patches),
                    (raw_outputs, raw_domain_optcuts),
                    (prepared_outputs, prepared_auto),
                )
            )
            topology_pair_complete = bool(
                unique_ids
                and expected_ids
                and exact_output_ids == expected_ids
                and set(prepared_inputs) == expected_ids
                and set(raw_outputs) == expected_ids
                and set(prepared_outputs) == expected_ids
            )
            paired_raw = [raw_outputs[patch_id] for patch_id in exact_ids]
            paired_prepared = [prepared_outputs[patch_id] for patch_id in exact_ids]
            topology_pair_methods = (
                {
                    "optcuts_without_topology_preparation": {
                        **quality_block(
                            paired_raw,
                            patch_gap=self.config.optcuts.patch_gap,
                            uv_key="uv_optcuts",
                        ),
                        "runtime": stage["optcuts_without_parameterization_topology_gate"],
                        "preprocessing_runtime": {"wall_sec": 0.0, "cpu_sec": 0.0},
                        "optcuts_runtime": stage["optcuts_without_parameterization_topology_gate"],
                    },
                    "optcuts_with_topology_preparation": {
                        **quality_block(
                            paired_prepared,
                            patch_gap=self.config.optcuts.patch_gap,
                            uv_key="uv_optcuts",
                        ),
                        "runtime": self._sum_stage_stats(
                            stage["parameterization_domain_preparation"],
                            stage["optcuts_automatic"],
                        ),
                        "preprocessing_runtime": stage["parameterization_domain_preparation"],
                        "optcuts_runtime": stage["optcuts_automatic"],
                    },
                }
                if exact_ids
                else {}
            )
            topology_pair_quality = {
                "status": "evaluated" if topology_pair_complete else "ineligible_or_incomplete",
                "complete": topology_pair_complete,
                "expected_patch_count": int(len(expected_ids)),
                "exact_pair_patch_count": int(len(exact_ids)),
                "exact_pair_patch_ids": exact_ids,
                "raw_success_patch_count": int(len(raw_outputs)),
                "prepared_success_patch_count": int(len(prepared_outputs)),
                "exact_source_face_match": bool(exact_input_ids == expected_ids),
                "unique_patch_ids": bool(unique_ids),
                "domain_signature": hashlib.sha256(
                    json.dumps(
                        {patch_id: face_domain_hash(raw_inputs[patch_id]) for patch_id in exact_ids},
                        sort_keys=True,
                    ).encode("utf-8")
                ).hexdigest(),
                "methods": topology_pair_methods,
                "rule": (
                    "structure-level paired ablation is eligible only when every extracted patch "
                    "succeeds in both arms and retains exactly the same source faces and geometry"
                ),
                "runtime_rule": (
                    "the prepared-arm wall and CPU times include parameterization-domain "
                    "preparation plus OptCuts; the raw arm has no corresponding preprocessing"
                ),
            }
            topology_ablation = {
                "status": "evaluated",
                "comparison_scope": (
                    "separate prespecified topology-preparation ablation; exact-domain pairs are "
                    "reported independently from the primary method comparison"
                ),
                "before_gate_patch_count": int(len(patches)),
                "after_gate_patch_count": int(len(prepared_patches)),
                "execution": raw_domain_diag,
                "quality_before_gate": raw_quality,
            }

        primary_standard_method_patches = {
            **patch_results,
            **{
                method: method_patches
                for method, method_patches in standard_optcuts_patches.items()
                if method != "optcuts_lscm_initialized"
            },
        }
        standard_common_ids = self._common_patch_ids(primary_standard_method_patches)
        standard_common = {
            method: self._filter_patch_ids(method_patches, standard_common_ids)
            for method, method_patches in primary_standard_method_patches.items()
        }
        comparison_signature = self._comparison_signature(standard_common_ids, prepared_patches)
        expected_hashes = {
            str(patch.metadata.get("patch_id", "unknown")): face_domain_hash(patch) for patch in prepared_patches
        }
        expected_ids = set(expected_hashes)
        standard_id_sets = {
            method: [str(patch.metadata.get("patch_id", "unknown")) for patch in method_patches]
            for method, method_patches in primary_standard_method_patches.items()
        }
        standard_unique_ids = all(len(ids) == len(set(ids)) for ids in standard_id_sets.values())
        exact_domain_match = all(
            all(
                face_domain_hash(patch) == expected_hashes.get(str(patch.metadata.get("patch_id", "unknown")))
                for patch in method_patches
            )
            for method_patches in standard_common.values()
        )
        domain_complete = bool(
            standard_unique_ids
            and standard_common_ids == expected_ids
            and all(set(ids) == expected_ids for ids in standard_id_sets.values())
            and exact_domain_match
        )

        qualities = {
            method: quality_block(
                standard_common[method],
                patch_gap=self.config.optcuts.patch_gap,
                uv_key="uv_optcuts" if _uses_optcuts_uv(method) else "uv",
            )
            for method in primary_standard_method_patches
        }
        metric_complete = all(
            np.isfinite(float(qualities[method]["distortion"]["mean"])) for method in primary_standard_method_patches
        )
        comparison_complete = bool(domain_complete and metric_complete)

        # Every baseline has its own exact, all-expected-patch comparison with
        # automatic OptCuts.  An unrelated method failure must not remove an
        # otherwise complete pair, while an arm failure remains visible in the
        # reliability endpoint.
        complete_method_quality: Dict[str, Dict[str, object]] = {}
        standard_arm_patches: Dict[str, List[object]] = {}
        method_domain_status: Dict[str, Dict[str, object]] = {}
        for method, method_patches in primary_standard_method_patches.items():
            ids = standard_id_sets[method]
            method_map = {str(patch.metadata.get("patch_id", "unknown")): patch for patch in method_patches}
            unique_ids = len(ids) == len(method_map)
            exact_ids = {
                patch_id
                for patch_id, patch in method_map.items()
                if patch_id in expected_hashes and face_domain_hash(patch) == expected_hashes[patch_id]
            }
            arm_domain_complete = bool(unique_ids and set(method_map) == expected_ids and exact_ids == expected_ids)
            if arm_domain_complete:
                ordered_arm_patches = [method_map[patch_id] for patch_id in sorted(expected_ids)]
                if domain_complete:
                    block = qualities[method]
                else:
                    block = quality_block(
                        ordered_arm_patches,
                        patch_gap=self.config.optcuts.patch_gap,
                        uv_key="uv_optcuts" if _uses_optcuts_uv(method) else "uv",
                    )
                complete_method_quality[method] = block
                standard_arm_patches[method] = ordered_arm_patches
            else:
                block = None
            method_domain_status[method] = {
                "domain_complete": arm_domain_complete,
                "unique_patch_ids": bool(unique_ids),
                "successful_patch_count": int(len(method_patches)),
                "exact_source_face_patch_count": int(len(exact_ids)),
                "metric_finite": bool(block is not None and np.isfinite(float(block["distortion"]["mean"]))),
                "globally_injective": bool(
                    block is not None
                    and (block.get("injectivity") or {}).get(
                        "all_patches_globally_injective",
                        False,
                    )
                ),
            }
            method_domain_status[method]["analysis_eligible"] = bool(
                method_domain_status[method]["domain_complete"] and method_domain_status[method]["metric_finite"]
            )
            method_domain_status[method]["usable"] = bool(
                method_domain_status[method]["analysis_eligible"] and method_domain_status[method]["globally_injective"]
            )

        standard_method_pair_quality: Dict[str, Dict[str, object]] = {}
        treatment = "optcuts_automatic"
        if treatment in primary_standard_method_patches:
            for baseline in PARAMETERIZATION_METHODS:
                if baseline not in primary_standard_method_patches:
                    continue
                pair_complete = bool(
                    method_domain_status[baseline]["domain_complete"]
                    and method_domain_status[treatment]["domain_complete"]
                    and method_domain_status[baseline]["metric_finite"]
                    and method_domain_status[treatment]["metric_finite"]
                )
                pair_ids = set(standard_id_sets[baseline]) & set(standard_id_sets[treatment])
                pair_key = f"{baseline}_vs_{treatment}"
                standard_method_pair_quality[pair_key] = {
                    "status": "evaluated" if pair_complete else "ineligible_or_incomplete",
                    "complete": pair_complete,
                    "baseline": baseline,
                    "treatment": treatment,
                    "expected_patch_count": int(len(expected_ids)),
                    "common_patch_count": int(len(pair_ids & expected_ids)),
                    "domain_signature": self._comparison_signature(expected_ids, prepared_patches),
                    "arms": {
                        baseline: dict(method_domain_status[baseline]),
                        treatment: dict(method_domain_status[treatment]),
                    },
                    "methods": (
                        {
                            baseline: self._paired_quality_projection(complete_method_quality[baseline]),
                            treatment: self._paired_quality_projection(complete_method_quality[treatment]),
                        }
                        if pair_complete
                        else {}
                    ),
                    "rule": (
                        "both arms must contain every expected patch exactly once, preserve "
                        "the exact source-face geometry, and have finite distortion metrics"
                    ),
                }

        initialization_common_ids: set[str] = set()
        initialization_common: Dict[str, List[object]] = {}
        initialization_comparison_signature = ""
        initialization_exact_domain_match = False
        initialization_domain_complete = False
        initialization_metric_complete = False
        initialization_pair_qualities: Dict[str, Dict[str, object]] = {}
        if "optcuts_lscm_initialized" in standard_optcuts_patches:
            initialization_methods = {
                "lscm": patch_results["lscm"],
                "optcuts_automatic": standard_optcuts_patches["optcuts_automatic"],
                "optcuts_lscm_initialized": standard_optcuts_patches["optcuts_lscm_initialized"],
            }
            initialization_common_ids = self._common_patch_ids(initialization_methods)
            initialization_common = {
                method: self._filter_patch_ids(method_patches, initialization_common_ids)
                for method, method_patches in initialization_methods.items()
            }
            initialization_comparison_signature = self._comparison_signature(
                initialization_common_ids,
                prepared_patches,
            )
            initialization_exact_domain_match = all(
                all(
                    face_domain_hash(patch) == expected_hashes.get(str(patch.metadata.get("patch_id", "unknown")))
                    for patch in method_patches
                )
                for method_patches in initialization_common.values()
            )
            initialization_id_lists = [
                [str(patch.metadata.get("patch_id", "unknown")) for patch in method_patches]
                for method_patches in initialization_methods.values()
            ]
            initialization_domain_complete = bool(
                initialization_common_ids == expected_ids
                and all(len(ids) == len(set(ids)) for ids in initialization_id_lists)
                and all(set(ids) == expected_ids for ids in initialization_id_lists)
                and initialization_exact_domain_match
            )
            initialization_pair_qualities = {
                method: quality_block(
                    method_patches,
                    patch_gap=self.config.optcuts.patch_gap,
                    uv_key="uv_optcuts" if _uses_optcuts_uv(method) else "uv",
                )
                for method, method_patches in initialization_common.items()
            }
            initialization_metric_complete = bool(initialization_pair_qualities) and all(
                np.isfinite(float(block["distortion"]["mean"])) for block in initialization_pair_qualities.values()
            )
            qualities["optcuts_lscm_initialized"] = initialization_pair_qualities["optcuts_lscm_initialized"]
        initialization_comparison_complete = bool(initialization_domain_complete and initialization_metric_complete)

        residue_aware_common_ids: set[str] = set()
        residue_aware_common: Dict[str, List[object]] = {}
        residue_aware_domain_complete = False
        residue_aware_metric_complete = False
        residue_aware_comparison_signature = ""
        residue_aware_exact_domain_match = False
        residue_aware_domain_qualities: Dict[str, Dict[str, object]] = {}
        residue_aware_arm_status: Dict[str, Dict[str, object]] = {}
        residue_aware_complete_method_quality: Dict[str, Dict[str, object]] = {}
        residue_aware_arm_patches: Dict[str, List[object]] = {}
        if residue_aware_method_patches:
            residue_aware_domain_methods = {
                **{
                    RESIDUE_AWARE_BASELINE[method]: standard_optcuts_patches[RESIDUE_AWARE_BASELINE[method]]
                    for method in residue_aware_method_patches
                },
                **residue_aware_method_patches,
            }
            residue_aware_common_ids = self._common_patch_ids(residue_aware_domain_methods)
            residue_aware_common = {
                method: self._filter_patch_ids(method_patches, residue_aware_common_ids)
                for method, method_patches in residue_aware_domain_methods.items()
            }
            residue_aware_comparison_signature = self._comparison_signature(residue_aware_common_ids, prepared_patches)
            residue_aware_exact_domain_match = all(
                all(
                    face_domain_hash(patch) == expected_hashes.get(str(patch.metadata.get("patch_id", "unknown")))
                    for patch in method_patches
                )
                for method_patches in residue_aware_common.values()
            )
            residue_aware_id_lists = [
                [str(patch.metadata.get("patch_id", "unknown")) for patch in method_patches]
                for method_patches in residue_aware_domain_methods.values()
            ]
            residue_aware_domain_complete = bool(
                residue_aware_common_ids == expected_ids
                and all(len(ids) == len(set(ids)) for ids in residue_aware_id_lists)
                and all(set(ids) == expected_ids for ids in residue_aware_id_lists)
                and residue_aware_exact_domain_match
            )
            residue_aware_domain_qualities = {
                method: quality_block(
                    residue_aware_common[method],
                    patch_gap=self.config.optcuts.patch_gap,
                    uv_key="uv_optcuts",
                )
                for method in residue_aware_domain_methods
            }
            for method, method_patches in residue_aware_domain_methods.items():
                method_ids = [str(patch.metadata.get("patch_id", "unknown")) for patch in method_patches]
                method_map = {str(patch.metadata.get("patch_id", "unknown")): patch for patch in method_patches}
                unique_ids = len(method_ids) == len(method_map)
                exact_ids = {
                    patch_id
                    for patch_id, patch in method_map.items()
                    if patch_id in expected_hashes and face_domain_hash(patch) == expected_hashes[patch_id]
                }
                arm_domain_complete = bool(unique_ids and set(method_map) == expected_ids and exact_ids == expected_ids)
                ordered_arm_patches = (
                    [method_map[patch_id] for patch_id in sorted(expected_ids)] if arm_domain_complete else []
                )
                arm_quality = (
                    quality_block(
                        ordered_arm_patches,
                        patch_gap=self.config.optcuts.patch_gap,
                        uv_key="uv_optcuts",
                    )
                    if arm_domain_complete
                    else None
                )
                metric_finite = bool(arm_quality is not None and np.isfinite(float(arm_quality["distortion"]["mean"])))
                globally_injective = bool(
                    arm_quality is not None
                    and (arm_quality.get("injectivity") or {}).get(
                        "all_patches_globally_injective",
                        False,
                    )
                )
                residue_aware_arm_status[method] = {
                    "domain_complete": arm_domain_complete,
                    "unique_patch_ids": bool(unique_ids),
                    "successful_patch_count": int(len(method_patches)),
                    "exact_source_face_patch_count": int(len(exact_ids)),
                    "metric_finite": metric_finite,
                    "globally_injective": globally_injective,
                    "analysis_eligible": bool(arm_domain_complete and metric_finite),
                    "usable": bool(arm_domain_complete and metric_finite and globally_injective),
                }
                if arm_quality is not None:
                    residue_aware_complete_method_quality[method] = arm_quality
                    residue_aware_arm_patches[method] = ordered_arm_patches
            for method in residue_aware_method_patches:
                qualities[method] = residue_aware_domain_qualities[method]
            residue_aware_metric_complete = all(
                np.isfinite(float(residue_aware_domain_qualities[method]["distortion"]["mean"]))
                for method in residue_aware_domain_methods
            )
        residue_aware_injectivity_complete = bool(residue_aware_arm_status) and all(
            bool(status.get("globally_injective", False)) for status in residue_aware_arm_status.values()
        )
        residue_aware_comparison_complete = bool(residue_aware_domain_complete and residue_aware_metric_complete)

        stage_wall, stage_cpu = time.perf_counter(), _worker_cpu_time()
        fragmentation_by_method = {
            method: residue_fragmentation_report(
                standard_common[method],
                source_residue_labels,
                uv_key="uv_optcuts" if _uses_optcuts_uv(method) else "uv",
                interaction_weights=interaction_weights,
                objective_weights=objective_weights,
            )
            for method in primary_standard_method_patches
        }
        if "optcuts_lscm_initialized" in initialization_common:
            fragmentation_by_method["optcuts_lscm_initialized"] = residue_fragmentation_report(
                initialization_common["optcuts_lscm_initialized"],
                source_residue_labels,
                uv_key="uv_optcuts",
                interaction_weights=interaction_weights,
                objective_weights=objective_weights,
            )
        fragmentation_by_method.update(
            {
                method: residue_fragmentation_report(
                    residue_aware_common[method],
                    source_residue_labels,
                    uv_key="uv_optcuts",
                    interaction_weights=interaction_weights,
                    objective_weights=objective_weights,
                )
                for method in residue_aware_method_patches
                if method in residue_aware_common
            }
        )
        residue_aware_domain_fragmentation: Dict[str, Dict[str, object]] = {}
        if residue_aware_method_patches:
            for method in residue_aware_common:
                if method in residue_aware_method_patches:
                    residue_aware_domain_fragmentation[method] = fragmentation_by_method[method]
                else:
                    residue_aware_domain_fragmentation[method] = residue_fragmentation_report(
                        residue_aware_common[method],
                        source_residue_labels,
                        uv_key="uv_optcuts",
                        interaction_weights=interaction_weights,
                        objective_weights=objective_weights,
                    )

        independent_optcuts_arm_quality: Dict[str, Dict[str, object]] = {}
        optcuts_arm_sources = {
            "optcuts_automatic": (
                method_domain_status.get("optcuts_automatic", {}),
                complete_method_quality.get("optcuts_automatic"),
                standard_arm_patches.get("optcuts_automatic"),
                standard_common_ids == expected_ids,
            ),
            "residue_aware_optcuts": (
                residue_aware_arm_status.get("residue_aware_optcuts", {}),
                residue_aware_complete_method_quality.get("residue_aware_optcuts"),
                residue_aware_arm_patches.get("residue_aware_optcuts"),
                residue_aware_common_ids == expected_ids,
            ),
        }
        for method, (status, quality, arm_patches, common_is_full) in optcuts_arm_sources.items():
            if not status:
                continue
            arm_quality = dict(quality) if isinstance(quality, dict) else None
            if arm_quality is not None and arm_patches is not None:
                footprint_report = (
                    fragmentation_by_method[method]
                    if common_is_full and method in fragmentation_by_method
                    else residue_fragmentation_report(
                        arm_patches,
                        source_residue_labels,
                        uv_key="uv_optcuts",
                        interaction_weights=interaction_weights,
                        objective_weights=objective_weights,
                    )
                )
                arm_quality["residue_footprint_fragmentation"] = footprint_report
            independent_optcuts_arm_quality[method] = {
                **status,
                "domain_signature": self._comparison_signature(expected_ids, prepared_patches),
                "quality": arm_quality,
                "rule": (
                    "quality uses every expected patch exactly once on unchanged source-face geometry "
                    "and is independent of failures in unrelated comparator methods"
                ),
            }
        for method, footprint_report in fragmentation_by_method.items():
            qualities[method]["residue_footprint_fragmentation"] = footprint_report
        stage["residue_footprint_fragmentation"] = self._stage_stats(stage_wall, stage_cpu)

        stage_wall, stage_cpu = time.perf_counter(), _worker_cpu_time()
        atlas_reference_method = next(
            method for method in ("optcuts_automatic", "harmonic", "lscm") if method in standard_common
        )
        atlas_reference_uv_key = "uv_optcuts" if _uses_optcuts_uv(atlas_reference_method) else "uv"
        atlas_reference_uv, _atlas_transforms, _atlas_packing = pack_mesh_charts(
            standard_common[atlas_reference_method],
            key=atlas_reference_uv_key,
            gap=self.config.optcuts.patch_gap,
        )
        atlas_map, patch_coverages = rasterize_feature_maps(
            standard_common[atlas_reference_method],
            size=self.config.raster_size,
            return_patch_coverage=True,
            uv_arrays=atlas_reference_uv,
        )
        atlas_trainability = atlas_trainability_metrics(atlas_map, patch_coverages)
        atlas_trainability["reference_method"] = atlas_reference_method
        stage["feature_rasterization"] = self._stage_stats(stage_wall, stage_cpu)

        energy_lscm = avg_energy(standard_common["lscm"], uv_key="uv")
        seam_lscm = avg_seam_length(standard_common["lscm"], uv_key="uv")
        optcuts_scoring_patches = {
            method: method_patches for method, method_patches in standard_common.items() if _uses_optcuts_uv(method)
        }
        if "optcuts_lscm_initialized" in initialization_common:
            optcuts_scoring_patches["optcuts_lscm_initialized"] = initialization_common["optcuts_lscm_initialized"]
        optcuts_scoring_patches.update(
            {
                method: method_patches
                for method, method_patches in residue_aware_common.items()
                if _uses_optcuts_uv(method)
            }
        )
        optcuts_energy = {
            method: avg_energy(method_patches, uv_key="uv_optcuts")
            for method in active_optcuts
            if (method_patches := optcuts_scoring_patches.get(method)) is not None
        }
        optcuts_seam = {
            method: avg_seam_length(method_patches, uv_key="uv_optcuts")
            for method in active_optcuts
            if (method_patches := optcuts_scoring_patches.get(method)) is not None
        }
        residue_aware_ablation: Dict[str, object] = {"status": "disabled"}
        residue_aware_pair_quality: Dict[str, object] = {
            "status": "disabled",
            "methods": {},
        }
        if residue_aware_method_patches:
            comparisons = {}
            if residue_aware_comparison_complete:
                for treatment in residue_aware_method_patches:
                    baseline = RESIDUE_AWARE_BASELINE[treatment]
                    initialization = "automatic"
                    comparisons[initialization] = {
                        "baseline": baseline,
                        "treatment": treatment,
                        "distortion_mean_baseline": residue_aware_domain_qualities[baseline]["distortion"]["mean"],
                        "distortion_mean_treatment": residue_aware_domain_qualities[treatment]["distortion"]["mean"],
                        "seam_baseline": avg_seam_length(residue_aware_common[baseline], uv_key="uv_optcuts"),
                        "seam_treatment": avg_seam_length(residue_aware_common[treatment], uv_key="uv_optcuts"),
                        "interaction_weighted_fragmentation_baseline": residue_aware_domain_fragmentation[baseline][
                            "interaction_weighted_fragmentation"
                        ],
                        "interaction_weighted_fragmentation_treatment": residue_aware_domain_fragmentation[treatment][
                            "interaction_weighted_fragmentation"
                        ],
                        "objective_weighted_fragmentation_baseline": residue_aware_domain_fragmentation[baseline][
                            "objective_weighted_fragmentation"
                        ],
                        "objective_weighted_fragmentation_treatment": residue_aware_domain_fragmentation[treatment][
                            "objective_weighted_fragmentation"
                        ],
                    }
            residue_aware_ablation = self._residue_aware_ablation_record(
                comparison_complete=residue_aware_comparison_complete,
                residue_fragmentation_weight=self.config.optcuts.residue_fragmentation_weight,
                comparisons=comparisons,
            )
            residue_aware_pair_quality = {
                "status": "evaluated" if residue_aware_comparison_complete else "incomplete_comparison",
                "complete": bool(residue_aware_comparison_complete),
                "expected_patch_count": int(len(expected_ids)),
                "common_patch_count": int(len(residue_aware_common_ids)),
                "domain_signature": residue_aware_comparison_signature,
                "arms": residue_aware_arm_status,
                "methods": (
                    {
                        method: {
                            **residue_aware_domain_qualities[method],
                            "residue_footprint_fragmentation": residue_aware_domain_fragmentation[method],
                        }
                        for method in residue_aware_common
                    }
                    if residue_aware_comparison_complete
                    else {}
                ),
                "rule": (
                    "both arms must return every expected patch exactly once on unchanged "
                    "source-face geometry with finite metrics; local flips and global overlap "
                    "remain in the efficacy domain and are reported as separate geometry QC"
                ),
            }
        initialization_ablation: Dict[str, object] = {"status": "not_evaluated"}
        initialization_pair_quality: Dict[str, object] = {
            "status": "not_evaluated",
            "methods": {},
        }
        if "optcuts_lscm_initialized" in standard_optcuts_patches:
            energy_automatic = avg_energy(
                initialization_common["optcuts_automatic"],
                uv_key="uv_optcuts",
            )
            energy_initialized = avg_energy(
                initialization_common["optcuts_lscm_initialized"],
                uv_key="uv_optcuts",
            )
            initialization_ablation = {
                "status": ("evaluated" if initialization_comparison_complete else "incomplete_comparison"),
                "baseline": "optcuts_automatic",
                "treatment": "optcuts_lscm_initialized",
                "same_domain": bool(initialization_comparison_complete),
                "same_optcuts_settings": True,
                "distortion_mean_automatic": initialization_pair_qualities["optcuts_automatic"]["distortion"]["mean"],
                "distortion_mean_lscm_initialized": initialization_pair_qualities["optcuts_lscm_initialized"][
                    "distortion"
                ]["mean"],
                "energy_automatic": energy_automatic,
                "energy_lscm_initialized": energy_initialized,
                "seam_automatic": avg_seam_length(initialization_common["optcuts_automatic"], uv_key="uv_optcuts"),
                "seam_lscm_initialized": avg_seam_length(
                    initialization_common["optcuts_lscm_initialized"], uv_key="uv_optcuts"
                ),
            }
            initialization_pair_quality = {
                "status": initialization_ablation["status"],
                "domain_signature": initialization_comparison_signature,
                "methods": initialization_pair_qualities,
                "rule": (
                    "free-boundary LSCM is supplied verbatim only when globally injective; "
                    "all values use the exact LSCM/direct/provided-OptCuts source-face intersection"
                ),
            }
        patch_records = self._patch_retention_records(
            prepared_patches,
            atoms_a,
            coords_a,
            atoms_b,
            coords_b,
            extracted_patches=patches,
            topology_components=topology_manager.component_provenance,
            preparation=preparation,
            job_metadata=job_metadata,
            contact_pairs_by_atom=contact_pairs_by_atom,
        )
        per_face_sample_records = self._per_face_sample_records(
            prepared_patches,
            {
                **patch_results,
                **standard_optcuts_patches,
                **residue_aware_method_patches,
            },
        )
        provenance_records = (
            self._provenance_records(prepared_patches) if bool(job_metadata.get("emit_provenance", True)) else []
        )

        result = {
            "pdb": os.path.basename(pdb_path),
            "input_sha256": str(job_metadata["input_sha256"]),
            "interaction_sha256": job_metadata.get("prolif_sha256") or None,
            "status": "ok" if comparison_complete else "incomplete_comparison",
            "chain_selection": {
                "chain_a": chain_a,
                "chain_b": chain_b,
                "mode": job_metadata.get("selection_mode", "configured"),
                "details": job_metadata.get("selection_details", {}),
            },
            **_result_identity_metadata(job_metadata),
            "benchmark_purpose": self.config.benchmark_purpose.strip().lower(),
            "execution_profile": self.config.execution_profile.strip().lower(),
            "topology_ablation_configured": bool(self.config.include_topology_ablation),
            "afdb_complex_confidence": _afdb_complex_confidence(job_metadata),
            "paired_geometry_qc": {
                "contact_cutoff_angstrom": _optional_finite_float(job_metadata.get("paired_contact_cutoff_angstrom")),
                "predicted_contact_count_total": _optional_nonnegative_int(
                    job_metadata.get("paired_predicted_contact_count_total")
                ),
                "contact_recall_fnat": _optional_finite_float(job_metadata.get("paired_contact_recall_fnat")),
                "contact_precision": _optional_finite_float(job_metadata.get("paired_contact_precision")),
                "contact_jaccard": _optional_finite_float(job_metadata.get("paired_contact_jaccard")),
                "experimental_contact_mapping_coverage": _optional_finite_float(
                    job_metadata.get("paired_experimental_contact_mapping_coverage")
                ),
                "interface_residue_a_mapping_coverage": _optional_finite_float(
                    job_metadata.get("paired_interface_residue_a_mapping_coverage")
                ),
                "interface_residue_b_mapping_coverage": _optional_finite_float(
                    job_metadata.get("paired_interface_residue_b_mapping_coverage")
                ),
                "interface_ligand_ca_mapping_coverage": _optional_finite_float(
                    job_metadata.get("paired_interface_ligand_ca_mapping_coverage")
                ),
                "interface_ligand_ca_rmsd_angstrom": _optional_finite_float(
                    job_metadata.get("paired_interface_ligand_ca_rmsd_angstrom")
                ),
                "cross_chain_clash_atom_fraction": _optional_finite_float(
                    job_metadata.get("paired_cross_chain_clash_atom_fraction")
                ),
            },
            "patch_count": int(len(patches)),
            "prepared_patch_count": int(len(prepared_patches)),
            "mesh_stats": {"vertex_count": int(len(mesh_a.vertices)), "face_count": int(len(mesh_a.faces))},
            "surface_generation": surface_generator.last_report,
            "comparison_domain": {
                "complete": bool(comparison_complete),
                "domain_complete": bool(domain_complete),
                "metric_complete": bool(metric_complete),
                "exact_source_face_match": bool(exact_domain_match),
                "unique_patch_ids": bool(standard_unique_ids),
                "expected_patch_count": int(len(prepared_patches)),
                "common_patch_count": int(len(standard_common_ids)),
                "common_patch_ids": sorted(standard_common_ids),
                "signature": comparison_signature,
                "rule": (
                    "intersection across the five standalone parameterizations and automatic "
                    "OptCuts; the feasibility-limited LSCM initialization diagnostic is separate"
                ),
            },
            "standard_method_pair_quality": standard_method_pair_quality,
            "initialization_comparison_domain": {
                "enabled": "optcuts_lscm_initialized" in standard_optcuts_patches,
                "complete": bool(initialization_comparison_complete),
                "domain_complete": bool(initialization_domain_complete),
                "metric_complete": bool(initialization_metric_complete),
                "exact_source_face_match": bool(initialization_exact_domain_match),
                "expected_patch_count": int(len(prepared_patches)),
                "common_patch_count": int(len(initialization_common_ids)),
                "common_patch_ids": sorted(initialization_common_ids),
                "signature": initialization_comparison_signature,
                "rule": (
                    "exact intersection of free-boundary LSCM, automatic OptCuts, and verbatim "
                    "LSCM-initialized OptCuts after the global-injectivity gate"
                ),
            },
            "residue_aware_comparison_domain": {
                "enabled": bool(residue_aware_method_patches),
                "complete": bool(residue_aware_comparison_complete),
                "domain_complete": bool(residue_aware_domain_complete),
                "metric_complete": bool(residue_aware_metric_complete),
                "injectivity_complete": bool(residue_aware_injectivity_complete),
                "exact_source_face_match": bool(residue_aware_exact_domain_match),
                "expected_patch_count": int(len(prepared_patches)),
                "common_patch_count": int(len(residue_aware_common_ids)),
                "common_patch_ids": sorted(residue_aware_common_ids),
                "signature": residue_aware_comparison_signature,
                "rule": (
                    "separate exact-domain finite-output intersection across standard and residue-aware "
                    "OptCuts methods; injectivity is reported independently as geometry QC"
                ),
            },
            "independent_optcuts_arm_quality": independent_optcuts_arm_quality,
            "method_execution": method_execution,
            "topology_extraction": topology_manager.last_report,
            "topology_preparation": preparation,
            "topology_preprocessing_ablation": topology_ablation,
            "topology_preprocessing_pair_quality": topology_pair_quality,
            "patch_records": patch_records,
            "per_face_sample_records": per_face_sample_records,
            "provenance_records": provenance_records,
            **qualities,
            "timing": {
                "stages": stage,
                "parameterization": self._parameterization_timing_block(stage),
                "end_to_end": {
                    "wall_sec": float(time.perf_counter() - overall_wall),
                    "cpu_sec": float(_worker_cpu_time() - overall_cpu),
                    "cpu_scope": "worker process plus waited-for child processes",
                },
                "gpu": {"available": False, "note": "No GPU backend is used."},
            },
            "memory": {"stage_sample_rss_mb": self._memory_rss_mb()},
            "topology_optimization": {
                "symmetric_dirichlet_energy": {
                    "lscm": energy_lscm,
                    **optcuts_energy,
                    "lscm_to_optcuts_automatic_excess_over_identity_improvement_rate": improvement_rate(
                        energy_lscm,
                        optcuts_energy.get("optcuts_automatic", float("nan")),
                        reference=2.0,
                    ),
                },
                "normalized_internal_seam_length": {
                    "lscm": seam_lscm,
                    **optcuts_seam,
                },
            },
            "residue_footprint_fragmentation": {
                "interaction_residue_source": interaction_source,
                "contact_weight_definition": contact_weight_definition,
                "methods": fragmentation_by_method,
            },
            "optcuts_initialization_ablation": initialization_ablation,
            "optcuts_initialization_pair_quality": initialization_pair_quality,
            "residue_aware_optcuts_ablation": residue_aware_ablation,
            "residue_aware_pair_quality": residue_aware_pair_quality,
            "atlas_trainability": atlas_trainability,
        }
        self._log_thread(f"Finished processing {os.path.basename(pdb_path)}")
        return result

    def _run_operational_optcuts(
        self,
        *,
        pdb_path: str,
        chain_a: str,
        chain_b: str,
        job_metadata: Dict[str, object],
        mesh_a,
        surface_report: Dict[str, object],
        topology_report: Dict[str, object],
        patches,
        prepared_patches,
        preparation: Dict[str, object],
        atoms_a,
        coords_a,
        atoms_b,
        coords_b,
        stage: Dict[str, object],
        overall_wall: float,
        overall_cpu: float,
    ) -> Dict[str, object]:
        """Measure the public single-method pipeline without comparison-only work."""

        method = self.config.resolved_optcuts_variants()[0]
        residue_aware = method in RESIDUE_AWARE_OPTCUTS_METHODS
        source_residue_labels = None
        interaction_weights = None
        objective_weights = None
        interaction_source = None
        contact_weight_definition = None
        if residue_aware:
            stage_wall, stage_cpu = time.perf_counter(), _worker_cpu_time()
            (
                source_residue_labels,
                interaction_weights,
                objective_weights,
                interaction_source,
                contact_weight_definition,
            ) = self._residue_objective(
                job_metadata=job_metadata,
                chain_a=chain_a,
                chain_b=chain_b,
                atoms_a=atoms_a,
                coords_a=coords_a,
                atoms_b=atoms_b,
                coords_b=coords_b,
            )
            stage["contact_weight_preparation"] = self._stage_stats(stage_wall, stage_cpu)
        optimizer = OptCutsUVOptimizer(self.config.optcuts, cancel_event=self.cancel_event)
        optimizer.preflight_binary()
        stage_wall, stage_cpu = time.perf_counter(), _worker_cpu_time()
        outputs, diagnostics = self._run_optcuts(
            prepared_patches,
            initialization="automatic",
            optimizer=optimizer,
            source_residue_labels=source_residue_labels,
            residue_weights=objective_weights,
        )
        stage[method] = self._stage_stats(stage_wall, stage_cpu)

        expected = {
            str(patch.metadata.get("patch_id", "unknown")): face_domain_hash(patch) for patch in prepared_patches
        }
        observed = {str(patch.metadata.get("patch_id", "unknown")): face_domain_hash(patch) for patch in outputs}
        unique_ids = len(observed) == len(outputs)
        exact_ids = {patch_id for patch_id, domain_hash in observed.items() if expected.get(patch_id) == domain_hash}
        complete = bool(unique_ids and set(observed) == set(expected) and exact_ids == set(expected))
        execution_certificates = {
            str(execution.get("patch_id", "unknown")): execution
            for execution in diagnostics.get("executions", [])
            if isinstance(execution, dict)
        }
        certified_ids = {
            patch_id
            for patch_id, execution in execution_certificates.items()
            if bool(
                (execution.get("output_uv_injectivity") or {}).get(
                    "globally_injective",
                    False,
                )
            )
            and bool(
                (execution.get("output_distortion_constraint") or {}).get(
                    "satisfied",
                    False,
                )
            )
        }
        certificate_complete = bool(
            len(execution_certificates) == len(expected)
            and set(execution_certificates) == set(expected)
            and certified_ids == set(expected)
        )
        scientifically_usable = bool(complete and certificate_complete)
        signature = self._comparison_signature(set(expected), prepared_patches)
        elapsed = {
            "wall_sec": float(time.perf_counter() - overall_wall),
            "cpu_sec": float(_worker_cpu_time() - overall_cpu),
            "cpu_scope": "worker process plus waited-for child processes",
        }
        result = {
            "pdb": os.path.basename(pdb_path),
            "input_sha256": str(job_metadata["input_sha256"]),
            "interaction_sha256": job_metadata.get("prolif_sha256") or None,
            "status": "ok" if scientifically_usable else "failed",
            **(
                {}
                if scientifically_usable
                else {
                    "error": (
                        "Operational OptCuts did not preserve every prepared source-face domain."
                        if not complete
                        else "Operational OptCuts did not provide a globally injective, constraint-satisfying "
                        "certificate for every prepared patch."
                    )
                }
            ),
            "chain_selection": {
                "chain_a": chain_a,
                "chain_b": chain_b,
                "mode": job_metadata.get("selection_mode", "configured"),
                "details": job_metadata.get("selection_details", {}),
            },
            **_result_identity_metadata(job_metadata),
            "benchmark_purpose": "performance",
            "execution_profile": "operational_optcuts",
            "operational_method": method,
            "topology_ablation_configured": False,
            "patch_count": int(len(patches)),
            "prepared_patch_count": int(len(prepared_patches)),
            "mesh_stats": {
                "vertex_count": int(len(mesh_a.vertices)),
                "face_count": int(len(mesh_a.faces)),
            },
            "surface_generation": dict(surface_report),
            "topology_extraction": dict(topology_report),
            "topology_preparation": dict(preparation),
            "execution_domain": {
                "complete": complete,
                "scientifically_usable": scientifically_usable,
                "exact_source_face_match": bool(exact_ids == set(expected)),
                "unique_patch_ids": bool(unique_ids),
                "expected_patch_count": int(len(expected)),
                "successful_patch_count": int(len(outputs)),
                "exact_source_face_patch_count": int(len(exact_ids)),
                "signature": signature,
                "rule": ("one automatic OptCuts output for every prepared patch on unchanged source-face geometry"),
            },
            "execution_certificate": {
                "complete": certificate_complete,
                "scientifically_usable": scientifically_usable,
                "expected_patch_count": int(len(expected)),
                "certificate_patch_count": int(len(execution_certificates)),
                "certified_patch_count": int(len(certified_ids)),
                "certified_patch_ids": sorted(certified_ids),
                "rule": (
                    "every exact-domain patch must have an independently recomputed globally injective UV map "
                    "and satisfy the requested raw-scale OptCuts distortion constraint"
                ),
            },
            "comparison_domain": {
                "complete": False,
                "domain_complete": False,
                "metric_complete": False,
                "exact_source_face_match": False,
                "expected_patch_count": 0,
                "common_patch_count": 0,
                "common_patch_ids": [],
                "signature": None,
                "rule": "not evaluated in the operational runtime profile",
            },
            "method_execution": {method: diagnostics},
            "patch_records": [],
            "per_face_sample_records": [],
            "provenance_records": [],
            **(
                {
                    "residue_footprint_fragmentation": {
                        "interaction_residue_source": interaction_source,
                        "contact_weight_definition": contact_weight_definition,
                        "interaction_weights": interaction_weights,
                        "objective_weights": objective_weights,
                    }
                }
                if residue_aware
                else {}
            ),
            "timing": {
                "stages": stage,
                "parameterization": {},
                "end_to_end": elapsed,
                "measurement_scope": (
                    f"structure loading through {method} and packed UV construction; "
                    "quality comparison metrics and ablations are excluded"
                ),
                "gpu": {"available": False, "note": "No GPU backend is used."},
            },
            "memory": {"stage_sample_rss_mb": self._memory_rss_mb()},
        }
        self._log_thread(f"Finished operational timing for {os.path.basename(pdb_path)}")
        return result

    def _failed_single_result(
        self,
        *,
        pdb_path: str,
        chain_a: str,
        chain_b: str,
        job_metadata: Dict[str, object],
        error: str,
        mesh_a,
        surface_report: Dict[str, object],
        topology_report: Dict[str, object],
        preparation: Dict[str, object],
        patch_records: List[Dict[str, object]],
        stage: Dict[str, object],
        overall_wall: float,
        overall_cpu: float,
        raw_patch_count: int,
    ) -> Dict[str, object]:
        """Return partial evidence when every interface component is rejected."""

        return {
            "pdb": os.path.basename(pdb_path),
            "input_sha256": str(job_metadata["input_sha256"]),
            "interaction_sha256": job_metadata.get("prolif_sha256") or None,
            "status": "failed",
            "error": error,
            "chain_selection": {
                "chain_a": chain_a,
                "chain_b": chain_b,
                "mode": job_metadata.get("selection_mode", "configured"),
                "details": job_metadata.get("selection_details", {}),
            },
            **_result_identity_metadata(job_metadata),
            "benchmark_purpose": self.config.benchmark_purpose.strip().lower(),
            "execution_profile": self.config.execution_profile.strip().lower(),
            "topology_ablation_configured": bool(self.config.include_topology_ablation),
            "afdb_complex_confidence": _afdb_complex_confidence(job_metadata),
            "paired_geometry_qc": {
                "contact_cutoff_angstrom": _optional_finite_float(job_metadata.get("paired_contact_cutoff_angstrom")),
                "predicted_contact_count_total": _optional_nonnegative_int(
                    job_metadata.get("paired_predicted_contact_count_total")
                ),
                "contact_recall_fnat": _optional_finite_float(job_metadata.get("paired_contact_recall_fnat")),
                "contact_precision": _optional_finite_float(job_metadata.get("paired_contact_precision")),
                "contact_jaccard": _optional_finite_float(job_metadata.get("paired_contact_jaccard")),
                "experimental_contact_mapping_coverage": _optional_finite_float(
                    job_metadata.get("paired_experimental_contact_mapping_coverage")
                ),
                "interface_residue_a_mapping_coverage": _optional_finite_float(
                    job_metadata.get("paired_interface_residue_a_mapping_coverage")
                ),
                "interface_residue_b_mapping_coverage": _optional_finite_float(
                    job_metadata.get("paired_interface_residue_b_mapping_coverage")
                ),
                "interface_ligand_ca_mapping_coverage": _optional_finite_float(
                    job_metadata.get("paired_interface_ligand_ca_mapping_coverage")
                ),
                "interface_ligand_ca_rmsd_angstrom": _optional_finite_float(
                    job_metadata.get("paired_interface_ligand_ca_rmsd_angstrom")
                ),
                "cross_chain_clash_atom_fraction": _optional_finite_float(
                    job_metadata.get("paired_cross_chain_clash_atom_fraction")
                ),
            },
            "patch_count": int(raw_patch_count),
            "prepared_patch_count": 0,
            "mesh_stats": {
                "vertex_count": int(len(mesh_a.vertices)),
                "face_count": int(len(mesh_a.faces)),
            },
            "surface_generation": dict(surface_report),
            "topology_extraction": dict(topology_report),
            "topology_preparation": dict(preparation),
            "patch_records": patch_records,
            "method_execution": {},
            "comparison_domain": {
                "complete": False,
                "domain_complete": False,
                "metric_complete": False,
                "exact_source_face_match": False,
                "expected_patch_count": 0,
                "common_patch_count": 0,
                "common_patch_ids": [],
                "signature": None,
                "rule": "No component survived the shared preparation domain.",
            },
            "standard_method_pair_quality": {},
            "residue_aware_comparison_domain": {
                "enabled": bool(self.config.optcuts.residue_fragmentation_weight > 0.0),
                "complete": False,
                "domain_complete": False,
                "metric_complete": False,
                "exact_source_face_match": False,
                "expected_patch_count": 0,
                "common_patch_count": 0,
                "common_patch_ids": [],
                "signature": None,
            },
            "topology_preprocessing_pair_quality": {
                "status": "not_evaluated",
                "complete": False,
                "methods": {},
            },
            "timing": {
                "stages": stage,
                "end_to_end": {
                    "wall_sec": float(time.perf_counter() - overall_wall),
                    "cpu_sec": float(_worker_cpu_time() - overall_cpu),
                    "cpu_scope": "worker process plus waited-for child processes",
                },
            },
            "memory": {"stage_sample_rss_mb": self._memory_rss_mb()},
        }

    def _prepare_parameterization_domains(self, patches, parameterizer: Parameterizer):
        prepared = []
        records = []
        for index, patch in enumerate(patches):
            candidate = patch.copy()
            candidate.metadata["patch_id"] = str(patch.metadata.get("patch_id", f"patch_{index:04d}"))
            output, info = parameterizer.prepare_patch(candidate, return_info=True)
            report_info = dict(info)
            provenance = dict(report_info.get("provenance", {}))
            if provenance:
                report_info["provenance"] = {
                    "vertex_count": provenance.get("vertex_count"),
                    "face_count": provenance.get("face_count"),
                    "area": provenance.get("area"),
                    "source_vertex_count": len(provenance.get("source_vertex_ids", [])),
                    "source_face_count": len(provenance.get("source_face_ids", [])),
                    "source_atom_indices": provenance.get("source_atom_indices", []),
                    "history": provenance.get("history", []),
                    "full_mapping_artifact": self.config.provenance_filename,
                }
            records.append({"patch_id": candidate.metadata["patch_id"], **report_info})
            if output is not None:
                prepared.append(candidate)
        return prepared, {
            "attempted": int(len(patches)),
            "success": int(len(prepared)),
            "failed": int(len(patches) - len(prepared)),
            "patches": records,
        }

    def _parameterize_patches(self, patches, method: str, parameterizer: Optional[Parameterizer] = None):
        parameterizer = parameterizer or Parameterizer(config=self.config.parameterization)
        output, wall_times, cpu_times = [], [], []
        failures = []
        for patch in patches:
            self._check_cancelled()
            patch_copy = patch.copy()
            patch_id = str(patch_copy.metadata.get("patch_id", "unknown"))
            wall_start, cpu_start = time.perf_counter(), _worker_cpu_time()
            uv, info = parameterizer.flatten_patch(patch_copy, method=method, return_info=True)
            wall_times.append(float(time.perf_counter() - wall_start))
            cpu_times.append(float(_worker_cpu_time() - cpu_start))
            if uv is None:
                failures.append({"patch_id": patch_id, "reason": info.get("failure_reason", "unknown")})
                continue
            set_uv_layout(patch_copy, uv, key="uv")
            patch_copy.metadata["parameterization_method"] = method
            patch_copy.metadata["parameterization_info"] = info
            output.append(patch_copy)
        return (
            output,
            wall_times,
            cpu_times,
            {
                "attempted": int(len(patches)),
                "success": int(len(output)),
                "failure_count": int(len(failures)),
                "failures": failures,
                "domain_hashes": [face_domain_hash(patch) for patch in output],
            },
        )

    def _run_optcuts(
        self,
        patches,
        *,
        initialization: str,
        optimizer: OptCutsUVOptimizer | None = None,
        source_residue_labels=None,
        residue_weights=None,
    ):
        optimizer = optimizer or OptCutsUVOptimizer(self.config.optcuts, cancel_event=self.cancel_event)
        output = []
        failures = []
        executions = []
        arm_started = time.perf_counter()
        arm_budget_sec = float(self.config.optcuts.timeout_sec)
        invoked_patch_count = 0
        for patch_index, patch in enumerate(patches):
            self._check_cancelled()
            patch_copy = patch.copy()
            patch_id = str(patch_copy.metadata.get("patch_id", "unknown"))
            remaining_sec = arm_budget_sec - (time.perf_counter() - arm_started)
            if remaining_sec <= 0.0:
                failures.extend(
                    {
                        "patch_id": str(item.metadata.get("patch_id", "unknown")),
                        "reason": "OptCuts method-arm time budget was exhausted before this patch.",
                        "failure_type": "arm_budget_exhausted",
                        "invoked": False,
                    }
                    for item in patches[patch_index:]
                )
                break
            invoked_patch_count += 1
            try:
                optimizer.optimize_patches(
                    [patch_copy],
                    initialization=initialization,
                    pack=False,
                    build_report=False,
                    source_residue_labels=source_residue_labels,
                    residue_weights=residue_weights,
                    timeout_sec=remaining_sec,
                )
                output.append(patch_copy)
                executions.append({"patch_id": patch_id, **patch_copy.metadata.get("optcuts_execution", {})})
            except (OSError, RuntimeError, ValueError, subprocess.SubprocessError) as exc:
                self._check_cancelled()
                reason = str(exc)
                failures.append(
                    {
                        "patch_id": patch_id,
                        "reason": reason,
                        "failure_type": "timeout" if "timed out" in reason.lower() else "execution_failure",
                        "invoked": True,
                    }
                )
        optimization_elapsed_sec = float(time.perf_counter() - arm_started)
        if output:
            packed, transforms, report = pack_mesh_charts(
                output,
                key="uv_optcuts",
                gap=self.config.optcuts.patch_gap,
            )
            apply_packed_uv(output, packed, transforms, key="uv_global")
        else:
            report = {"status": "empty", "chart_count": 0}
        return output, {
            "attempted": int(len(patches)),
            "invoked": int(invoked_patch_count),
            "not_invoked": int(len(patches) - invoked_patch_count),
            "success": int(len(output)),
            "failure_count": int(len(failures)),
            "failures": failures,
            "initialization": initialization,
            "fallback_used": False,
            "method_arm_time_budget_sec": arm_budget_sec,
            "method_arm_optimization_wall_sec": optimization_elapsed_sec,
            "method_arm_budget_exhausted": bool(optimization_elapsed_sec >= arm_budget_sec),
            "time_budget_scope": (
                "shared across all external OptCuts solver invocations for this method and structure"
            ),
            "executions": executions,
            "packing": report,
            "domain_hashes": [face_domain_hash(patch) for patch in output],
        }

    @staticmethod
    def _common_patch_ids(method_patches: Dict[str, list]) -> set[str]:
        sets = [
            {str(patch.metadata.get("patch_id", "unknown")) for patch in patches} for patches in method_patches.values()
        ]
        return set.intersection(*sets) if sets else set()

    @staticmethod
    def _filter_patch_ids(patches, patch_ids: set[str]):
        return [patch for patch in patches if str(patch.metadata.get("patch_id", "unknown")) in patch_ids]

    @staticmethod
    def _paired_quality_projection(block: Dict[str, object]) -> Dict[str, object]:
        """Keep the exact-pair endpoints without duplicating atlas diagnostics."""

        keys = (
            "valid_patch_count",
            "scored_face_count",
            "scored_area_3d",
            "domain_hashes",
            "distortion",
            "symmetric_dirichlet",
            "angle_distortion",
            "area_distortion",
            "flip_rate",
            "seam",
            "injectivity",
        )
        return {key: block[key] for key in keys if key in block}

    @staticmethod
    def _residue_aware_ablation_record(
        *,
        comparison_complete: bool,
        residue_fragmentation_weight: float,
        comparisons: Dict[str, object],
    ) -> Dict[str, object]:
        """Expose TopoPPI efficacy values only for a complete exact pair."""

        complete = bool(comparison_complete)
        return {
            "status": "evaluated" if complete else "incomplete_comparison",
            "same_domain": complete,
            "residue_fragmentation_weight": float(residue_fragmentation_weight),
            "residue_weight_definition": ("one plus the number of distinct contacting chain-B residues"),
            "comparisons": dict(comparisons) if complete else {},
            "efficacy_values_available": complete,
            "incomplete_pair_rule": (
                "partial-domain efficacy values are not emitted; arm-level completion, "
                "metric, and injectivity status remain in residue_aware_pair_quality"
            ),
        }

    @staticmethod
    def _comparison_signature(common_ids: set[str], prepared_patches) -> str:
        mapping = {
            str(patch.metadata.get("patch_id", "unknown")): face_domain_hash(patch)
            for patch in prepared_patches
            if str(patch.metadata.get("patch_id", "unknown")) in common_ids
        }
        return hashlib.sha256(json.dumps(mapping, sort_keys=True).encode("utf-8")).hexdigest()

    def _per_face_sample_records(
        self,
        reference_patches,
        method_patches: Dict[str, list],
    ) -> List[Dict[str, object]]:
        """Return a deterministic, source-traceable face sample for audit."""

        method_maps = {
            method: {str(patch.metadata.get("patch_id", "unknown")): patch for patch in patches}
            for method, patches in method_patches.items()
        }
        reference_map = {str(patch.metadata.get("patch_id", "unknown")): patch for patch in reference_patches}
        records: List[Dict[str, object]] = []
        sample_size = int(self.config.per_face_sample_size_per_patch)
        for patch_id in sorted(reference_map):
            reference = reference_map[patch_id]
            face_count = len(reference.faces)
            if face_count == 0:
                continue
            source_ids = np.asarray(
                reference.metadata.get("source_face_ids", np.arange(face_count)),
                dtype=np.int64,
            )
            ordered = np.argsort(source_ids, kind="stable")
            if len(ordered) > sample_size:
                positions = np.linspace(0, len(ordered) - 1, num=sample_size, dtype=np.int64)
                selected = ordered[positions]
            else:
                selected = ordered

            method_samples = {}
            for method, patch_map in method_maps.items():
                patch = patch_map.get(patch_id)
                if patch is None:
                    continue
                uv_key = "uv_optcuts" if _uses_optcuts_uv(method) else "uv"
                uv = as_corner_uv(patch, key=uv_key)
                distortion, _weights = UVAtlasMetrics.distortion_samples(patch, uv)
                symmetric_dirichlet, _weights = UVAtlasMetrics.symmetric_dirichlet_samples(patch, uv)
                angle, _weights = UVAtlasMetrics.angle_distortion_samples(patch, uv)
                area, _weights = UVAtlasMetrics.area_distortion_samples(patch, uv)
                flipped, _weights = UVAtlasMetrics.flip_samples(patch, uv)
                method_source_ids = np.asarray(
                    patch.metadata.get("source_face_ids", np.arange(len(patch.faces))),
                    dtype=np.int64,
                )
                source_to_index = {int(source_id): index for index, source_id in enumerate(method_source_ids)}
                method_samples[method] = (
                    distortion,
                    symmetric_dirichlet,
                    angle,
                    area,
                    flipped,
                    source_to_index,
                )

            areas = np.asarray(reference.area_faces, dtype=np.float64)
            for sampling_rank, face_index in enumerate(selected):
                row: Dict[str, object] = {
                    "patch_id": patch_id,
                    "source_face_id": int(source_ids[face_index]),
                    "face_area_3d": float(areas[face_index]),
                    "sampling_rank": int(sampling_rank),
                    "sampling_rule": "source_face_id_sorted_evenly_spaced_v1",
                }
                source_face_id = int(source_ids[face_index])
                for method, samples in method_samples.items():
                    distortion, symmetric_dirichlet, angle, area, flipped, source_to_index = samples
                    method_face_index = source_to_index.get(source_face_id)
                    if method_face_index is None:
                        continue
                    row[f"{method}_distortion"] = float(distortion[method_face_index])
                    row[f"{method}_symmetric_dirichlet"] = float(symmetric_dirichlet[method_face_index])
                    row[f"{method}_angle_distortion_rad"] = float(angle[method_face_index])
                    row[f"{method}_area_distortion"] = float(area[method_face_index])
                    row[f"{method}_flipped_after_global_reflection"] = int(bool(flipped[method_face_index]))
                records.append(row)
        return records

    @staticmethod
    def _provenance_records(patches) -> List[Dict[str, object]]:
        records = []
        for patch in patches:
            patch_id = str(patch.metadata.get("patch_id", "unknown"))
            source_faces = np.asarray(
                patch.metadata.get("source_face_ids", np.arange(len(patch.faces))),
                dtype=np.int64,
            )
            source_vertices = np.asarray(
                patch.metadata.get("source_vertex_ids", np.arange(len(patch.vertices))),
                dtype=np.int64,
            )
            source_atoms = np.asarray(
                patch.metadata.get("source_atom_indices", np.full(len(patch.vertices), -1)),
                dtype=np.int64,
            )
            records.extend(
                {
                    "patch_id": patch_id,
                    "entity": "face",
                    "final_index": int(index),
                    "source_id": int(source_id),
                    "source_atom_index": "",
                }
                for index, source_id in enumerate(source_faces)
            )
            records.extend(
                {
                    "patch_id": patch_id,
                    "entity": "vertex",
                    "final_index": int(index),
                    "source_id": int(source_id),
                    "source_atom_index": int(source_atoms[index]) if index < len(source_atoms) else -1,
                }
                for index, source_id in enumerate(source_vertices)
            )
        return records

    def _patch_retention_records(
        self,
        prepared_patches,
        atoms_a,
        coords_a,
        atoms_b,
        coords_b,
        *,
        extracted_patches=None,
        topology_components=None,
        preparation=None,
        job_metadata: Dict[str, object],
        contact_pairs_by_atom=None,
    ):
        records = []
        extracted_patches = list(extracted_patches or [])
        topology_components = list(topology_components or [])
        preparation_records = {
            str(item.get("patch_id", "unknown")): dict(item) for item in dict(preparation or {}).get("patches", [])
        }
        prepared_map = {str(patch.metadata.get("patch_id", "unknown")): patch for patch in prepared_patches}

        def _empty_summary() -> Dict[str, object]:
            return {
                "face_count": 0,
                "vertex_count": 0,
                "area": 0.0,
                "source_face_ids": [],
                "source_vertex_ids": [],
                "source_atom_indices": [],
            }

        descriptors = []
        for extracted in extracted_patches:
            patch_id = str(extracted.metadata.get("patch_id", "unknown"))
            prepared = prepared_map.get(patch_id)
            raw_summary = dict(extracted.metadata.get("topology_component_before") or provenance_summary(extracted))
            topology_summary = provenance_summary(extracted)
            final_summary = provenance_summary(prepared) if prepared is not None else _empty_summary()
            info = (
                dict(prepared.metadata.get("parameterization_preparation", {}))
                if prepared is not None
                else preparation_records.get(patch_id, {})
            )
            descriptors.append(
                {
                    "patch_id": patch_id,
                    "source_face_hash": face_domain_hash(prepared if prepared is not None else extracted),
                    "status": "prepared" if prepared is not None else "rejected",
                    "rejection_stage": None if prepared is not None else "parameterization_topology_gate",
                    "failure_reason": info.get("failure_reason"),
                    "parameterization_attempted": True,
                    "raw": raw_summary,
                    "topology": topology_summary,
                    "final": final_summary,
                }
            )
        for component in topology_components:
            if component.get("status") == "accepted":
                continue
            descriptors.append(
                {
                    "patch_id": str(component.get("patch_id") or f"component_{len(descriptors):04d}"),
                    "source_face_hash": "",
                    "status": "rejected",
                    "rejection_stage": "topology_extraction",
                    "failure_reason": component.get("reason") or "component_dropped",
                    "parameterization_attempted": False,
                    "raw": dict(component.get("before_sanitation") or _empty_summary()),
                    "topology": dict(component.get("after_sanitation") or _empty_summary()),
                    "final": _empty_summary(),
                }
            )

        hotspot_tokens = self._parse_manifest_list(job_metadata.get("hotspot_residues_a"))
        confidence_metric = str(job_metadata.get("confidence_metric") or "").strip().lower()
        confidence_threshold = self._finite_float_or_default(
            job_metadata.get("confidence_threshold"),
            70.0,
        )
        declared_interactions, interaction_provenance = self._load_declared_interactions(
            job_metadata.get("prolif_file"),
            expected_chain_a=str(job_metadata.get("chain_a") or ""),
            expected_chain_b=str(job_metadata.get("chain_b") or ""),
            expected_input_sha256=str(job_metadata.get("input_sha256") or ""),
            known_file_sha256=str(job_metadata.get("prolif_sha256") or ""),
            require_bindings=bool(self.config.formal_mode),
            atoms_a=atoms_a,
            atoms_b=atoms_b,
        )
        all_raw_atoms = {
            int(value) for descriptor in descriptors for value in descriptor["raw"].get("source_atom_indices", [])
        }
        if contact_pairs_by_atom is None:
            tree_b = cKDTree(coords_b) if len(coords_b) else None
            contact_pairs_by_atom = self._geometric_contact_pairs_by_atom(
                all_raw_atoms,
                atoms_a,
                coords_a,
                atoms_b,
                tree_b,
            )
        for descriptor in descriptors:
            raw_summary = descriptor["raw"]
            topology_summary = descriptor["topology"]
            final_summary = descriptor["final"]
            raw_atoms = set(int(value) for value in raw_summary.get("source_atom_indices", []))
            topology_atoms = set(int(value) for value in topology_summary.get("source_atom_indices", []))
            final_atoms = set(int(value) for value in final_summary.get("source_atom_indices", []))

            raw_residues = self._residue_labels(raw_atoms, atoms_a)
            topology_residues = self._residue_labels(topology_atoms, atoms_a)
            final_residues = self._residue_labels(final_atoms, atoms_a)
            raw_sequences = self._residue_sequences(raw_atoms, atoms_a)
            topology_sequences = self._residue_sequences(topology_atoms, atoms_a)
            final_sequences = self._residue_sequences(final_atoms, atoms_a)
            raw_confidence = self._atom_confidence_values(raw_atoms, atoms_a, confidence_metric)
            topology_confidence = self._atom_confidence_values(
                topology_atoms,
                atoms_a,
                confidence_metric,
            )
            final_confidence = self._atom_confidence_values(final_atoms, atoms_a, confidence_metric)

            raw_contacts = self._merge_contact_pairs(raw_atoms, contact_pairs_by_atom)
            topology_contacts = self._merge_contact_pairs(topology_atoms, contact_pairs_by_atom)
            final_contacts = self._merge_contact_pairs(final_atoms, contact_pairs_by_atom)

            hotspots_on_raw = hotspot_tokens & self._hotspot_aliases(raw_atoms, atoms_a)
            hotspots_after = hotspots_on_raw & self._hotspot_aliases(final_atoms, atoms_a)
            hotspots_after_topology = hotspots_on_raw & self._hotspot_aliases(topology_atoms, atoms_a)
            interactions_on_raw = {
                identifier for identifier, sequence in declared_interactions.items() if sequence in raw_sequences
            }
            interactions_after_topology = {
                identifier
                for identifier in interactions_on_raw
                if declared_interactions[identifier] in topology_sequences
            }
            interactions_after = {
                identifier for identifier in interactions_on_raw if declared_interactions[identifier] in final_sequences
            }

            raw_face_count = int(raw_summary.get("face_count", 0))
            raw_vertex_count = int(raw_summary.get("vertex_count", 0))
            raw_source_vertices = set(int(value) for value in raw_summary.get("source_vertex_ids", []))
            raw_area = float(raw_summary.get("area", 0.0))
            topology_face_count = int(topology_summary.get("face_count", 0))
            topology_vertex_count = int(topology_summary.get("vertex_count", 0))
            topology_source_vertices = set(int(value) for value in topology_summary.get("source_vertex_ids", []))
            topology_area = float(topology_summary.get("area", 0.0))
            final_face_count = int(final_summary.get("face_count", 0))
            final_vertex_count = int(final_summary.get("vertex_count", 0))
            final_source_vertices = set(int(value) for value in final_summary.get("source_vertex_ids", []))
            final_area = float(final_summary.get("area", 0.0))
            parameterization_attempted = bool(descriptor["parameterization_attempted"])
            records.append(
                {
                    "patch_id": descriptor["patch_id"],
                    "source_face_hash": descriptor["source_face_hash"],
                    "retention_status": descriptor["status"],
                    "rejection_stage": descriptor["rejection_stage"],
                    "failure_reason": descriptor["failure_reason"],
                    "retention_denominator": "raw_interface_component_before_sanitation",
                    "face_count_before": raw_face_count,
                    "face_count_after_topology_sanitation": topology_face_count,
                    "face_count_after": final_face_count,
                    "materialized_vertex_count_before": raw_vertex_count,
                    "materialized_vertex_count_after_topology_sanitation": topology_vertex_count,
                    "materialized_vertex_count_after": final_vertex_count,
                    "source_vertex_count_before": int(len(raw_source_vertices)),
                    "source_vertex_count_after_topology_sanitation": int(len(topology_source_vertices)),
                    "source_vertex_count_after": int(len(final_source_vertices)),
                    "area_before": raw_area,
                    "area_after_topology_sanitation": topology_area,
                    "area_after": final_area,
                    "source_atom_count_before": int(len(raw_atoms)),
                    "source_atom_count_after_topology_sanitation": int(len(topology_atoms)),
                    "source_atom_count_after": int(len(final_atoms)),
                    "face_retention_ratio": float(final_face_count / raw_face_count)
                    if raw_face_count
                    else float("nan"),
                    "materialized_vertex_count_ratio": float(final_vertex_count / raw_vertex_count)
                    if raw_vertex_count
                    else float("nan"),
                    "source_vertex_retention_ratio": float(len(final_source_vertices) / len(raw_source_vertices))
                    if raw_source_vertices
                    else float("nan"),
                    "area_retention_ratio": float(final_area / raw_area) if raw_area > 0.0 else float("nan"),
                    "source_atom_retention_ratio": float(len(final_atoms) / len(raw_atoms))
                    if raw_atoms
                    else float("nan"),
                    "topology_face_retention_ratio": float(topology_face_count / raw_face_count)
                    if raw_face_count
                    else float("nan"),
                    "topology_materialized_vertex_count_ratio": float(topology_vertex_count / raw_vertex_count)
                    if raw_vertex_count
                    else float("nan"),
                    "topology_source_vertex_retention_ratio": float(
                        len(topology_source_vertices) / len(raw_source_vertices)
                    )
                    if raw_source_vertices
                    else float("nan"),
                    "topology_area_retention_ratio": float(topology_area / raw_area)
                    if raw_area > 0.0
                    else float("nan"),
                    "topology_source_atom_retention_ratio": float(len(topology_atoms) / len(raw_atoms))
                    if raw_atoms
                    else float("nan"),
                    "parameterization_face_retention_ratio": float(final_face_count / topology_face_count)
                    if parameterization_attempted and topology_face_count
                    else float("nan"),
                    "parameterization_materialized_vertex_count_ratio": float(
                        final_vertex_count / topology_vertex_count
                    )
                    if parameterization_attempted and topology_vertex_count
                    else float("nan"),
                    "parameterization_source_vertex_retention_ratio": float(
                        len(final_source_vertices) / len(topology_source_vertices)
                    )
                    if parameterization_attempted and topology_source_vertices
                    else float("nan"),
                    "parameterization_area_retention_ratio": float(final_area / topology_area)
                    if parameterization_attempted and topology_area > 0.0
                    else float("nan"),
                    "parameterization_source_atom_retention_ratio": float(len(final_atoms) / len(topology_atoms))
                    if parameterization_attempted and topology_atoms
                    else float("nan"),
                    "residue_count_before": int(len(raw_residues)),
                    "residue_count_after_topology_sanitation": int(len(topology_residues)),
                    "residue_count_after": int(len(final_residues)),
                    "topology_residue_retention_ratio": float(len(topology_residues) / len(raw_residues))
                    if raw_residues
                    else float("nan"),
                    "parameterization_residue_retention_ratio": float(len(final_residues) / len(topology_residues))
                    if topology_residues
                    else float("nan"),
                    "residue_retention_ratio": float(len(final_residues) / len(raw_residues))
                    if raw_residues
                    else float("nan"),
                    "removed_residues": sorted(raw_residues - final_residues),
                    "geometric_contact_definition": (
                        f"A/B heavy-atom residue pair distance <= {self.config.contact_distance_angstrom:g} Angstrom"
                    ),
                    "geometric_contact_pair_count_before": int(len(raw_contacts)),
                    "geometric_contact_pair_count_after_topology_sanitation": int(len(topology_contacts)),
                    "geometric_contact_pair_count_after": int(len(final_contacts)),
                    "geometric_contact_pair_retention_ratio": float(len(final_contacts) / len(raw_contacts))
                    if raw_contacts
                    else float("nan"),
                    "topology_geometric_contact_pair_retention_ratio": float(len(topology_contacts) / len(raw_contacts))
                    if raw_contacts
                    else float("nan"),
                    "parameterization_geometric_contact_pair_retention_ratio": float(
                        len(final_contacts) / len(topology_contacts)
                    )
                    if topology_contacts
                    else float("nan"),
                    "removed_geometric_contact_pairs": sorted(
                        "--".join(pair) for pair in raw_contacts - final_contacts
                    ),
                    "declared_hotspot_source": "benchmark_manifest.hotspot_residues_a",
                    "declared_hotspot_count_global": int(len(hotspot_tokens)),
                    "declared_hotspot_count_on_patch_before": int(len(hotspots_on_raw)),
                    "declared_hotspot_count_after_topology_sanitation": int(len(hotspots_after_topology)),
                    "declared_hotspot_count_after": int(len(hotspots_after)),
                    "declared_hotspot_retention_ratio": float(len(hotspots_after) / len(hotspots_on_raw))
                    if hotspots_on_raw
                    else float("nan"),
                    "topology_declared_hotspot_retention_ratio": float(
                        len(hotspots_after_topology) / len(hotspots_on_raw)
                    )
                    if hotspots_on_raw
                    else float("nan"),
                    "parameterization_declared_hotspot_retention_ratio": float(
                        len(hotspots_after) / len(hotspots_after_topology)
                    )
                    if hotspots_after_topology
                    else float("nan"),
                    "removed_declared_hotspots": sorted(hotspots_on_raw - hotspots_after),
                    "declared_interaction_provenance": interaction_provenance,
                    "declared_interaction_count_on_patch_before": int(len(interactions_on_raw)),
                    "declared_interaction_count_after_topology_sanitation": int(len(interactions_after_topology)),
                    "declared_interaction_count_after": int(len(interactions_after)),
                    "declared_interaction_retention_ratio": float(len(interactions_after) / len(interactions_on_raw))
                    if interactions_on_raw
                    else float("nan"),
                    "topology_declared_interaction_retention_ratio": float(
                        len(interactions_after_topology) / len(interactions_on_raw)
                    )
                    if interactions_on_raw
                    else float("nan"),
                    "parameterization_declared_interaction_retention_ratio": float(
                        len(interactions_after) / len(interactions_after_topology)
                    )
                    if interactions_after_topology
                    else float("nan"),
                    "removed_declared_interactions": sorted(interactions_on_raw - interactions_after),
                    "confidence_metric": confidence_metric or "not_declared",
                    "confidence_source": job_metadata.get("confidence_source") or None,
                    "confidence_threshold": confidence_threshold if confidence_metric else float("nan"),
                    "confidence_atom_count_before": int(len(raw_confidence)),
                    "confidence_atom_count_after_topology_sanitation": int(len(topology_confidence)),
                    "confidence_atom_count_after": int(len(final_confidence)),
                    "confidence_atom_retention_ratio": float(len(final_confidence) / len(raw_confidence))
                    if len(raw_confidence)
                    else float("nan"),
                    "topology_confidence_atom_retention_ratio": float(len(topology_confidence) / len(raw_confidence))
                    if len(raw_confidence)
                    else float("nan"),
                    "parameterization_confidence_atom_retention_ratio": float(
                        len(final_confidence) / len(topology_confidence)
                    )
                    if len(topology_confidence)
                    else float("nan"),
                    "confidence_mean_before": self._mean_or_nan(raw_confidence),
                    "confidence_mean_after_topology_sanitation": self._mean_or_nan(topology_confidence),
                    "confidence_mean_after": self._mean_or_nan(final_confidence),
                    "low_confidence_atom_fraction_before": self._fraction_below(
                        raw_confidence,
                        confidence_threshold,
                    ),
                    "low_confidence_atom_fraction_after_topology_sanitation": self._fraction_below(
                        topology_confidence,
                        confidence_threshold,
                    ),
                    "low_confidence_atom_fraction_after": self._fraction_below(
                        final_confidence,
                        confidence_threshold,
                    ),
                }
            )
        return records

    @staticmethod
    def _finite_float_or_default(value: object, default: float) -> float:
        try:
            number = float(value)
        except (TypeError, ValueError):
            return float(default)
        return number if np.isfinite(number) else float(default)

    @staticmethod
    def _atom_confidence_values(atom_indices, atoms, metric: str) -> np.ndarray:
        """Read AlphaFold-style pLDDT stored in the atom B-factor field.

        Experimental B factors are deliberately not interpreted as confidence;
        the manifest must explicitly declare a pLDDT metric.
        """

        normalized = metric.replace("-", "_")
        if normalized not in {"plddt", "plddt_bfactor", "b_factor_plddt"}:
            return np.empty(0, dtype=np.float64)
        values = []
        for index in atom_indices:
            if not 0 <= int(index) < len(atoms):
                continue
            try:
                value = float(atoms[int(index)].get_bfactor())
            except (AttributeError, TypeError, ValueError):
                continue
            if np.isfinite(value):
                values.append(value)
        return np.asarray(values, dtype=np.float64)

    @staticmethod
    def _mean_or_nan(values: np.ndarray) -> float:
        return float(np.mean(values)) if len(values) else float("nan")

    @staticmethod
    def _fraction_below(values: np.ndarray, threshold: float) -> float:
        return float(np.mean(values < threshold)) if len(values) else float("nan")

    _atom_residue_label = staticmethod(atom_residue_label)

    @staticmethod
    def _atom_residue_sequence(atom) -> str:
        residue = atom.get_parent()
        insertion = str(residue.id[2]).strip()
        return f"{residue.id[1]}{insertion}".upper()

    def _residue_labels(self, atom_indices, atoms) -> set[str]:
        return {self._atom_residue_label(atoms[index]) for index in atom_indices if 0 <= index < len(atoms)}

    def _residue_sequences(self, atom_indices, atoms) -> set[str]:
        return {self._atom_residue_sequence(atoms[index]) for index in atom_indices if 0 <= index < len(atoms)}

    @staticmethod
    def _parse_manifest_list(value: object) -> set[str]:
        if value is None:
            return set()
        if isinstance(value, (list, tuple, set)):
            raw_values = list(value)
        else:
            text = str(value).strip()
            if not text:
                return set()
            try:
                parsed = json.loads(text)
                raw_values = parsed if isinstance(parsed, list) else [parsed]
            except json.JSONDecodeError:
                raw_values = re.split(r"[;,|]", text)
        return {str(item).strip().upper().replace(" ", "") for item in raw_values if str(item).strip()}

    def _hotspot_aliases(self, atom_indices, atoms) -> set[str]:
        aliases: set[str] = set()
        for atom_index in atom_indices:
            if not 0 <= atom_index < len(atoms):
                continue
            atom = atoms[atom_index]
            residue = atom.get_parent()
            chain = residue.get_parent()
            sequence = self._atom_residue_sequence(atom)
            residue_name = str(residue.get_resname()).upper()
            aliases.update(
                {
                    sequence,
                    f"{chain.id}:{sequence}".upper(),
                    f"{residue_name}:{sequence}",
                    f"{chain.id}:{residue_name}:{sequence}".upper(),
                }
            )
        return aliases

    def _load_declared_interactions(
        self,
        path_value: object,
        *,
        expected_chain_a: str,
        expected_chain_b: str,
        expected_input_sha256: str = "",
        known_file_sha256: str = "",
        require_bindings: bool = False,
        atoms_a=None,
        atoms_b=None,
    ) -> tuple[Dict[str, str], Dict[str, object]]:
        path = str(path_value or "").strip()
        if not path:
            return {}, {"status": "not_declared"}
        with open(path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
        if not isinstance(payload, dict):
            raise ValueError(f"Interaction JSON root must be an object: {path}")
        engine = str(payload.get("engine") or "").strip().lower()
        if engine and engine != "prolif":
            raise ValueError(f"Interaction JSON engine must be 'prolif', got {engine!r}.")
        file_chain_a = str(payload.get("chain_a") or "").strip()
        file_chain_b = str(payload.get("chain_b") or "").strip()
        declared_source_sha256 = str(payload.get("input_sha256") or payload.get("source_sha256") or "").strip().lower()
        if require_bindings and (not file_chain_a or not file_chain_b or not declared_source_sha256):
            raise ValueError("Formal interaction JSON must declare chain_a, chain_b, and source_sha256/input_sha256.")
        if file_chain_a and expected_chain_a and file_chain_a != expected_chain_a:
            raise ValueError(f"Interaction file chain_a={file_chain_a} does not match {expected_chain_a}.")
        if file_chain_b and expected_chain_b and file_chain_b != expected_chain_b:
            raise ValueError(f"Interaction file chain_b={file_chain_b} does not match {expected_chain_b}.")
        if (
            declared_source_sha256
            and expected_input_sha256
            and declared_source_sha256 != expected_input_sha256.strip().lower()
        ):
            raise ValueError(
                "Interaction file source checksum does not match the benchmark structure: "
                f"expected {expected_input_sha256}, got {declared_source_sha256}."
            )
        raw_records = payload.get("interactions", [])
        valid_a = {self._atom_residue_sequence(atom) for atom in atoms_a} if atoms_a is not None else None
        valid_b = {self._atom_residue_sequence(atom) for atom in atoms_b} if atoms_b is not None else None
        interactions: Dict[str, str] = {}
        for index, item in enumerate(raw_records):
            if not isinstance(item, dict):
                continue
            sequence_a = residue_sequence_token(item.get("res_a_seq"))
            if sequence_a is None:
                continue
            sequence_b = residue_sequence_token(item.get("res_b_seq"))
            if sequence_b is None:
                continue
            if valid_a is not None and sequence_a not in valid_a:
                continue
            if valid_b is not None and sequence_b not in valid_b:
                continue
            interaction_type = str(item.get("interaction") or "unspecified").strip() or "unspecified"
            identifier = f"{sequence_a}--{sequence_b}--{interaction_type}--record{index}"
            interactions[identifier] = sequence_a
        return interactions, {
            "status": "loaded",
            "path": os.path.abspath(path),
            "sha256": known_file_sha256 or sha256_file(path),
            "record_count": int(len(interactions)),
            "declared_source_sha256": declared_source_sha256 or None,
            "source_checksum_matches": bool(
                not declared_source_sha256
                or not expected_input_sha256
                or declared_source_sha256 == expected_input_sha256.strip().lower()
            ),
        }

    def _residue_objective(
        self,
        *,
        job_metadata: Dict[str, object],
        chain_a: str,
        chain_b: str,
        atoms_a,
        coords_a,
        atoms_b,
        coords_b,
    ) -> tuple[List[str], Dict[str, float], Dict[str, float], str, str]:
        interaction_path = str(job_metadata.get("prolif_file") or "").strip()
        if interaction_path:
            cached_partner_map = job_metadata.get("interaction_partner_map")
            partner_map = (
                cached_partner_map
                if isinstance(cached_partner_map, dict)
                else load_prolif_partner_map(
                    interaction_path,
                    atoms_a,
                    atoms_b,
                    expected_chain_a=chain_a,
                    expected_chain_b=chain_b,
                    expected_source_sha256=str(job_metadata.get("input_sha256") or ""),
                )
            )
            if not partner_map:
                raise ValueError("ProLIF did not yield any interaction residue pairs for the selected chains.")
            source = "prolif"
            definition = "distinct Chain-B residues paired with each Chain-A residue in ProLIF records"
        else:
            if self.config.formal_mode:
                raise ValueError("Formal residue-aware benchmarks require declared ProLIF evidence.")
            partner_map = geometric_contact_partner_map(
                coords_a,
                atoms_a,
                coords_b,
                atoms_b,
                distance_cutoff=float(self.config.contact_distance_angstrom),
            )
            source = "geometric_fallback"
            definition = (
                "distinct Chain-B residues with any heavy-atom pair at distance <= "
                f"{self.config.contact_distance_angstrom:g} Angstrom (explicit geometric fallback)"
            )
        source_labels = [str(value) for value in source_atom_residue_labels(atoms_a)]
        interaction_weights = contact_partner_degrees(partner_map)
        objective_weights = residue_aware_residue_weights(source_labels, interaction_weights)
        return source_labels, interaction_weights, objective_weights, source, definition

    def _geometric_contact_pairs_by_atom(self, atom_indices, atoms_a, coords_a, atoms_b, tree_b):
        pairs_by_atom = {int(index): set() for index in atom_indices}
        if tree_b is None:
            return pairs_by_atom
        for atom_index in atom_indices:
            if not 0 <= atom_index < len(coords_a):
                continue
            residue_a = self._atom_residue_label(atoms_a[atom_index])
            for partner_index in tree_b.query_ball_point(
                coords_a[atom_index],
                r=float(self.config.contact_distance_angstrom),
            ):
                pairs_by_atom[atom_index].add((residue_a, self._atom_residue_label(atoms_b[partner_index])))
        return pairs_by_atom

    @staticmethod
    def _merge_contact_pairs(atom_indices, pairs_by_atom):
        pairs = set()
        for atom_index in atom_indices:
            pairs.update(pairs_by_atom.get(atom_index, ()))
        return pairs

    def _validate_manifest_record(
        self,
        record: Dict[str, object],
        *,
        actual_sha256: str,
    ) -> None:
        expected_sha256 = (
            str(record.get("input_sha256") or record.get("sha256") or record.get("file_sha256") or "").strip().lower()
        )
        if self.config.formal_mode and not expected_sha256:
            raise ValueError("Formal manifest record is missing input_sha256")
        if expected_sha256 and actual_sha256.lower() != expected_sha256:
            raise ValueError(f"Input checksum mismatch: expected {expected_sha256}, got {actual_sha256}")
        if not self.config.formal_mode:
            return

        record_id = str(record.get("record_id") or record.get("id") or "").strip()
        if not record_id:
            raise ValueError("Formal manifest record is missing record_id/id")

        cluster_id = str(record.get("cluster_id") or record.get("interface_cluster_id") or "").strip()
        if not cluster_id:
            raise ValueError("Formal manifest record is missing cluster_id/interface_cluster_id")
        missing_dependence_metadata = [
            name
            for name in ("family_id", "sequence_cluster_a", "sequence_cluster_b")
            if not str(record.get(name) or "").strip()
        ]
        if missing_dependence_metadata:
            raise ValueError(
                "Formal manifest record is missing inferential-dependence metadata: "
                + ", ".join(missing_dependence_metadata)
            )
        required_metadata = {
            "dataset_source": record.get("dataset_source") or record.get("source"),
            "source_accession": record.get("source_accession") or record.get("accession"),
            "license_or_terms": record.get("license_or_terms") or record.get("license"),
            "structure_type": record.get("structure_type"),
        }
        missing = [name for name, value in required_metadata.items() if not str(value or "").strip()]
        if missing:
            raise ValueError("Formal manifest record is missing: " + ", ".join(missing))
        raw_analysis_split = str(record.get("analysis_split") or "").strip()
        if not raw_analysis_split:
            raise ValueError("Formal manifest record is missing analysis_split")
        analysis_split = raw_analysis_split.lower()
        if analysis_split not in {"development", "test", "exploratory"}:
            raise ValueError("analysis_split must be development, test, or exploratory")
        missing_split_metadata = [
            name
            for name in ("analysis_split_component_id", "analysis_split_basis")
            if not str(record.get(name) or "").strip()
        ]
        if missing_split_metadata:
            raise ValueError(
                "Formal manifest record is missing split-dependence metadata: " + ", ".join(missing_split_metadata)
            )
        structure_type = str(required_metadata["structure_type"]).strip().lower()
        if structure_type not in FORMAL_STRUCTURE_TYPES:
            raise ValueError(
                "Formal manifest record has unsupported structure_type; expected one of "
                + ", ".join(sorted(FORMAL_STRUCTURE_TYPES))
            )
        if structure_type in PREDICTED_STRUCTURE_TYPES:
            missing_prediction_dependencies = [
                name for name in INFERENCE_DEPENDENCY_FIELDS if not str(record.get(name) or "").strip()
            ]
            if missing_prediction_dependencies:
                raise ValueError(
                    "Formal predicted record is missing prediction-dependence metadata: "
                    + ", ".join(missing_prediction_dependencies)
                )
            inference_a = str(record["inference_sequence_cluster_a"]).strip()
            inference_b = str(record["inference_sequence_cluster_b"]).strip()
            if str(record["inference_family_id"]).strip() != inference_family_id(inference_a, inference_b):
                raise ValueError("Formal predicted record has an invalid inference_family_id")
            if str(record["inference_dependency_basis"]).strip() != INFERENCE_DEPENDENCY_BASIS:
                raise ValueError("Formal predicted record has an invalid inference_dependency_basis")
        paired_record_id = str(record.get("paired_record_id") or "").strip()
        if structure_type == "afdb_monomer_replacement" and not paired_record_id:
            raise ValueError("AFDB monomer-replacement record is missing paired_record_id")
        if paired_record_id:
            if (
                structure_type in PREDICTED_STRUCTURE_TYPES
                and not str(record.get("paired_experimental_record_id") or "").strip()
            ):
                raise ValueError("Paired predicted record is missing paired_experimental_record_id")
            stratum = str(record.get("paired_geometry_stratum") or "").strip()
            allowed_strata = {"high_fidelity", "moderate_fidelity", "geometry_stress_test"}
            if stratum not in allowed_strata:
                raise ValueError(
                    "Paired record has invalid paired_geometry_stratum; expected one of "
                    + ", ".join(sorted(allowed_strata))
                )
            contact_cutoff = _optional_finite_float(record.get("paired_contact_cutoff_angstrom"))
            if contact_cutoff is None or contact_cutoff <= 0.0:
                raise ValueError(
                    "Paired record has invalid paired_contact_cutoff_angstrom; expected a positive finite value"
                )
            try:
                predicted_contact_count = int(record["paired_predicted_contact_count_total"])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(
                    "Paired record has invalid paired_predicted_contact_count_total; expected a non-negative integer"
                ) from exc
            if predicted_contact_count < 0:
                raise ValueError(
                    "Paired record has invalid paired_predicted_contact_count_total; expected a non-negative integer"
                )
            for side in ("a", "b"):
                count_name = f"paired_alignment_{side}_optimal_correspondence_count"
                try:
                    correspondence_count = int(record[count_name])
                except (KeyError, TypeError, ValueError):
                    correspondence_count = 0
                if correspondence_count < 1:
                    raise ValueError(f"Paired record has invalid {count_name}; expected a positive integer")
                consensus_name = f"paired_alignment_{side}_selected_pair_consensus_fraction"
                consensus = _optional_finite_float(record.get(consensus_name))
                if consensus is None or not 0.0 <= consensus <= 1.0:
                    raise ValueError(f"Paired record has invalid {consensus_name}; expected a value in [0, 1]")
            optional_bounded_fields = (
                "paired_contact_recall_fnat",
                "paired_contact_precision",
                "paired_contact_jaccard",
                "paired_experimental_contact_mapping_coverage",
                "paired_interface_residue_a_mapping_coverage",
                "paired_interface_residue_b_mapping_coverage",
                "paired_interface_ligand_ca_mapping_coverage",
            )
            for name in optional_bounded_fields:
                value = _optional_finite_float(record.get(name))
                if value is not None and not 0.0 <= value <= 1.0:
                    raise ValueError(f"Paired record has invalid {name}; expected a value in [0, 1]")
            clash_fraction = _optional_finite_float(record.get("paired_cross_chain_clash_atom_fraction"))
            if clash_fraction is None or not 0.0 <= clash_fraction <= 1.0:
                raise ValueError(
                    "Paired record has invalid paired_cross_chain_clash_atom_fraction; expected a value in [0, 1]"
                )
            rmsd = _optional_finite_float(record.get("paired_interface_ligand_ca_rmsd_angstrom"))
            if rmsd is not None and rmsd < 0.0:
                raise ValueError(
                    "Paired record has invalid paired_interface_ligand_ca_rmsd_angstrom; "
                    "expected a non-negative finite value"
                )
        if structure_type not in PREDICTED_STRUCTURE_TYPES:
            return
        missing_confidence = [
            name for name in ("confidence_metric", "confidence_source") if not str(record.get(name) or "").strip()
        ]
        if missing_confidence:
            raise ValueError("Predicted-structure manifest record is missing: " + ", ".join(missing_confidence))
        confidence_metric = str(record["confidence_metric"]).strip().lower().replace("-", "_")
        if confidence_metric not in {"plddt", "plddt_bfactor", "b_factor_plddt"}:
            raise ValueError(
                "Unsupported predicted-structure confidence_metric; use plddt/plddt_bfactor/b_factor_plddt"
            )
        threshold = self._finite_float_or_default(record.get("confidence_threshold"), 70.0)
        if not 0.0 <= threshold <= 100.0:
            raise ValueError("confidence_threshold must be between 0 and 100 for pLDDT")
        if structure_type in {"afdb", "afdb_monomer_replacement"}:
            if not str(record.get("afdb_model_id") or record.get("model_id") or "").strip():
                raise ValueError("AFDB manifest record is missing afdb_model_id")
        if structure_type == "afdb":
            for name in ("afdb_iptm", "afdb_ipsae"):
                value = self._finite_float_or_default(record.get(name), float("nan"))
                if not np.isfinite(value) or not 0.0 <= value <= 1.0:
                    raise ValueError(f"AFDB manifest record has invalid {name}; expected a value in [0, 1]")

    def _select_chains(self, loader: PDBLoader, record: Dict[str, object], chain_ids: list[str]):
        mode = self.config.chain_selection_mode
        if mode == "auto_contact":
            chain_a, chain_b, details = loader.select_contact_chain_pair(
                min_chain_residues=self.config.min_chain_residues,
                distance_cutoff=self.config.topology.distance_cutoff,
            )
        elif mode == "manifest":
            chain_a = str(record.get("chain_a") or "").strip()
            chain_b = str(record.get("chain_b") or "").strip()
            details = {"manifest_record": record}
        else:
            chain_a = str(self.config.chain_a).strip()
            chain_b = str(self.config.chain_b).strip()
            details = {}
        if chain_a == chain_b or chain_a not in chain_ids or chain_b not in chain_ids:
            raise ValueError(f"Selected chains {chain_a}/{chain_b} are invalid; available={chain_ids}")
        residue_a = loader.get_chain_residue_count(chain_a)
        residue_b = loader.get_chain_residue_count(chain_b)
        if min(residue_a, residue_b) < self.config.min_chain_residues:
            raise ValueError(
                f"Selected chains too short: {chain_a}={residue_a}, {chain_b}={residue_b}; "
                f"minimum={self.config.min_chain_residues}"
            )
        return chain_a, chain_b, details

    def _confidence_preflight(
        self,
        record: Dict[str, object],
        surface_atoms,
        partner_atoms,
    ) -> Dict[str, object]:
        structure_type = str(record.get("structure_type") or "experimental").strip().lower()
        if not self.config.formal_mode or structure_type not in PREDICTED_STRUCTURE_TYPES:
            return {"status": "not_declared"}
        metric = str(record.get("confidence_metric") or "").strip().lower()
        per_chain = {}
        residue_arrays = []
        atom_arrays = []
        for label, atoms in (("surface_chain", surface_atoms), ("partner_chain", partner_atoms)):
            atom_values = self._atom_confidence_values(range(len(atoms)), atoms, metric)
            if (
                not len(atom_values)
                or len(atom_values) != len(atoms)
                or np.any((atom_values < 0.0) | (atom_values > 100.0))
            ):
                raise ValueError(
                    "Predicted structure does not provide finite 0-100 pLDDT values "
                    f"for every retained {label.replace('_', '-')} atom B-factor"
                )
            residue_values = np.asarray(residue_plddt_values(atoms), dtype=np.float64)
            per_chain[label] = {
                "summary_unit": "residue",
                "atom_count": int(len(atom_values)),
                "residue_count": int(len(residue_values)),
                "minimum": float(np.min(residue_values)),
                "mean": float(np.mean(residue_values)),
                "maximum": float(np.max(residue_values)),
                "atom_weighted_mean": float(np.mean(atom_values)),
            }
            residue_arrays.append(residue_values)
            atom_arrays.append(atom_values)
        residue_values = np.concatenate(residue_arrays)
        atom_values = np.concatenate(atom_arrays)
        return {
            "status": "validated_plddt_bfactor",
            "summary_unit": "residue",
            "atom_count": int(len(atom_values)),
            "residue_count": int(len(residue_values)),
            "minimum": float(np.min(residue_values)),
            "mean": float(np.mean(residue_values)),
            "maximum": float(np.max(residue_values)),
            "atom_weighted_mean": float(np.mean(atom_values)),
            "chains": per_chain,
        }

    def _interaction_job_metadata(
        self,
        record: Dict[str, object],
        *,
        structure_path: str,
        chain_a: str,
        chain_b: str,
        input_sha256: str,
        atoms_a,
        atoms_b,
    ) -> Dict[str, object]:
        interaction_value = record.get("prolif_file") or ""
        if not interaction_value:
            requires_interactions = self.config.execution_profile.strip().lower() != "operational_optcuts" or bool(
                set(self.config.resolved_optcuts_variants()) & set(RESIDUE_AWARE_OPTCUTS_METHODS)
            )
            if self.config.formal_mode and requires_interactions:
                raise ValueError(
                    "Formal benchmarks with residue-aware objectives or metrics require prolif_file and prolif_sha256"
                )
            return {}
        interaction_path = Path(str(interaction_value))
        if not interaction_path.is_absolute():
            base = (
                Path(self.config.manifest_path).resolve().parent
                if self.config.manifest_path
                else Path(structure_path).parent
            )
            interaction_path = base / interaction_path
        if not interaction_path.is_file():
            raise ValueError(f"Declared interaction file does not exist: {interaction_path}")
        expected_sha256 = str(record.get("prolif_sha256") or "").strip().lower()
        if self.config.formal_mode and not expected_sha256:
            raise ValueError("Formal manifest interaction file is missing prolif_sha256")
        actual_sha256 = sha256_file(interaction_path)
        if expected_sha256 and actual_sha256.lower() != expected_sha256:
            raise ValueError(f"Interaction checksum mismatch: expected {expected_sha256}, got {actual_sha256}")
        partner_map = load_prolif_partner_map(
            str(interaction_path),
            atoms_a,
            atoms_b,
            expected_chain_a=chain_a,
            expected_chain_b=chain_b,
            expected_source_sha256=input_sha256,
            require_bindings=bool(self.config.formal_mode),
        )
        if not partner_map:
            raise ValueError("Interaction file contains no residue pair that resolves to both selected chains")
        return {
            "prolif_file": str(interaction_path.resolve()),
            "prolif_sha256": actual_sha256,
            "interaction_partner_map": partner_map,
        }

    def _manifest_entry_skip(self, filename: str, record: Dict[str, object]) -> Dict[str, object] | None:
        if not record:
            return {
                "pdb": filename,
                "available_chains": [],
                "reason": "Structure is absent from the explicit manifest",
                "fatal_integrity_error": bool(self.config.formal_mode),
            }
        if self._manifest_record_is_included(record):
            return None
        return {
            "pdb": filename,
            "available_chains": [],
            "reason": str(record.get("exclusion_reason") or "Excluded by dataset manifest"),
        }

    def _missing_manifest_entries(
        self,
        structure_files: List[str],
        manifest: Dict[str, Dict[str, object]],
    ) -> List[Dict[str, object]]:
        present_names = set(structure_files) | {Path(name).stem for name in structure_files}
        declared_records = {id(record): record for record in manifest.values()}.values()
        missing = []
        for record in declared_records:
            declared_name = str(record.get("_manifest_filename") or "").strip()
            if not declared_name or declared_name in present_names or Path(declared_name).stem in present_names:
                continue
            included = self._manifest_record_is_included(record)
            missing.append(
                {
                    "pdb": declared_name,
                    "available_chains": [],
                    "reason": (
                        "Included structure declared by manifest is missing from input folder"
                        if included
                        else "Excluded manifest record has no local structure (expected)"
                    ),
                    "fatal_integrity_error": bool(included),
                }
            )
        return missing

    def _prepare_benchmark_jobs(self, structure_files: List[str]) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
        manifest = self._load_manifest_records() if self.config.chain_selection_mode == "manifest" else {}
        if self.config.formal_mode and manifest:
            self._validate_manifest_cohort(manifest)
        accepted_jobs: List[Dict[str, object]] = []
        accepted_record_ids: set[str] = set()
        skipped_files: List[Dict[str, object]] = []
        for filename in structure_files:
            self._check_cancelled()
            path = os.path.join(self.config.input_folder, filename)
            record = manifest.get(filename) or manifest.get(Path(filename).stem) or {}
            mode = self.config.chain_selection_mode
            if mode == "manifest" and (skip := self._manifest_entry_skip(filename, record)):
                skipped_files.append(skip)
                continue
            try:
                loader = PDBLoader(path)
                chain_ids = loader.get_protein_chain_ids()
            except Exception as exc:
                skipped_files.append(
                    {
                        "pdb": filename,
                        "reason": f"Failed to parse structure: {exc}",
                        "fatal_integrity_error": bool(self.config.formal_mode and mode == "manifest"),
                    }
                )
                continue
            if len(chain_ids) < 2:
                skipped_files.append(
                    {
                        "pdb": filename,
                        "reason": f"Need >=2 protein chains, found {len(chain_ids)}",
                        "fatal_integrity_error": bool(self.config.formal_mode and mode == "manifest"),
                    }
                )
                continue

            actual_sha256 = sha256_file(path)
            try:
                if mode == "manifest":
                    self._validate_manifest_record(record, actual_sha256=actual_sha256)
                chain_a, chain_b, selection_details = self._select_chains(loader, record, chain_ids)
                selected_coords, selected_atoms = loader.get_chain_data(chain_a)
                _partner_coords, partner_atoms = loader.get_chain_data(chain_b)
                confidence_preflight = self._confidence_preflight(
                    record,
                    selected_atoms,
                    partner_atoms,
                )
            except (OSError, RuntimeError, ValueError) as exc:
                skipped_files.append(
                    {
                        "pdb": filename,
                        "available_chains": chain_ids,
                        "reason": str(exc),
                        "fatal_integrity_error": bool(self.config.formal_mode and mode == "manifest"),
                    }
                )
                continue
            surface_estimate = SurfaceGenerator.estimate_grid(selected_coords, self.config.surface)
            if surface_estimate.get("status") != "ok":
                skipped_files.append(
                    {
                        "pdb": filename,
                        "available_chains": chain_ids,
                        "reason": (f"Surface grid preflight failed: {surface_estimate.get('status', 'unknown')}"),
                        "surface_grid_estimate": surface_estimate,
                        "fatal_integrity_error": bool(self.config.formal_mode),
                    }
                )
                continue

            job = {
                "pdb": filename,
                "input_sha256": actual_sha256,
                "chain_a": chain_a,
                "chain_b": chain_b,
                "selection_mode": mode,
                "selection_details": selection_details,
                "available_chains": chain_ids,
                "cluster_id": str(record.get("cluster_id") or record.get("interface_cluster_id") or "").strip(),
                "family_id": str(record.get("family_id") or "").strip(),
                "sequence_cluster_a": str(record.get("sequence_cluster_a") or "").strip(),
                "sequence_cluster_b": str(record.get("sequence_cluster_b") or "").strip(),
                "inference_sequence_cluster_a": str(record.get("inference_sequence_cluster_a") or "").strip(),
                "inference_sequence_cluster_b": str(record.get("inference_sequence_cluster_b") or "").strip(),
                "inference_family_id": str(record.get("inference_family_id") or "").strip(),
                "inference_dependency_basis": str(record.get("inference_dependency_basis") or "").strip(),
                "analysis_split": str(record.get("analysis_split") or "test").strip().lower(),
                "analysis_split_component_id": str(record.get("analysis_split_component_id") or "").strip(),
                "analysis_split_basis": str(record.get("analysis_split_basis") or "").strip(),
                "chain_a_residue_count": int(loader.get_chain_residue_count(chain_a)),
                "chain_b_residue_count": int(loader.get_chain_residue_count(chain_b)),
                "candidate_chain_pair_count": _first_present(record.get("candidate_chain_pair_count"), default=""),
                "selected_atom_contact_fraction": _first_present(
                    record.get("selected_atom_contact_fraction"), default=""
                ),
                "selected_residue_contact_fraction": _first_present(
                    record.get("selected_residue_contact_fraction"), default=""
                ),
                "hotspot_residues_a": record.get("hotspot_residues_a") or "",
                "dataset_source": record.get("dataset_source") or record.get("source") or "",
                "source_accession": record.get("source_accession") or record.get("accession") or "",
                "license_or_terms": record.get("license_or_terms") or record.get("license") or "",
                "structure_type": str(record.get("structure_type") or "experimental").strip().lower(),
                "structure_method": record.get("structure_method") or record.get("experimental_method") or "",
                "resolution_angstrom": record.get("resolution_angstrom") or record.get("resolution") or "",
                "experimental_methods_json": record.get("experimental_methods_json") or "",
                "experimental_method_group": record.get("experimental_method_group") or "",
                "experimental_method_contains_nmr": record.get("experimental_method_contains_nmr", ""),
                "pdbbind_index_resolution_angstrom": record.get("pdbbind_index_resolution_angstrom", ""),
                "rcsb_resolution_combined_angstrom_json": record.get("rcsb_resolution_combined_angstrom_json", ""),
                "rcsb_experiment_metadata_source": record.get("rcsb_experiment_metadata_source", ""),
                "confidence_metric": record.get("confidence_metric") or "",
                "confidence_stratum": record.get("confidence_stratum") or "",
                "afdb_ipsae_stratum": record.get("afdb_ipsae_stratum") or "",
                "confidence_source": record.get("confidence_source") or "",
                # Zero is a valid pLDDT threshold and must not be collapsed
                # to the empty/default value by truthiness.
                "confidence_threshold": record.get("confidence_threshold", ""),
                "afdb_model_id": _first_present(record.get("afdb_model_id"), record.get("model_id"), default=""),
                "afdb_iptm": _first_present(record.get("afdb_iptm"), record.get("iptm"), default=""),
                "afdb_ipsae": _first_present(record.get("afdb_ipsae"), record.get("ipsae"), default=""),
                "afdb_pdockq": _first_present(record.get("afdb_pdockq"), record.get("pdockq"), default=""),
                "afdb_pdockq2": _first_present(record.get("afdb_pdockq2"), record.get("pdockq2"), default=""),
                "afdb_lis": _first_present(record.get("afdb_lis"), record.get("lis"), default=""),
                "paired_record_id": record.get("paired_record_id") or "",
                "paired_experimental_record_id": record.get("paired_experimental_record_id") or "",
                "paired_geometry_stratum": record.get("paired_geometry_stratum") or "",
                "paired_contact_cutoff_angstrom": record.get("paired_contact_cutoff_angstrom", ""),
                "paired_predicted_contact_count_total": record.get("paired_predicted_contact_count_total", ""),
                "paired_contact_recall_fnat": record.get("paired_contact_recall_fnat", ""),
                "paired_contact_precision": record.get("paired_contact_precision", ""),
                "paired_contact_jaccard": record.get("paired_contact_jaccard", ""),
                "paired_experimental_contact_mapping_coverage": record.get(
                    "paired_experimental_contact_mapping_coverage", ""
                ),
                "paired_interface_residue_a_mapping_coverage": record.get(
                    "paired_interface_residue_a_mapping_coverage", ""
                ),
                "paired_interface_residue_b_mapping_coverage": record.get(
                    "paired_interface_residue_b_mapping_coverage", ""
                ),
                "paired_interface_ligand_ca_mapping_coverage": record.get(
                    "paired_interface_ligand_ca_mapping_coverage", ""
                ),
                "paired_interface_ligand_ca_rmsd_angstrom": record.get("paired_interface_ligand_ca_rmsd_angstrom", ""),
                "paired_cross_chain_clash_atom_fraction": record.get("paired_cross_chain_clash_atom_fraction", ""),
                "paired_alignment_a_optimal_correspondence_count": record.get(
                    "paired_alignment_a_optimal_correspondence_count", ""
                ),
                "paired_alignment_b_optimal_correspondence_count": record.get(
                    "paired_alignment_b_optimal_correspondence_count", ""
                ),
                "paired_alignment_a_selected_pair_consensus_fraction": record.get(
                    "paired_alignment_a_selected_pair_consensus_fraction", ""
                ),
                "paired_alignment_b_selected_pair_consensus_fraction": record.get(
                    "paired_alignment_b_selected_pair_consensus_fraction", ""
                ),
                "confidence_preflight": confidence_preflight,
                "surface_grid_estimate": surface_estimate,
                "surface_estimate_status": surface_estimate.get("status"),
                "surface_requested_voxel_count": surface_estimate.get("requested_voxel_count"),
                "surface_effective_voxel_count": surface_estimate.get("effective_voxel_count"),
                "surface_effective_resolution_angstrom": surface_estimate.get("effective_target_resolution_angstrom"),
                "surface_dense_field_bytes_lower_bound": surface_estimate.get(
                    "estimated_dense_field_bytes_lower_bound"
                ),
                "manifest_record_id": record.get("record_id") or record.get("id") or "",
            }

            try:
                job.update(
                    self._interaction_job_metadata(
                        record,
                        structure_path=path,
                        chain_a=chain_a,
                        chain_b=chain_b,
                        input_sha256=actual_sha256,
                        atoms_a=selected_atoms,
                        atoms_b=partner_atoms,
                    )
                )
            except (OSError, ValueError, json.JSONDecodeError) as exc:
                skipped_files.append(
                    {
                        "pdb": filename,
                        "available_chains": chain_ids,
                        "reason": f"Invalid declared interaction file: {exc}",
                        "fatal_integrity_error": True,
                    }
                )
                continue
            manifest_record_id = str(job.get("manifest_record_id") or "").strip()
            if self.config.formal_mode and manifest_record_id in accepted_record_ids:
                skipped_files.append(
                    {
                        "pdb": filename,
                        "available_chains": chain_ids,
                        "reason": (
                            "Formal manifest record_id is matched by more than one local structure: "
                            f"{manifest_record_id}"
                        ),
                        "fatal_integrity_error": True,
                    }
                )
                continue
            accepted_jobs.append(job)
            if manifest_record_id:
                accepted_record_ids.add(manifest_record_id)

        if self.config.chain_selection_mode == "manifest":
            skipped_files.extend(self._missing_manifest_entries(structure_files, manifest))

        return accepted_jobs, {
            "total_files": int(len(structure_files)),
            "accepted_files": int(len(accepted_jobs)),
            "skipped_files": int(len(skipped_files)),
            "integrity_error_count": int(sum(bool(item.get("fatal_integrity_error")) for item in skipped_files)),
            "accepted": accepted_jobs,
            "skipped": skipped_files,
            "rules": [
                "Input suffix must be .pdb, .cif, or .mmcif.",
                f"Each selected chain must contain >= {self.config.min_chain_residues} recognized amino-acid residues.",
                f"Chain selection mode: {self.config.chain_selection_mode}.",
                "Every accepted input and chain pair is written to benchmark_manifest.csv.",
                "Formal mode requires per-input SHA-256 and a cluster identifier.",
                "Formal included records require source, accession, license/terms, and structure type metadata.",
                "Unlisted or invalid inputs are fatal integrity errors in formal manifest mode.",
                "Surface voxel/memory estimates are computed without allocating a dense grid.",
            ],
        }

    @staticmethod
    def _manifest_record_is_included(record: Dict[str, object]) -> bool:
        value = str(record.get("include", record.get("status", "included"))).strip().lower()
        return value not in {"0", "false", "no", "exclude", "excluded", "skip", "skipped"}

    def _validate_manifest_cohort(self, manifest: Dict[str, Dict[str, object]]) -> None:
        """Validate inferential grouping and split invariants across included rows."""

        records = [
            record
            for record in {id(record): record for record in manifest.values()}.values()
            if self._manifest_record_is_included(record)
        ]
        record_ids: set[str] = set()
        pair_to_family: dict[tuple[str, str], str] = {}
        family_to_pair: dict[str, tuple[str, str]] = {}
        inference_pair_to_family: dict[tuple[str, str], str] = {}
        inference_family_to_pair: dict[str, tuple[str, str]] = {}
        split_by_sequence_cluster: dict[str, str] = {}
        split_by_component: dict[str, str] = {}
        split_by_family: dict[str, str] = {}
        split_by_inference_sequence_cluster: dict[str, str] = {}
        split_by_inference_family: dict[str, str] = {}
        split_by_analysis_component: dict[str, str] = {}
        component_by_sequence_cluster: dict[str, str] = {}
        component_by_family: dict[str, str] = {}
        component_by_homology_group: dict[str, str] = {}
        component_by_inference_sequence_cluster: dict[str, str] = {}
        component_by_inference_family: dict[str, str] = {}
        cohort_split_basis: dict[str, str] = {}

        def require_one(mapping: dict, key, value, label: str) -> None:
            previous = mapping.setdefault(key, value)
            if previous != value:
                raise ValueError(f"Formal manifest assigns {label} {key!r} to both {previous!r} and {value!r}")

        for record in records:
            split = str(record.get("analysis_split") or "").strip().lower()
            record_id = str(record.get("record_id") or record.get("id") or "").strip()
            if not record_id:
                raise ValueError("Formal manifest contains an included record without record_id/id")
            if record_id in record_ids:
                raise ValueError(f"Duplicate formal manifest record_id: {record_id}")
            record_ids.add(record_id)

            analysis_component = str(record.get("analysis_split_component_id") or "").strip()
            analysis_basis = str(record.get("analysis_split_basis") or "").strip()
            if self.config.formal_mode and (not analysis_component or not analysis_basis):
                raise ValueError("Formal manifest has incomplete split-dependence metadata")
            if analysis_component:
                require_one(
                    split_by_analysis_component,
                    analysis_component,
                    split,
                    "analysis-split dependency component",
                )
            if analysis_basis:
                require_one(cohort_split_basis, "cohort", analysis_basis, "analysis_split_basis")

            inference_values = {name: str(record.get(name) or "").strip() for name in INFERENCE_DEPENDENCY_FIELDS}
            if any(inference_values.values()):
                missing = [name for name, value in inference_values.items() if not value]
                if missing:
                    raise ValueError(
                        "Formal manifest has incomplete prediction-dependence metadata: " + ", ".join(missing)
                    )
                inference_pair = tuple(
                    sorted(
                        (
                            inference_values["inference_sequence_cluster_a"],
                            inference_values["inference_sequence_cluster_b"],
                        )
                    )
                )
                inference_family = inference_values["inference_family_id"]
                if inference_family != inference_family_id(*inference_pair):
                    raise ValueError(f"Formal manifest has invalid inference_family_id {inference_family!r}")
                if inference_values["inference_dependency_basis"] != INFERENCE_DEPENDENCY_BASIS:
                    raise ValueError("Formal manifest has an invalid inference_dependency_basis")
                require_one(
                    inference_pair_to_family,
                    inference_pair,
                    inference_family,
                    "unordered prediction-dependency partner pair",
                )
                require_one(
                    inference_family_to_pair,
                    inference_family,
                    inference_pair,
                    "inference_family_id",
                )
                require_one(
                    split_by_inference_sequence_cluster,
                    inference_values["inference_sequence_cluster_a"],
                    split,
                    "prediction-dependency sequence cluster",
                )
                require_one(
                    split_by_inference_sequence_cluster,
                    inference_values["inference_sequence_cluster_b"],
                    split,
                    "prediction-dependency sequence cluster",
                )
                require_one(
                    split_by_inference_family,
                    inference_family,
                    split,
                    "prediction-dependency family",
                )
                require_one(
                    component_by_inference_sequence_cluster,
                    inference_values["inference_sequence_cluster_a"],
                    analysis_component,
                    "prediction-dependency sequence-cluster component",
                )
                require_one(
                    component_by_inference_sequence_cluster,
                    inference_values["inference_sequence_cluster_b"],
                    analysis_component,
                    "prediction-dependency sequence-cluster component",
                )
                require_one(
                    component_by_inference_family,
                    inference_family,
                    analysis_component,
                    "prediction-dependency family component",
                )

            sequence_a = str(record.get("sequence_cluster_a") or "").strip()
            sequence_b = str(record.get("sequence_cluster_b") or "").strip()
            family = str(record.get("family_id") or "").strip()
            component = str(record.get("cluster_id") or record.get("interface_cluster_id") or "").strip()
            if not sequence_a or not sequence_b or not family or not component:
                continue

            pair = tuple(sorted((sequence_a, sequence_b)))
            require_one(pair_to_family, pair, family, "unordered partner pair")
            require_one(family_to_pair, family, pair, "family_id")
            require_one(split_by_sequence_cluster, sequence_a, split, "sequence cluster")
            require_one(split_by_sequence_cluster, sequence_b, split, "sequence cluster")
            require_one(split_by_component, component, split, "homology component")
            require_one(split_by_family, family, split, "interaction family")
            if analysis_component:
                require_one(
                    component_by_sequence_cluster,
                    sequence_a,
                    analysis_component,
                    "sequence-cluster component",
                )
                require_one(
                    component_by_sequence_cluster,
                    sequence_b,
                    analysis_component,
                    "sequence-cluster component",
                )
                require_one(component_by_family, family, analysis_component, "interaction-family component")
                require_one(component_by_homology_group, component, analysis_component, "homology-group component")

    def _load_manifest_records(self) -> Dict[str, Dict[str, object]]:
        path = Path(self.config.manifest_path)
        if path.suffix.lower() == ".json":
            with open(path, "r", encoding="utf-8") as handle:
                payload = json.load(handle)
            records = payload.get("records", payload.get("files", [])) if isinstance(payload, dict) else payload
        else:
            with open(path, "r", encoding="utf-8-sig", newline="") as handle:
                records = list(csv.DictReader(handle))
        if not isinstance(records, list):
            raise ValueError("Benchmark manifest records/files must be a list.")
        result = {}
        for raw in records:
            if not isinstance(raw, dict):
                raise ValueError("Every benchmark manifest record must be an object/CSV row.")
            record = dict(raw)
            name = str(record.get("pdb") or record.get("file") or record.get("filename") or "").strip()
            if not name:
                if self.config.formal_mode:
                    raise ValueError("Formal benchmark manifest contains a record without pdb/file/filename.")
                continue
            record["_manifest_filename"] = name
            if name in result:
                raise ValueError(f"Duplicate benchmark manifest record: {name}")
            result[name] = record
            stem = Path(name).stem
            if stem in result and result[stem] is not record:
                raise ValueError(f"Ambiguous benchmark manifest stem: {stem}")
            result[stem] = record
        return result

    def _load_resume_state(
        self, jobs: List[Dict[str, object]]
    ) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
        if not bool(self.config.resume) or not os.path.exists(self._checkpoint_path):
            return [], jobs
        try:
            with open(self._checkpoint_path, "r", encoding="utf-8") as handle:
                checkpoint = json.load(handle)
        except (OSError, json.JSONDecodeError) as exc:
            if self.config.formal_mode:
                raise ValueError(f"Could not read the matching formal resume checkpoint: {exc}") from exc
            self.log(f"[Benchmark] Failed to load checkpoint: {exc}")
            return [], jobs
        if not isinstance(checkpoint, dict) or not isinstance(checkpoint.get("files"), list):
            if self.config.formal_mode:
                raise ValueError("Formal resume checkpoint has an invalid root or files list.")
            self.log("[Benchmark] Checkpoint has an invalid root or files list; ignoring it.")
            return [], jobs
        if str(checkpoint.get("config_fingerprint") or "") != self._checkpoint_fingerprint:
            self.log("[Benchmark] Checkpoint fingerprint differs; ignoring old results.")
            return [], jobs

        jobs_by_name = {str(job.get("pdb")): job for job in jobs}
        checkpoint_names = []
        duplicate_names = set()
        seen_names = set()
        for item in checkpoint["files"]:
            if not isinstance(item, dict) or not item.get("pdb"):
                if self.config.formal_mode:
                    raise ValueError("Formal resume checkpoint contains a record without a structure name.")
                continue
            name = str(item["pdb"])
            checkpoint_names.append(name)
            if name in seen_names:
                duplicate_names.add(name)
            seen_names.add(name)
        if duplicate_names:
            raise ValueError(
                "Resume checkpoint contains duplicate structure records: " + ", ".join(sorted(duplicate_names))
            )
        unknown_names = sorted(set(checkpoint_names) - set(jobs_by_name))
        if self.config.formal_mode and unknown_names:
            raise ValueError(
                "Resume checkpoint contains structures outside the current formal cohort: " + ", ".join(unknown_names)
            )
        completed_records = []
        mismatched_names = []
        for item in checkpoint["files"]:
            if self._resume_record_matches(item, jobs_by_name):
                completed_records.append(item)
            elif isinstance(item, dict) and str(item.get("pdb") or "") in jobs_by_name:
                mismatched_names.append(str(item["pdb"]))
        if self.config.formal_mode and mismatched_names:
            raise ValueError(
                "Resume checkpoint records do not match current input identity: " + ", ".join(sorted(mismatched_names))
            )
        completed_names = {str(item["pdb"]) for item in completed_records}
        remaining = [job for job in jobs if str(job["pdb"]) not in completed_names]
        return completed_records, remaining

    @staticmethod
    def _resume_record_matches(item, jobs_by_name) -> bool:
        if not isinstance(item, dict) or not item.get("pdb"):
            return False
        if str(item.get("status") or "") not in {"ok", "failed", "incomplete_comparison"}:
            return False
        name = str(item["pdb"])
        job = jobs_by_name.get(name)
        selection = item.get("chain_selection")
        if job is None or not isinstance(selection, dict):
            return False
        identity_matches = all(
            str(item.get(field) or "") == str(job.get(field) or "")
            for field in (*RESULT_IDENTITY_FIELDS, "structure_type")
        )
        return bool(
            str(item.get("input_sha256") or "").lower() == str(job.get("input_sha256") or "").lower()
            and str(item.get("interaction_sha256") or "").lower() == str(job.get("prolif_sha256") or "").lower()
            and str(selection.get("chain_a") or "") == str(job.get("chain_a") or "")
            and str(selection.get("chain_b") or "") == str(job.get("chain_b") or "")
            and identity_matches
        )

    def _save_checkpoint(self, results: List[Dict[str, object]]) -> bool:
        payload = {
            "created_at": _utc_now(),
            "topoppi_version": __version__,
            "config_fingerprint": self._checkpoint_fingerprint,
            "files": results,
        }
        try:
            dump_json_atomic(payload, self._checkpoint_path, indent=None)
            return True
        except OSError as exc:
            self.log(f"[Benchmark] Failed to save checkpoint: {exc}")
            return False

    def _memory_rss_mb(self) -> float:
        if self._proc is None:
            return 0.0
        try:
            return float(self._proc.memory_info().rss) / (1024.0 * 1024.0)
        except psutil.Error:
            return 0.0

    def _environment_metadata(self) -> Dict[str, object]:
        package_names = (
            "topoppi",
            "numpy",
            "scipy",
            "trimesh",
            "scikit-image",
            "libigl",
            "biopython",
            "shapely",
            "psutil",
        )
        versions = {}
        for name in package_names:
            try:
                versions[name] = importlib_metadata.version(name)
            except importlib_metadata.PackageNotFoundError:
                versions[name] = None

        repo_root = Path(__file__).resolve().parents[3]
        git_commit, git_dirty = git_worktree_state(repo_root)

        return {
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "logical_cpu_count": os.cpu_count(),
            "total_ram_bytes": int(psutil.virtual_memory().total) if psutil is not None else None,
            "package_versions": versions,
            "git_commit": git_commit,
            "git_worktree_dirty": git_dirty,
            "invocation": [sys.executable, *sys.argv],
            "thread_environment": {
                name: os.environ.get(name)
                for name in (
                    "OMP_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS",
                    "NUMEXPR_NUM_THREADS",
                )
            },
            "cpu_affinity": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
        }

    @staticmethod
    def _stage_stats(start_wall: float, start_cpu: float) -> Dict[str, float]:
        return {
            "wall_sec": float(time.perf_counter() - start_wall),
            "cpu_sec": float(_worker_cpu_time() - start_cpu),
        }

    @staticmethod
    def _from_timing_list(wall: List[float], cpu: List[float]) -> Dict[str, float]:
        if not wall:
            return {"wall_sec": float("inf"), "cpu_sec": float("inf")}
        return {"wall_sec": float(np.sum(wall)), "cpu_sec": float(np.sum(cpu))}

    @staticmethod
    def _sum_stage_stats(*stages: Dict[str, float]) -> Dict[str, float]:
        return {
            "wall_sec": float(sum(float(stage.get("wall_sec", 0.0)) for stage in stages)),
            "cpu_sec": float(sum(float(stage.get("cpu_sec", 0.0)) for stage in stages)),
        }

    def _parameterization_timing_block(self, stages: Dict[str, Dict[str, float]]) -> Dict[str, float]:
        result = {}
        for method in self._configured_methods():
            if method not in PARAMETERIZATION_METHODS:
                continue
            block = stages.get(f"{method}_parameterization", {})
            result[f"{method}_total_wall_sec"] = float(block.get("wall_sec", float("inf")))
            result[f"{method}_total_cpu_sec"] = float(block.get("cpu_sec", float("inf")))
        return result

    def _log_thread(self, message: str) -> None:
        thread = threading.current_thread()
        self.log(f"[Benchmark][Supervisor {thread.name}:{thread.ident}] {message}")

    def _safe_progress(self, completed: int, total: int, message: str) -> None:
        try:
            self.progress(int(completed), int(total), str(message))
        except Exception:
            pass
