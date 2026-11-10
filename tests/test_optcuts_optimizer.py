import os
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import numpy as np
import trimesh
from PIL import Image

from topoppi.atlas.metrics import UVAtlasMetrics
from topoppi.atlas.uv import set_uv_layout
from topoppi.config import OptCutsConfig, ParameterizationConfig
from topoppi.mesh.parameterization import Parameterizer
from topoppi.mesh.provenance import OPTCUTS_GEOMETRY_VERTEX_IDS, initialize_provenance
from topoppi.optimization.optcuts.joint_optimizer import OptCutsUVOptimizer, resolve_optcuts_binary

OPTCUTS_BIN = Path(__file__).parents[1] / "tools" / "OptCuts" / "OptCuts_bin"
RESIDUE_AWARE_OPTCUTS_BIN = Path(__file__).parents[1] / "tools" / "OptCuts" / "OptCuts_bin"


class OptCutsOptimizerTests(unittest.TestCase):
    @staticmethod
    def _fake_optcuts_binary(folder: str, *, advertise_footprints: bool) -> Path:
        binary = Path(folder) / "fake_optcuts"
        marker = "print('residue footprint energy enabled')" if advertise_footprints else ""
        binary.write_text(
            "#!/usr/bin/env python3\n"
            "import os, shutil, sys\n"
            "assert os.path.basename(sys.argv[2]) == sys.argv[2], sys.argv[2]\n"
            "os.makedirs('output', exist_ok=True)\n"
            "shutil.copy2(sys.argv[2], 'output/finalResult_mesh.obj')\n"
            f"{marker}\n",
            encoding="utf-8",
        )
        binary.chmod(0o755)
        return binary

    def test_binary_preflight_hashes_once_per_optimizer(self):
        with tempfile.TemporaryDirectory() as tmp:
            binary = Path(tmp) / "OptCuts_bin"
            binary.write_bytes(b"binary")
            binary.chmod(0o755)
            optimizer = OptCutsUVOptimizer(replace(OptCutsConfig(), optcuts_bin=str(binary)))
            with patch(
                "topoppi.optimization.optcuts.joint_optimizer.sha256_file",
                return_value="a" * 64,
            ) as digest:
                first = optimizer.preflight_binary()
                second = optimizer.preflight_binary()

        self.assertEqual(first, second)
        digest.assert_called_once_with(str(binary.resolve()))

    def test_residue_capability_is_checked_once_per_optimizer(self):
        with tempfile.TemporaryDirectory() as tmp:
            binary = Path(tmp) / "OptCuts_bin"
            binary.write_bytes(b"residue footprint energy enabled")
            binary.chmod(0o755)
            optimizer = OptCutsUVOptimizer(
                replace(
                    OptCutsConfig(),
                    optcuts_bin=str(binary),
                    residue_fragmentation_weight=1.0,
                )
            )
            with patch(
                "topoppi.optimization.optcuts.joint_optimizer.supports_residue_footprint_energy",
                return_value=True,
            ) as capability:
                optimizer.preflight_binary()
                optimizer.preflight_binary()

        capability.assert_called_once_with(str(binary.resolve()))

    def test_missing_binary_error_explains_how_to_install_or_select_optcuts(self):
        optimizer = OptCutsUVOptimizer(replace(OptCutsConfig(), optcuts_bin="missing-optcuts"))
        with patch.dict(
            "topoppi.optimization.optcuts.joint_optimizer.os.environ",
            {"TOPOPPI_OPTCUTS_BIN": "missing-native-optcuts"},
        ):
            with self.assertRaises(RuntimeError) as raised:
                optimizer.preflight_binary()

        message = str(raised.exception)
        self.assertIn("topoppi-install-optcuts", message)
        self.assertIn("TOPOPPI_OPTCUTS_BIN", message)
        self.assertIn("native OptCuts executable", message)

    def test_default_binary_is_found_in_source_checkout_outside_working_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            checkout = Path(tmp) / "checkout"
            filename = "OptCuts_bin.exe" if os.name == "nt" else "OptCuts_bin"
            binary = checkout / "tools" / "OptCuts" / filename
            binary.parent.mkdir(parents=True)
            binary.write_bytes(b"binary")
            module = checkout / "src" / "topoppi" / "optimization" / "optcuts" / "joint_optimizer.py"
            with (
                patch.dict("topoppi.optimization.optcuts.joint_optimizer.os.environ", {}, clear=True),
                patch("topoppi.optimization.optcuts.joint_optimizer.shutil.which", return_value=None),
                patch("topoppi.optimization.optcuts.joint_optimizer.__file__", str(module)),
                patch("topoppi.optimization.optcuts.joint_optimizer.sys.executable", str(Path(tmp) / "python")),
                patch("topoppi.optimization.optcuts.joint_optimizer.default_install_dir", return_value=Path(tmp) / "bin"),
            ):
                self.assertEqual(resolve_optcuts_binary(OptCutsConfig()), str(binary))
                missing = replace(OptCutsConfig(), optcuts_bin=str(Path(tmp) / "missing"))
                self.assertIsNone(resolve_optcuts_binary(missing))

    def test_explicit_binary_overrides_environment_and_reports_its_own_missing_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            environmental, explicit = Path(tmp) / "environmental", Path(tmp) / "explicit"
            environmental.write_bytes(b"binary")
            explicit.write_bytes(b"binary")
            with patch.dict(os.environ, {"TOPOPPI_OPTCUTS_BIN": str(environmental)}):
                self.assertEqual(resolve_optcuts_binary(OptCutsConfig()), str(environmental))
                config = replace(OptCutsConfig(), optcuts_bin=str(explicit))
                self.assertEqual(resolve_optcuts_binary(config), str(explicit))
                explicit.unlink()
                self.assertIsNone(resolve_optcuts_binary(config))
                with self.assertRaisesRegex(RuntimeError, str(explicit)):
                    OptCutsUVOptimizer(config).preflight_binary()

    def test_automatic_discovery_prefers_residue_capability_and_respects_explicit_choice(self):
        with tempfile.TemporaryDirectory() as tmp:
            filename = "OptCuts_bin.exe" if os.name == "nt" else "OptCuts_bin"
            plain = Path(tmp) / "bin" / filename
            plain.parent.mkdir()
            plain.write_bytes(b"ordinary OptCuts")
            checkout = Path(tmp) / "checkout"
            patched = checkout / "tools" / "OptCuts" / filename
            patched.parent.mkdir(parents=True)
            patched.write_bytes(b"residue footprint energy enabled")
            module = checkout / "src" / "topoppi" / "optimization" / "optcuts" / "joint_optimizer.py"
            config = replace(OptCutsConfig(), residue_fragmentation_weight=20.0)
            with (
                patch.dict("topoppi.optimization.optcuts.joint_optimizer.os.environ", {}, clear=True),
                patch("topoppi.optimization.optcuts.joint_optimizer.shutil.which", return_value=str(plain)),
                patch("topoppi.optimization.optcuts.joint_optimizer.__file__", str(module)),
                patch("topoppi.optimization.optcuts.joint_optimizer.sys.executable", str(plain.parent / "python")),
                patch("topoppi.optimization.optcuts.joint_optimizer.default_install_dir", return_value=plain.parent),
            ):
                self.assertEqual(resolve_optcuts_binary(config), str(patched))
                self.assertEqual(resolve_optcuts_binary(replace(config, optcuts_bin=str(plain))), str(plain))
                with patch.dict("topoppi.optimization.optcuts.joint_optimizer.os.environ", {"TOPOPPI_OPTCUTS_BIN": str(plain)}):
                    self.assertEqual(resolve_optcuts_binary(config), str(plain))

    def test_default_binary_is_found_in_installer_directory_without_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            filename = "OptCuts_bin.exe" if os.name == "nt" else "OptCuts_bin"
            binary = Path(tmp) / filename
            binary.write_bytes(b"binary")
            with (
                patch.dict("topoppi.optimization.optcuts.joint_optimizer.os.environ", {}, clear=True),
                patch("topoppi.optimization.optcuts.joint_optimizer.shutil.which", return_value=None),
                patch("topoppi.optimization.optcuts.joint_optimizer.sys.executable", str(Path(tmp) / "other" / "python")),
                patch("topoppi.optimization.optcuts.joint_optimizer.default_install_dir", return_value=Path(tmp)),
            ):
                self.assertEqual(resolve_optcuts_binary(OptCutsConfig()), str(binary))
                with patch.dict("topoppi.optimization.optcuts.joint_optimizer.os.environ", {"TOPOPPI_OPTCUTS_BIN": "missing-choice"}):
                    self.assertIsNone(resolve_optcuts_binary(OptCutsConfig()))

    @unittest.skipIf(os.name == "nt", "POSIX executable permission")
    def test_preflight_reports_nonexecutable_binary_before_optimization(self):
        with tempfile.TemporaryDirectory() as tmp:
            binary = Path(tmp) / "OptCuts_bin"
            binary.write_bytes(b"binary")
            binary.chmod(0o644)
            optimizer = OptCutsUVOptimizer(replace(OptCutsConfig(), optcuts_bin=str(binary)))
            with self.assertRaisesRegex(RuntimeError, "not executable"):
                optimizer.preflight_binary()

    def test_frame_resize_targets_the_long_edge(self):
        image = Image.new("RGB", (100, 50))

        resized = OptCutsUVOptimizer._resize_image_if_needed(image, min_long_edge=200)

        self.assertEqual(resized.size, (200, 100))

    def test_benchmark_mode_can_skip_duplicate_quality_report(self):
        mesh = trimesh.Trimesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2]]),
            process=False,
        )
        corners = np.asarray(mesh.vertices[:, :2])[mesh.faces]
        set_uv_layout(mesh, corners)
        optimizer = OptCutsUVOptimizer()
        with (
            patch.object(
                optimizer,
                "_run_optcuts_for_patch",
                return_value=(corners, {"status": "ok"}),
            ),
            patch.object(optimizer, "_build_report") as build_report,
        ):
            optimizer.optimize_patches(
                [mesh],
                initialization="provided",
                pack=False,
                build_report=False,
            )

        build_report.assert_not_called()
        self.assertEqual(optimizer.get_last_report(), {})
        self.assertIn("uv_optcuts", mesh.metadata)
        self.assertNotIn("uv_global", mesh.metadata)

    def test_per_call_timeout_override_reaches_the_solver_invocation(self):
        mesh = trimesh.Trimesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2]]),
            process=False,
        )
        corners = np.asarray(mesh.vertices[:, :2])[mesh.faces]
        optimizer = OptCutsUVOptimizer()
        with patch.object(
            optimizer,
            "_run_optcuts_for_patch",
            return_value=(corners, {"status": "ok"}),
        ) as run_patch:
            optimizer.optimize_patches(
                [mesh],
                initialization="automatic",
                pack=False,
                build_report=False,
                timeout_sec=12.5,
            )

        self.assertEqual(run_patch.call_args.kwargs["timeout_sec"], 12.5)

    def test_bijective_run_rejects_a_globally_overlapping_initial_map(self):
        mesh = trimesh.Trimesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [2.0, 0.0, 0.0],
                    [3.0, 0.0, 0.0],
                    [2.0, 1.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [3, 4, 5]]),
            process=False,
        )
        repeated = np.repeat(
            np.array([[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]]),
            2,
            axis=0,
        )
        set_uv_layout(mesh, repeated)
        optimizer = OptCutsUVOptimizer()

        with self.assertRaisesRegex(RuntimeError, "not globally injective"):
            optimizer.optimize_patches(
                [mesh],
                initialization="provided",
                pack=False,
                build_report=False,
            )

    def test_provided_global_reflection_is_made_positive_and_recorded(self):
        mesh = trimesh.Trimesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2]]),
            process=False,
        )
        reflected = np.asarray(mesh.vertices[:, :2]).copy()
        reflected[:, 0] *= -1.0
        set_uv_layout(mesh, reflected)
        with tempfile.TemporaryDirectory() as tmp:
            binary = self._fake_optcuts_binary(tmp, advertise_footprints=False)
            optimizer = OptCutsUVOptimizer(replace(OptCutsConfig(), optcuts_bin=str(binary)))
            optimizer.optimize_patches(
                [mesh],
                initialization="provided",
                pack=False,
                build_report=False,
            )

        execution = mesh.metadata["optcuts_execution"]
        self.assertEqual(
            execution["provided_uv_transform"],
            "global_u_reflection_for_optcuts_positive_orientation",
        )
        self.assertNotEqual(
            execution["source_initial_uv_checksum"],
            execution["initial_uv_checksum"],
        )
        self.assertFalse(execution["initial_uv_injectivity"]["global_reflection_required_for_positive_orientation"])

    def test_command_preserves_full_float_precision(self):
        mesh = trimesh.Trimesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2]]),
            process=False,
        )
        set_uv_layout(mesh, np.asarray(mesh.vertices[:, :2]))
        lambda_init = float(np.nextafter(0.9, 1.0))
        distortion_bound = float(np.nextafter(4.2, 5.0))
        with tempfile.TemporaryDirectory() as tmp:
            binary = self._fake_optcuts_binary(tmp, advertise_footprints=False)
            optimizer = OptCutsUVOptimizer(
                replace(
                    OptCutsConfig(),
                    optcuts_bin=str(binary),
                    optcuts_lambda_init=lambda_init,
                    optcuts_distortion_bound=distortion_bound,
                )
            )
            optimizer.optimize_patches(
                [mesh],
                initialization="provided",
                pack=False,
                build_report=False,
            )

        command = mesh.metadata["optcuts_execution"]["command"]
        self.assertEqual(float(command[3]), lambda_init)
        self.assertEqual(float(command[6]), distortion_bound)

    def test_output_must_independently_satisfy_the_optcuts_distortion_bound(self):
        mesh = trimesh.Trimesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2]]),
            process=False,
        )
        stretched = np.asarray([[0.0, 0.0], [10.0, 0.0], [0.0, 1.0]])
        set_uv_layout(mesh, stretched)
        with tempfile.TemporaryDirectory() as tmp:
            binary = self._fake_optcuts_binary(tmp, advertise_footprints=False)
            optimizer = OptCutsUVOptimizer(replace(OptCutsConfig(), optcuts_bin=str(binary)))

            with self.assertRaisesRegex(RuntimeError, "outside the requested distortion"):
                optimizer.optimize_patches(
                    [mesh],
                    initialization="provided",
                    pack=False,
                    build_report=False,
                )

    def test_uncollected_objective_trace_is_not_reported_as_zero(self):
        mesh = trimesh.Trimesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2]]),
            process=False,
        )
        corner_uv = np.asarray(mesh.vertices[:, :2])[mesh.faces]
        set_uv_layout(mesh, corner_uv, key="uv_optcuts")
        set_uv_layout(mesh, corner_uv, key="uv_global")

        report = OptCutsUVOptimizer()._build_report(
            [mesh],
            iteration_time=1.25,
            packing_report={"status": "disabled"},
        )

        stability = report["stability_efficiency"]
        self.assertIsNone(stability["objective_history"])
        self.assertIsNone(stability["objective_drop"])
        self.assertEqual(stability["objective_trace_status"], "not_collected_from_optcuts")

    def test_manual_obj_uv_parser(self):
        obj_text = """\
v 0 0 0
v 1 0 0
v 0 1 0
vt 0 0
vt 1 0
vt 0 1
f 1/1 2/2 3/3
"""
        with tempfile.TemporaryDirectory() as tmp:
            obj_path = Path(tmp) / "uv.obj"
            obj_path.write_text(obj_text, encoding="utf-8")

            parsed = OptCutsUVOptimizer._parse_obj_uv(str(obj_path))

        self.assertEqual(parsed.corner_uv.shape, (1, 3, 2))
        self.assertEqual(float(parsed.corner_uv[0, 1, 0]), 1.0)

    def test_manual_parser_preserves_multiple_uvs_for_one_vertex(self):
        obj_text = """\
v 0 0 0
v 1 0 0
v 1 1 0
v 0 1 0
vt 0 0
vt 1 0
vt 1 1
vt 4 5
vt 1 1
vt 0 1
f 1/1 2/2 3/3
f 1/4 3/5 4/6
"""
        with tempfile.TemporaryDirectory() as tmp:
            obj_path = Path(tmp) / "seam.obj"
            obj_path.write_text(obj_text, encoding="utf-8")
            parsed = OptCutsUVOptimizer._parse_obj_uv(str(obj_path))

        self.assertEqual(parsed.corner_uv.shape, (2, 3, 2))
        np.testing.assert_array_equal(parsed.corner_uv[0, 0], [0.0, 0.0])
        np.testing.assert_array_equal(parsed.corner_uv[1, 0], [4.0, 5.0])

    def test_obj_writer_emits_face_corner_texture_indices_only_when_initialized(self):
        mesh = trimesh.Trimesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2]]),
            process=False,
        )
        corners = np.array([[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]])
        with tempfile.TemporaryDirectory() as tmp:
            initialized = Path(tmp) / "initialized.obj"
            automatic = Path(tmp) / "automatic.obj"
            OptCutsUVOptimizer._write_obj_with_uv(mesh, str(initialized), corners)
            OptCutsUVOptimizer._write_obj_with_uv(mesh, str(automatic), None)
            initialized_text = initialized.read_text(encoding="utf-8")
            automatic_text = automatic.read_text(encoding="utf-8")

        self.assertEqual(initialized_text.count("\nvt "), 3)
        self.assertIn("f 1/1 2/2 3/3", initialized_text)
        self.assertNotIn("\nvt ", automatic_text)
        self.assertIn("f 1 2 3", automatic_text)

    def test_obj_writer_shares_texture_indices_across_continuous_edges(self):
        mesh = trimesh.Trimesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2], [0, 2, 3]]),
            process=False,
        )
        corners = np.asarray(mesh.vertices[:, :2])[mesh.faces]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "continuous.obj"
            OptCutsUVOptimizer._write_obj_with_uv(mesh, str(path), corners)
            text = path.read_text(encoding="utf-8")

        self.assertEqual(text.count("\nvt "), 4)
        self.assertIn("f 1/1 2/2 3/3", text)
        self.assertIn("f 1/1 3/3 4/4", text)

    def test_obj_writer_keeps_distinct_indices_across_a_real_seam(self):
        mesh = trimesh.Trimesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2], [0, 2, 3]]),
            process=False,
        )
        corners = np.asarray(mesh.vertices[:, :2])[mesh.faces]
        corners[1] += np.array([3.0, 0.0])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "seam.obj"
            OptCutsUVOptimizer._write_obj_with_uv(mesh, str(path), corners)
            text = path.read_text(encoding="utf-8")

        self.assertEqual(text.count("\nvt "), 6)
        self.assertIn("f 1/1 2/2 3/3", text)
        self.assertIn("f 1/5 3/6 4/4", text)

    def test_obj_writer_encodes_diskification_cuts_as_uv_seams_on_shared_geometry(self):
        mesh = trimesh.Trimesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [1.0, 1.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [4, 5, 3]]),
            process=False,
        )
        mesh.metadata["source_vertex_ids"] = np.array([10, 11, 12, 13, 10, 12])
        mesh.metadata[OPTCUTS_GEOMETRY_VERTEX_IDS] = np.array([0, 1, 2, 3, 0, 2])
        corners = np.asarray(mesh.vertices[:, :2])[mesh.faces]
        corners[1] += np.array([3.0, 0.0])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "topology_cut.obj"
            metadata = OptCutsUVOptimizer._write_obj_with_uv(mesh, str(path), corners)
            text = path.read_text(encoding="utf-8")

        self.assertEqual(text.count("\nv "), 4)
        self.assertEqual(text.count("\nvt "), 6)
        self.assertIn("f 1/1 2/2 3/3", text)
        self.assertIn("f 1/5 3/6 4/4", text)
        self.assertEqual(metadata["collapsed_vertex_copy_count"], 2)
        self.assertEqual(metadata["preserved_topology_vertex_copy_count"], 0)
        self.assertEqual(metadata["source_vertex_ids"], [10, 11, 12, 13])
        self.assertEqual(metadata["footprint_topology_vertex_ids"], [0, 1, 2, 3])

    def test_automatic_obj_writer_restores_pre_diskification_topology(self):
        mesh = trimesh.Trimesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [1.0, 1.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [4, 5, 3]]),
            process=False,
        )
        mesh.metadata["source_vertex_ids"] = np.array([10, 11, 12, 13, 10, 12])
        mesh.metadata[OPTCUTS_GEOMETRY_VERTEX_IDS] = np.array([0, 1, 2, 3, 0, 2])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "automatic.obj"
            OptCutsUVOptimizer._write_obj_with_uv(mesh, str(path), None)
            text = path.read_text(encoding="utf-8")

        self.assertEqual(text.count("\nv "), 4)
        self.assertIn("f 1 2 3", text)
        self.assertIn("f 1 3 4", text)

    def test_obj_writer_preserves_repaired_vertex_fan_copies(self):
        mesh = trimesh.Trimesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, -1.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [3, 4, 5]]),
            process=False,
        )
        mesh.metadata["source_vertex_ids"] = np.array([10, 11, 12, 10, 13, 14])
        mesh.metadata[OPTCUTS_GEOMETRY_VERTEX_IDS] = np.arange(6, dtype=np.int64)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "repaired_fans.obj"
            metadata = OptCutsUVOptimizer._write_obj_with_uv(mesh, str(path), None)
            text = path.read_text(encoding="utf-8")

        self.assertEqual(text.count("\nv "), 6)
        self.assertEqual(metadata["collapsed_vertex_copy_count"], 0)
        self.assertEqual(metadata["preserved_topology_vertex_copy_count"], 1)
        self.assertEqual(metadata["footprint_topology_vertex_ids"], list(range(6)))

    def test_obj_writer_keeps_base_texture_indices_aligned_with_vertices(self):
        mesh = trimesh.Trimesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[2, 0, 3], [0, 2, 1]]),
            process=False,
        )
        corners = np.asarray(mesh.vertices[:, :2])[mesh.faces]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "reordered_vertices.obj"
            OptCutsUVOptimizer._write_obj_with_uv(mesh, str(path), corners)
            text = path.read_text(encoding="utf-8")

        self.assertIn("f 3/3 1/1 4/4", text)
        self.assertIn("f 1/1 3/3 2/2", text)

    def test_output_face_reordering_is_aligned_to_frozen_domain(self):
        mesh = trimesh.Trimesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2], [0, 2, 3]]),
            process=False,
        )
        obj_text = """\
v 0 0 0
v 1 0 0
v 1 1 0
v 0 1 0
vt 0 0
vt 1 1
vt 0 1
vt 1 1
vt 1 0
vt 0 0
f 1/1 3/2 4/3
f 3/4 2/5 1/6
"""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "reordered.obj"
            path.write_text(obj_text, encoding="utf-8")
            parsed = OptCutsUVOptimizer._parse_obj_uv(str(path))
            aligned = OptCutsUVOptimizer._align_output_corners(mesh, parsed)

        expected = np.asarray(mesh.vertices[:, :2])[mesh.faces]
        np.testing.assert_array_equal(aligned, expected)

    def test_output_alignment_handles_cut_vertices_with_identical_coordinates(self):
        mesh = trimesh.Trimesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [-1.0, 0.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [3, 4, 5]]),
            process=False,
        )
        obj_text = """\
v 0 0 0
v 1 0 0
v 0 1 0
v 0 0 0
v 0 1 0
v -1 0 0
vt 10 10
vt 10 11
vt 9 10
vt 0 0
vt 1 0
vt 0 1
f 4/1 5/2 6/3
f 1/4 2/5 3/6
"""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "duplicated_seam_vertices.obj"
            path.write_text(obj_text, encoding="utf-8")
            parsed = OptCutsUVOptimizer._parse_obj_uv(str(path))
            aligned = OptCutsUVOptimizer._align_output_corners(mesh, parsed)

        np.testing.assert_array_equal(aligned[0], [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        np.testing.assert_array_equal(aligned[1], [[10.0, 10.0], [10.0, 11.0], [9.0, 10.0]])

    def test_residue_aware_run_writes_sidecar_and_confirms_binary_capability(self):
        mesh = trimesh.Trimesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2]]),
            process=False,
        )
        mesh.metadata["source_atom_indices"] = np.arange(3, dtype=np.int64)
        set_uv_layout(mesh, np.asarray(mesh.vertices[:, :2]))
        with tempfile.TemporaryDirectory() as tmp:
            binary = self._fake_optcuts_binary(tmp, advertise_footprints=True)
            optimizer = OptCutsUVOptimizer(
                replace(
                    OptCutsConfig(),
                    optcuts_bin=str(binary),
                    residue_fragmentation_weight=0.75,
                )
            )
            optimizer.optimize_patches(
                [mesh],
                initialization="provided",
                pack=False,
                build_report=False,
                source_residue_labels=["A:GLY:1"] * 3,
                residue_weights={"A:GLY:1": 2.0},
            )

        objective = mesh.metadata["optcuts_execution"]["residue_aware_objective"]
        self.assertTrue(objective["enabled"])
        self.assertTrue(objective["capability_confirmed"])
        self.assertEqual(objective["residue_count"], 1)
        self.assertEqual(objective["residue_fragmentation_weight"], 0.75)
        self.assertEqual(len(objective["sidecar_sha256"]), 64)
        self.assertEqual(objective["final_objective_weighted_fragmentation"], 0.0)

    def test_residue_aware_run_rejects_an_unpatched_binary(self):
        mesh = trimesh.Trimesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2]]),
            process=False,
        )
        mesh.metadata["source_atom_indices"] = np.arange(3, dtype=np.int64)
        set_uv_layout(mesh, np.asarray(mesh.vertices[:, :2]))
        with tempfile.TemporaryDirectory() as tmp:
            binary = self._fake_optcuts_binary(tmp, advertise_footprints=False)
            optimizer = OptCutsUVOptimizer(
                replace(
                    OptCutsConfig(),
                    optcuts_bin=str(binary),
                    residue_fragmentation_weight=1.0,
                )
            )
            with self.assertRaisesRegex(RuntimeError, "does not support residue-aware layouts"):
                optimizer.optimize_patches(
                    [mesh],
                    initialization="provided",
                    pack=False,
                    source_residue_labels=["A:GLY:1"] * 3,
                )

    @unittest.skipUnless(
        RESIDUE_AWARE_OPTCUTS_BIN.is_file(),
        "TopoPPI OptCuts binary is not available",
    )
    def test_interaction_binary_matches_zero_fragmentation_control(self):
        base = trimesh.Trimesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2], [0, 2, 3]]),
            process=False,
        )
        initialize_provenance(base, stage="zero_fragmentation_control")
        base.metadata["source_atom_indices"] = np.arange(4, dtype=np.int64)
        set_uv_layout(base, np.asarray(base.vertices[:, :2]))
        executions = []
        outputs = []
        for weight in (0.0, 0.75):
            mesh = base.copy()
            optimizer = OptCutsUVOptimizer(
                replace(
                    OptCutsConfig(),
                    optcuts_bin=str(RESIDUE_AWARE_OPTCUTS_BIN),
                    optcuts_mode=OptCutsConfig().optcuts_headless_mode,
                    residue_fragmentation_weight=weight,
                    timeout_sec=60.0,
                )
            )
            optimizer.optimize_patches(
                [mesh],
                initialization="provided",
                pack=False,
                build_report=False,
                source_residue_labels=["A:GLY:1"] * 4,
                residue_weights={"A:GLY:1": 2.0},
            )
            executions.append(mesh.metadata["optcuts_execution"])
            outputs.append(np.asarray(mesh.metadata["uv_optcuts"]))

        self.assertEqual(executions[0]["input_obj_sha256"], executions[1]["input_obj_sha256"])
        self.assertEqual(executions[0]["output_uv_checksum"], executions[1]["output_uv_checksum"])
        np.testing.assert_array_equal(outputs[0], outputs[1])
        self.assertEqual(
            executions[1]["residue_aware_objective"]["final_objective_weighted_fragmentation"],
            0.0,
        )

    @unittest.skipUnless(
        RESIDUE_AWARE_OPTCUTS_BIN.is_file(),
        "TopoPPI OptCuts binary is not available",
    )
    def test_real_interaction_binary_consumes_nonidentity_source_vertex_map(self):
        mesh = trimesh.Trimesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2], [0, 2, 3]]),
            process=False,
        )
        mesh.metadata["source_vertex_ids"] = np.array([10, 11, 12, 13], dtype=np.int64)
        mesh.metadata["source_atom_indices"] = np.arange(4, dtype=np.int64)
        set_uv_layout(mesh, np.asarray(mesh.vertices[:, :2]))
        optimizer = OptCutsUVOptimizer(
            replace(
                OptCutsConfig(),
                optcuts_bin=str(RESIDUE_AWARE_OPTCUTS_BIN),
                optcuts_mode=OptCutsConfig().optcuts_headless_mode,
                residue_fragmentation_weight=0.5,
                timeout_sec=60.0,
            )
        )

        optimizer.optimize_patches(
            [mesh],
            initialization="provided",
            pack=False,
            source_residue_labels=["A:GLY:1"] * 4,
            residue_weights={"A:GLY:1": 2.0},
        )

        self.assertEqual(mesh.metadata["uv_optcuts"].shape, (2, 3, 2))
        objective = mesh.metadata["optcuts_execution"]["residue_aware_objective"]
        self.assertTrue(objective["enabled"])
        self.assertEqual(objective["sidecar_schema_version"], 2)

    @unittest.skipUnless(OPTCUTS_BIN.is_file(), "OptCuts binary is not available")
    def test_real_binary_recovers_and_tracks_a_diskification_seam(self):
        mesh = trimesh.Trimesh(
            vertices=np.array(
                [
                    [-2.0, -2.0, 0.0],
                    [2.0, -2.0, 0.0],
                    [2.0, 2.0, 0.0],
                    [-2.0, 2.0, 0.0],
                    [-1.0, -1.0, 0.0],
                    [1.0, -1.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [-1.0, 1.0, 0.0],
                ]
            ),
            faces=np.array(
                [
                    [0, 1, 5],
                    [0, 5, 4],
                    [1, 2, 6],
                    [1, 6, 5],
                    [2, 3, 7],
                    [2, 7, 6],
                    [3, 0, 4],
                    [3, 4, 7],
                ]
            ),
            process=False,
        )
        initialize_provenance(mesh, stage="annulus")
        parameterizer = Parameterizer(ParameterizationConfig(min_face_area=1e-10))
        prepared, report = parameterizer.prepare_patch(mesh, return_info=True)
        self.assertIsNotNone(prepared, msg=report)
        self.assertGreater(report["diskification_added_vertex_count"], 0)
        optimizer = OptCutsUVOptimizer(
            replace(
                OptCutsConfig(),
                optcuts_bin=str(OPTCUTS_BIN),
                optcuts_mode=OptCutsConfig().optcuts_headless_mode,
                timeout_sec=60.0,
            )
        )

        optimizer.optimize_patches(
            [prepared],
            initialization="automatic",
            pack=False,
            build_report=False,
        )

        execution = prepared.metadata["optcuts_execution"]
        self.assertEqual(
            execution["input_geometry"]["collapsed_diskification_vertex_copy_count"],
            report["diskification_added_vertex_count"],
        )
        self.assertTrue(
            UVAtlasMetrics.parameterization_injectivity_stats(
                prepared,
                prepared.metadata["uv_optcuts"],
            )["globally_injective"]
        )
        self.assertRegex(execution["stdout_tail"], r"evaluate edge merge, [1-9]")
        constraint = execution["output_distortion_constraint"]
        self.assertTrue(constraint["satisfied"])
        self.assertLessEqual(
            constraint["energy"],
            constraint["bound"] + constraint["numeric_tolerance"],
        )

    @unittest.skipUnless(OPTCUTS_BIN.is_file(), "OptCuts binary is not available")
    def test_real_binary_consumes_provided_uv_and_returns_corner_uv(self):
        mesh = trimesh.Trimesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2], [0, 2, 3]]),
            process=False,
        )
        initial_uv = np.asarray(mesh.vertices[:, :2])
        set_uv_layout(mesh, initial_uv)
        config = replace(
            OptCutsConfig(),
            optcuts_bin=str(OPTCUTS_BIN),
            optcuts_mode=OptCutsConfig().optcuts_headless_mode,
            timeout_sec=60.0,
        )
        optimizer = OptCutsUVOptimizer(config)
        optimizer.optimize_patches([mesh], initialization="provided")

        self.assertEqual(mesh.metadata["uv_optcuts"].shape, (2, 3, 2))
        execution = mesh.metadata["optcuts_execution"]
        self.assertEqual(execution["initialization"], "provided_uv")
        self.assertIsNotNone(execution["initial_uv_checksum"])
        self.assertTrue(execution["per_corner_uv_preserved"])
        self.assertTrue(execution["upstream_reference"]["matches_packaged_linux_binary"])
