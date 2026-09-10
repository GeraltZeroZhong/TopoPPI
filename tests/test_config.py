import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

from topoppi.config import (
    DEFAULT_GUI_CONFIG,
    DEFAULT_RESIDUE_FRAGMENTATION_WEIGHT,
    DEFAULT_RUN_CONFIG,
    BenchmarkConfig,
    OptCutsConfig,
    TopologyConfig,
    TopoPPIRunConfig,
    VisualizationConfig,
)
from topoppi.errors import ConfigurationError


class ConfigTests(unittest.TestCase):
    def test_surface_camera_accepts_finite_angles_and_positive_zoom(self):
        VisualizationConfig(view="surface", surface_projection="perspective", surface_elevation=-25.,
                            surface_azimuth=450., surface_zoom=1.5, show_mesh=False).validate()
        for field, value in (
            ("view", "volume"), ("surface_projection", "flat"), ("show_mesh", "false"),
            ("surface_elevation", float("nan")), ("surface_azimuth", float("inf")),
            ("surface_zoom", 0.), ("surface_zoom", -1.), ("surface_zoom", float("inf")),
        ):
            with self.subTest(field=field, value=value), self.assertRaisesRegex(ConfigurationError, field):
                replace(VisualizationConfig(), **{field: value}).validate()
        with self.assertRaisesRegex(ConfigurationError, "footprint_color"):
            VisualizationConfig(view="surface", footprint_color="invalid-color").validate()

    def test_topoppi_defaults_to_the_complete_residue_aware_method(self):
        self.assertEqual(DEFAULT_RESIDUE_FRAGMENTATION_WEIGHT, 20.0)
        self.assertEqual(DEFAULT_RUN_CONFIG.optcuts.residue_fragmentation_weight, 20.0)
        self.assertEqual(TopoPPIRunConfig().optcuts.residue_fragmentation_weight, 20.0)
        self.assertEqual(OptCutsConfig().residue_fragmentation_weight, 0.0)
        self.assertEqual(DEFAULT_RUN_CONFIG.visualization.residue_scope, "interaction")
        self.assertTrue(DEFAULT_RUN_CONFIG.visualization.color_by_interaction_type)
        self.assertFalse(DEFAULT_RUN_CONFIG.visualization.use_geometric_interaction_fallback)
        self.assertEqual(DEFAULT_RUN_CONFIG.topology.distance_cutoff, 4.0)
        self.assertEqual(DEFAULT_RUN_CONFIG.visualization.min_points, 1)
        self.assertEqual(DEFAULT_GUI_CONFIG.default_patch_cutoff, 4.0)
        self.assertEqual(DEFAULT_GUI_CONFIG.default_min_points, 1)

    def test_config_validates_input_path_and_numeric_ranges(self):
        with tempfile.TemporaryDirectory() as tmp:
            pdb = Path(tmp) / "x.pdb"
            pdb.write_text("END\n", encoding="utf-8")
            config = TopoPPIRunConfig(pdb_file=str(pdb), chain_a="A", chain_b="B")
            config.validate()

            bad = TopoPPIRunConfig(
                pdb_file=str(pdb),
                chain_a="A",
                chain_b="B",
                topology=config.topology.__class__(distance_cutoff=0),
            )
            with self.assertRaises(ConfigurationError):
                bad.validate()

            with self.assertRaises(ConfigurationError):
                replace(config, surface=replace(config.surface, sigma=float("nan"))).validate()

            with self.assertRaises(ConfigurationError):
                replace(
                    config,
                    optcuts=replace(config.optcuts, optcuts_lambda_init=1.0),
                ).validate()

            with self.assertRaises(ConfigurationError):
                replace(
                    config,
                    optcuts=replace(config.optcuts, optcuts_quick_lambda_init=1.0),
                ).validate()

            with self.assertRaises(ConfigurationError):
                replace(
                    config,
                    optcuts=replace(config.optcuts, optcuts_distortion_bound=4.0),
                ).validate()

            with self.assertRaises(ConfigurationError):
                replace(
                    config,
                    optcuts=replace(config.optcuts, optcuts_quick_distortion_bound=3.9),
                ).validate()

            with self.assertRaises(ConfigurationError):
                replace(
                    config,
                    optcuts=replace(config.optcuts, optcuts_initial_cut_option=2),
                ).validate()

            with self.assertRaises(ConfigurationError):
                replace(
                    config,
                    optcuts=replace(config.optcuts, optcuts_method_type=1),
                ).validate()

            with self.assertRaises(ConfigurationError):
                replace(
                    config,
                    optcuts=replace(config.optcuts, residue_fragmentation_weight=-0.1),
                ).validate()

            with self.assertRaises(ConfigurationError):
                replace(config, surface=replace(config.surface, max_voxels=100.5)).validate()

            with self.assertRaises(ConfigurationError):
                replace(config, chain_b="A").validate()

            with self.assertRaises(ConfigurationError):
                replace(config, output_file=str(pdb)).validate()

            with self.assertRaises(ConfigurationError):
                replace(
                    config,
                    visualization=VisualizationConfig(residue_scope="all"),
                ).validate()

            with self.assertRaises(ConfigurationError):
                replace(
                    config,
                    visualization=VisualizationConfig(
                        color_by_interaction_type="yes"  # type: ignore[arg-type]
                    ),
                ).validate()

            with self.assertRaises(ConfigurationError):
                replace(
                    config,
                    visualization=VisualizationConfig(min_points=0),
                ).validate()

            with self.assertRaises(ConfigurationError):
                TopologyConfig(max_edge_face_incidence=3).validate()

            with self.assertRaises(ConfigurationError):
                replace(
                    config,
                    parameterization=replace(
                        config.parameterization,
                        expected_euler_characteristic=0,
                    ),
                ).validate()

            with self.assertRaises(ConfigurationError):
                replace(
                    config,
                    parameterization=replace(
                        config.parameterization,
                        expected_boundary_loops=2,
                    ),
                ).validate()

    def test_benchmark_config_rejects_string_booleans_and_fractional_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ConfigurationError):
                BenchmarkConfig(tmp, str(Path(tmp) / "out"), formal_mode="false").validate()
            with self.assertRaises(ConfigurationError):
                BenchmarkConfig(tmp, str(Path(tmp) / "out"), repetitions=3.5).validate()
            with self.assertRaises(ConfigurationError):
                BenchmarkConfig(
                    tmp,
                    str(Path(tmp) / "out"),
                    checkpoint_interval_structures=0,
                ).validate()

    def test_formal_worker_budget_covers_all_scheduled_optcuts_arms(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.csv"
            manifest.write_text("pdb,chain_a,chain_b\n", encoding="utf-8")
            common = {
                "input_folder": tmp,
                "output_root": str(Path(tmp) / "out"),
                "chain_selection_mode": "manifest",
                "manifest_path": str(manifest),
                "formal_mode": True,
                "benchmark_purpose": "quality",
                "repetitions": 1,
                "warmup_runs": 0,
                "optcuts_variants": (
                    "optcuts_automatic",
                    "optcuts_lscm_initialized",
                    "residue_aware_optcuts",
                ),
                "include_topology_ablation": False,
                "optcuts": OptCutsConfig(
                    expected_binary_sha256="a" * 64,
                    residue_fragmentation_weight=1.0,
                    timeout_sec=180.0,
                ),
            }
            BenchmarkConfig(**common, worker_timeout_sec=600.0).validate()
            with self.assertRaisesRegex(ConfigurationError, "sum of all scheduled"):
                BenchmarkConfig(**common, worker_timeout_sec=540.0).validate()
            with self.assertRaises(ConfigurationError):
                BenchmarkConfig(
                    tmp,
                    str(Path(tmp) / "out"),
                    worker_memory_limit_mb=-1.0,
                ).validate()
            with self.assertRaises(ConfigurationError):
                BenchmarkConfig(
                    tmp,
                    str(Path(tmp) / "out"),
                    report_filename="../outside.json",
                ).validate()
