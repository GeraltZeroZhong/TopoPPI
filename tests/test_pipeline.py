import logging
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from topoppi.config import DEFAULT_RUN_CONFIG
from topoppi.errors import InputDataError, PipelineError
from topoppi.file_utils import sha256_file
from topoppi.io.io_loader import PDBLoader
from topoppi.pipeline import (
    _build_interaction_partner_map,
    _extract_patches,
    _load_chain_data,
    _prepare_output_path,
    _render_output,
    _resolve_prolif_file,
    layout_interface,
    run_interface_mapping,
)

FIXTURES = Path(__file__).parent / "fixtures"


class PipelineTests(unittest.TestCase):
    def test_shared_layout_defaults_are_standard_and_explicit_controls_are_preserved(self):
        from topoppi.config import OptCutsConfig

        prepared = SimpleNamespace(source={"path": "input.pdb"},
            partners=(SimpleNamespace(chains=("A",), layout_patches=lambda: []), SimpleNamespace(chains=("B",))),
            chain_data=lambda side: ([], [], [], []), contacts=lambda side: {})
        explicit = OptCutsConfig(residue_fragmentation_weight=0.0, optcuts_mode=10)
        for selected in (None, explicit):
            with (
                self.subTest(explicit=selected is not None),
                mock.patch("topoppi.pipeline.OptCutsUVOptimizer") as optimizer,
                mock.patch("topoppi.pipeline._parameterize_patches", return_value=[]),
                mock.patch("topoppi.pipeline._optimize_patches", return_value={}),
            ):
                layout_interface(prepared, optcuts=selected)
                effective = optimizer.call_args.args[0]
                if selected is None:
                    self.assertEqual(effective.residue_fragmentation_weight, 20.0)
                    self.assertEqual(effective.optcuts_mode, effective.optcuts_headless_mode)
                else:
                    self.assertIs(effective, explicit)

    def test_explicit_geometric_source_skips_prolif_generation(self):
        config = replace(DEFAULT_RUN_CONFIG, interaction_source="geometric")
        with mock.patch("topoppi.pipeline.generate_prolif_interactions") as generate:
            self.assertIsNone(_resolve_prolif_file(config, logging.getLogger(), input_sha256="unused"))
        generate.assert_not_called()
        loader = PDBLoader(FIXTURES / "tiny_complex.pdb")
        coords_a, atoms_a = loader.get_chain_data("A")
        coords_b, atoms_b = loader.get_chain_data("B")
        partners, source = _build_interaction_partner_map(
            (coords_a, atoms_a, coords_b, atoms_b), None, config, logging.getLogger(), input_sha256="unused",
        )
        self.assertEqual(source, "geometric")
        self.assertTrue(partners)

    def test_footprint_mode_retains_patches_without_interaction_markers(self):
        config = replace(DEFAULT_RUN_CONFIG, visualization=replace(DEFAULT_RUN_CONFIG.visualization,
                                                                  map_style="footprints", min_points=10))
        visualizer, patch = mock.Mock(), mock.Mock()
        visualizer.count_patch_interaction_residues.return_value = 0
        visualizer.last_report = {"displayed_residue_count": 3}
        with mock.patch("topoppi.pipeline.InterfaceVisualizer", return_value=visualizer), mock.patch("topoppi.pipeline.plt.close"):
            report = _render_output([patch], ([], [], [], []), None, config, logging.getLogger(),
                                    interaction_partner_map={}, interaction_source="geometric")
        self.assertEqual(visualizer.plot_patches.call_args.args[0], [patch])
        self.assertEqual(report["display_filter"]["hidden_patch_count"], 0)
        self.assertEqual(report["display_filter"]["policy"], "complete_footprints")

    def test_output_path_is_prepared_before_structure_work(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "new" / "results" / "interface.TIFF"
            config = replace(
                DEFAULT_RUN_CONFIG,
                pdb_file=str(FIXTURES / "tiny_complex.pdb"),
                output_file=str(output_path),
            )

            def stop_after_output_preparation(*_args):
                self.assertTrue(output_path.parent.is_dir())
                raise InputDataError("stop after output preparation")

            with (
                mock.patch(
                    "topoppi.pipeline._load_chain_data",
                    side_effect=stop_after_output_preparation,
                ),
                mock.patch("topoppi.pipeline.OptCutsUVOptimizer") as optimizer,
                mock.patch("topoppi.pipeline.generate_prolif_interactions") as generate_prolif,
            ):
                with self.assertRaisesRegex(InputDataError, "stop after output preparation"):
                    run_interface_mapping(config)

        optimizer.assert_not_called()
        generate_prolif.assert_not_called()

    def test_unsupported_output_extension_stops_before_structure_work(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            for filename in ("interface.jpg", "interface"):
                with self.subTest(filename=filename):
                    config = replace(
                        DEFAULT_RUN_CONFIG,
                        pdb_file=str(FIXTURES / "tiny_complex.pdb"),
                        output_file=str(Path(temp_dir) / filename),
                    )
                    with (
                        mock.patch("topoppi.pipeline._load_chain_data") as load_chain_data,
                        mock.patch("topoppi.pipeline.OptCutsUVOptimizer") as optimizer,
                        mock.patch("topoppi.pipeline.generate_prolif_interactions") as generate_prolif,
                    ):
                        with self.assertRaisesRegex(
                            PipelineError,
                            r"Choose a \.png, \.tif, \.tiff, \.svg, or \.pdf file",
                        ):
                            run_interface_mapping(config)

                    load_chain_data.assert_not_called()
                    optimizer.assert_not_called()
                    generate_prolif.assert_not_called()

    def test_output_parent_failure_is_actionable(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            parent_file = Path(temp_dir) / "occupied"
            parent_file.write_text("file", encoding="utf-8")

            with self.assertRaisesRegex(
                PipelineError,
                "Choose another output path or create the directory first",
            ):
                _prepare_output_path(str(parent_file / "interface.png"))

    def test_missing_chain_error_lists_available_protein_chains(self):
        config = replace(
            DEFAULT_RUN_CONFIG,
            pdb_file=str(FIXTURES / "tiny_complex.pdb"),
            chain_a="Z",
        )

        with self.assertRaisesRegex(
            InputDataError,
            r"Selected chain\(s\) Z were not found\. Available protein chains: A, B\.",
        ):
            _load_chain_data(config, logging.getLogger())

    def test_missing_chain_stops_before_prolif_and_optcuts(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            config = replace(
                DEFAULT_RUN_CONFIG,
                pdb_file=str(FIXTURES / "tiny_complex.pdb"),
                chain_a="Z",
                output_file=str(Path(temp_dir) / "interface.png"),
            )
            with (
                mock.patch("topoppi.pipeline.OptCutsUVOptimizer") as optimizer,
                mock.patch("topoppi.pipeline.sha256_file") as input_digest,
                mock.patch("topoppi.pipeline.generate_prolif_interactions") as generate_prolif,
            ):
                with self.assertRaisesRegex(
                    InputDataError,
                    r"Available protein chains: A, B",
                ):
                    run_interface_mapping(config)

        optimizer.assert_not_called()
        input_digest.assert_not_called()
        generate_prolif.assert_not_called()

    def test_no_interface_error_suggests_a_cutoff_from_the_nearest_face(self):
        manager = mock.Mock()
        manager.get_interface_patches.return_value = []
        manager.last_report = {
            "status": "no_interface_faces",
            "nearest_partner_distance_angstrom": 6.25,
        }
        with mock.patch("topoppi.pipeline.TopologyManager", return_value=manager):
            with self.assertRaisesRegex(
                PipelineError,
                r"nearest face was 6\.25 Å away.*increase --cutoff",
            ):
                _extract_patches(
                    mock.Mock(),
                    mock.Mock(),
                    DEFAULT_RUN_CONFIG,
                    logging.getLogger(),
                )

    def test_prolif_pairs_are_the_authoritative_single_run_interaction_map(self):
        loader = PDBLoader(FIXTURES / "tiny_complex.pdb")
        coords_a, atoms_a = loader.get_chain_data("A")
        coords_b, atoms_b = loader.get_chain_data("B")
        config = replace(
            DEFAULT_RUN_CONFIG,
            pdb_file=str(FIXTURES / "tiny_complex.pdb"),
        )

        partners, source = _build_interaction_partner_map(
            (coords_a, atoms_a, coords_b, atoms_b),
            str(FIXTURES / "prolif_interactions.json"),
            config,
            logging.getLogger(),
            input_sha256=sha256_file(config.pdb_file),
        )

        self.assertEqual(source, "prolif")
        self.assertEqual(partners, {"A:GLY:1": {"B:ALA:2": 1}})

    def test_geometry_is_used_only_for_explicit_interaction_fallback(self):
        config = replace(
            DEFAULT_RUN_CONFIG,
            visualization=replace(
                DEFAULT_RUN_CONFIG.visualization,
                use_geometric_interaction_fallback=True,
            ),
        )
        expected = {"A:GLY:1": {"B:GLY:1": 2}}
        with mock.patch(
            "topoppi.pipeline.geometric_contact_partner_map",
            return_value=expected,
        ) as geometric_map:
            partners, source = _build_interaction_partner_map(
                ([], [], [], []),
                None,
                config,
                logging.getLogger(),
                input_sha256="a" * 64,
            )

        self.assertEqual(source, "geometric_fallback")
        self.assertIs(partners, expected)
        geometric_map.assert_called_once()

    def test_prolif_generation_failure_is_reported_as_a_pipeline_error(self):
        with mock.patch(
            "topoppi.pipeline.generate_prolif_interactions",
            side_effect=RuntimeError("unsupported bond topology"),
        ):
            with self.assertRaisesRegex(PipelineError, "ProLIF interaction generation failed"):
                _resolve_prolif_file(
                    DEFAULT_RUN_CONFIG,
                    logging.getLogger(),
                    input_sha256="a" * 64,
                )

    def test_generated_prolif_is_written_beside_the_output_image(self):
        config = replace(DEFAULT_RUN_CONFIG, output_file="results/interface.png")
        with mock.patch(
            "topoppi.pipeline.generate_prolif_interactions",
            return_value="results/complex.A-B.prolif.json",
        ) as generate:
            result = _resolve_prolif_file(
                config,
                logging.getLogger(),
                input_sha256="a" * 64,
            )

        self.assertEqual(result, "results/complex.A-B.prolif.json")
        self.assertEqual(generate.call_args.kwargs["output_dir"], "results")

    def test_geometric_fallback_requires_explicit_configuration(self):
        config = replace(
            DEFAULT_RUN_CONFIG,
            visualization=replace(
                DEFAULT_RUN_CONFIG.visualization,
                use_geometric_interaction_fallback=True,
            ),
        )
        with mock.patch(
            "topoppi.pipeline.generate_prolif_interactions",
            side_effect=RuntimeError("missing ProLIF stack"),
        ):
            self.assertIsNone(
                _resolve_prolif_file(
                    config,
                    logging.getLogger(),
                    input_sha256="a" * 64,
                )
            )

    def test_run_manifest_report_keeps_resolved_optcuts_provenance(self):
        artifact = {
            "requested": "OptCuts_bin",
            "resolved": "/opt/topoppi/OptCuts_bin",
            "sha256": "a" * 64,
            "env_var": "TOPOPPI_OPTCUTS_BIN",
        }
        interaction_partners = {"A:GLY:1": {"B:ALA:2": 3}}
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = Path(temp_dir) / "complex.pdb"
            input_path.write_text("END\n", encoding="utf-8")
            config = replace(
                DEFAULT_RUN_CONFIG,
                pdb_file=str(input_path),
                output_file=str(Path(temp_dir) / "interface.png"),
            )
            optimizer = mock.Mock()
            optimizer.preflight_binary.return_value = artifact
            patch = mock.Mock()
            chain_data = ([], [], [], [])
            prepared = mock.Mock()
            prepared.partners = [mock.Mock(surface_report={}, topology_report={})]
            prepared.partners[0].layout_patches.return_value = [patch]
            prepared.metadata.return_value = {"analysis_id": "test"}
            with (
                mock.patch("topoppi.pipeline.OptCutsUVOptimizer", return_value=optimizer),
                mock.patch(
                    "topoppi.pipeline._resolve_prolif_file",
                    return_value="generated.prolif.json",
                ) as resolve_prolif_file,
                mock.patch("topoppi.pipeline._load_chain_data", return_value=chain_data),
                mock.patch(
                    "topoppi.pipeline._build_interaction_partner_map",
                    return_value=(interaction_partners, "prolif"),
                ) as build_interaction_partner_map,
                mock.patch("topoppi.pipeline.prepare_interface", return_value=prepared),
                mock.patch("topoppi.pipeline._parameterize_patches", return_value=[patch]),
                mock.patch(
                    "topoppi.pipeline._optimize_patches",
                    return_value={"status": "ok"},
                ) as optimize_patches,
                mock.patch(
                    "topoppi.pipeline._render_output",
                    return_value={"status": "ok", "residue_scope": "interaction"},
                ) as render_output,
                mock.patch("topoppi.pipeline._write_run_manifest"),
                mock.patch(
                    "topoppi.pipeline.sha256_file",
                    return_value="b" * 64,
                ) as input_digest,
            ):
                result = run_interface_mapping(config)

        input_digest.assert_called_once_with(str(input_path))
        self.assertEqual(resolve_prolif_file.call_args.kwargs["input_sha256"], "b" * 64)
        self.assertEqual(
            build_interaction_partner_map.call_args.kwargs["input_sha256"],
            "b" * 64,
        )
        self.assertEqual(result.input_sha256, "b" * 64)
        self.assertEqual(result.optimizer_report["optcuts_resolved"], artifact)
        self.assertEqual(result.visualization["residue_scope"], "interaction")
        self.assertIs(
            optimize_patches.call_args.kwargs["interaction_partner_map"],
            interaction_partners,
        )
        self.assertEqual(optimize_patches.call_args.kwargs["interaction_source"], "prolif")
        self.assertIs(
            render_output.call_args.kwargs["interaction_partner_map"],
            interaction_partners,
        )

    def test_render_output_closes_the_saved_headless_figure(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            config = replace(
                DEFAULT_RUN_CONFIG,
                output_file=str(Path(temp_dir) / "interface.png"),
                visualization=replace(DEFAULT_RUN_CONFIG.visualization, min_points=10),
            )
            chain_data = ([], [], [], [])
            figure = mock.Mock()
            visualizer = mock.Mock()
            patch = mock.Mock()
            visualizer.count_patch_interaction_residues.return_value = 10
            visualizer.plot_patches.return_value = figure
            visualizer.last_report = {
                "status": "ok",
                "residue_scope": "interaction",
                "patch_interaction_residue_count": 2,
            }
            with (
                mock.patch("topoppi.pipeline.InterfaceVisualizer", return_value=visualizer),
                mock.patch("topoppi.pipeline.plt.close") as close,
            ):
                report = _render_output(
                    [patch],
                    chain_data,
                    None,
                    config,
                    logging.getLogger(),
                    interaction_partner_map={},
                    interaction_source="geometric_fallback",
                )

        close.assert_called_once_with(figure)
        visualizer.plot_patches.assert_called_once()
        self.assertEqual(visualizer.plot_patches.call_args.args[0], [patch])
        self.assertEqual(report["patch_interaction_residue_count"], 2)
        self.assertEqual(report["interaction_residue_source"], "geometric_fallback")
        self.assertEqual(
            report["display_filter"],
            {
                "policy": "interaction_threshold",
                "min_points": 10,
                "optimized_patch_count": 1,
                "displayed_patch_count": 1,
                "hidden_patch_count": 0,
                "patch_interaction_residue_counts": [10],
            },
        )

    def test_render_output_filters_patches_after_optimization(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            config = replace(
                DEFAULT_RUN_CONFIG,
                output_file=str(Path(temp_dir) / "interface.png"),
                visualization=replace(DEFAULT_RUN_CONFIG.visualization, min_points=10),
            )
            chain_data = ([], [], [], [])
            kept_patch = mock.Mock()
            hidden_patch = mock.Mock()
            visualizer = mock.Mock()
            visualizer.count_patch_interaction_residues.side_effect = [12, 4]
            visualizer.plot_patches.return_value = None
            visualizer.last_report = {"status": "ok"}
            with mock.patch("topoppi.pipeline.InterfaceVisualizer", return_value=visualizer):
                report = _render_output(
                    [kept_patch, hidden_patch],
                    chain_data,
                    None,
                    config,
                    logging.getLogger(),
                    interaction_partner_map={},
                    interaction_source="prolif",
                )

        self.assertEqual(visualizer.plot_patches.call_args.args[0], [kept_patch])
        self.assertEqual(report["display_filter"]["displayed_patch_count"], 1)
        self.assertEqual(report["display_filter"]["hidden_patch_count"], 1)
        self.assertEqual(report["display_filter"]["patch_interaction_residue_counts"], [12, 4])

    def test_render_output_translates_expected_visualization_errors(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            config = replace(
                DEFAULT_RUN_CONFIG,
                output_file=str(Path(temp_dir) / "interface.png"),
            )
            chain_data = ([], [], [], [])
            with mock.patch(
                "topoppi.pipeline.InterfaceVisualizer",
                side_effect=ValueError("ProLIF chain_b mismatch: C != B"),
            ):
                with self.assertRaisesRegex(PipelineError, "Visualization failed"):
                    _render_output(
                        [],
                        chain_data,
                        "interactions.json",
                        config,
                        logging.getLogger(),
                        interaction_partner_map={},
                        interaction_source="prolif",
                    )


if __name__ == "__main__":
    unittest.main()
