"""CLI controls for computation-free 3D interface rendering."""

import io
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest import mock

import numpy as np
from test_atlas_io import make_atlas

from topoppi import cli
from topoppi.atlas.uv import as_corner_uv
from topoppi.visualization.atlas_io import load_atlas, save_atlas
from topoppi.visualization.visualizer import InterfaceVisualizer


class SurfaceParserTests(unittest.TestCase):
    def test_surface_view_and_camera_reach_new_run_configuration(self):
        with mock.patch("topoppi.cli.run_interface_mapping") as run:
            code = cli.main([
                "input.pdb", "--view", "surface", "--projection", "perspective", "--elevation", "45",
                "--azimuth", "-35", "--zoom", "1.6", "--no-mesh", "--map-style", "footprints",
            ])
        self.assertEqual(code, 0)
        config = run.call_args.args[0].visualization
        self.assertEqual(config.view, "surface")
        self.assertEqual(config.surface_projection, "perspective")
        self.assertEqual((config.surface_elevation, config.surface_azimuth, config.surface_zoom), (45., -35., 1.6))
        self.assertFalse(config.show_mesh)
        self.assertEqual(config.map_style, "footprints")
        self.assertEqual(config.residue_scope, "patch")

    def test_view_and_residue_representation_are_independent(self):
        for options, expected_scope in (([], "patch"), (["--residue-scope", "interaction"], "interaction")):
            with self.subTest(options=options), mock.patch("topoppi.cli.run_interface_mapping") as run:
                self.assertEqual(cli.main(["input.pdb", "--view", "surface", *options]), 0)
                config = run.call_args.args[0].visualization
                self.assertEqual(config.map_style, "markers")
                self.assertEqual(config.residue_scope, expected_scope)

    def test_nonfinite_angles_and_invalid_zoom_are_actionable_parser_errors(self):
        for prefix in (["input.pdb"], ["render", "saved.npz", "-o", "map.png"]):
            for option, value in (("--zoom", "0"), ("--zoom", "-1"), ("--zoom", "inf"),
                                  ("--elevation", "nan"), ("--azimuth", "inf")):
                with self.subTest(prefix=prefix, option=option, value=value):
                    error = io.StringIO()
                    with redirect_stderr(error), self.assertRaisesRegex(SystemExit, "2"):
                        cli.main([*prefix, option, value])
                    self.assertIn(option, error.getvalue())
                    self.assertIn("finite number", error.getvalue())


class SurfaceRenderTests(unittest.TestCase):
    def test_saved_marker_atlas_exports_surface_without_geometry_recomputation(self):
        patches, viz = make_atlas(style="markers")
        viz.interaction_partner_map = {}
        with tempfile.TemporaryDirectory() as tmp:
            source, target = Path(tmp) / "source.npz", Path(tmp) / "surface.npz"
            output = Path(tmp) / "surface.pdf"
            save_atlas(source, patches, viz, style_config={"residue_scope": "interaction"})
            with mock.patch("topoppi.cli.run_interface_mapping", side_effect=AssertionError("Recomputed geometry")):
                code = cli.main(["render", str(source), "--view", "surface", "-o", str(output),
                                 "--export-atlas", str(target)])
            self.assertEqual(code, 0)
            self.assertTrue(output.read_bytes().startswith(b"%PDF"))
            restored = load_atlas(target)
        self.assertEqual(restored.style["view"], "surface")
        self.assertEqual(restored.style["map_style"], "markers")
        self.assertEqual(restored.style["residue_scope"], "patch")
        np.testing.assert_array_equal(restored.patches[0].vertices, patches[0].vertices)
        np.testing.assert_array_equal(as_corner_uv(restored.patches[0]), as_corner_uv(patches[0]))

    def test_saved_surface_camera_and_mesh_survive_unspecified_render_options(self):
        patches, viz = make_atlas()
        with tempfile.TemporaryDirectory() as tmp:
            source, first, second = (Path(tmp) / name for name in ("source.npz", "first.npz", "second.npz"))
            output = Path(tmp) / "surface.png"
            save_atlas(source, patches, viz)
            self.assertEqual(cli.main([
                "render", str(source), "--view", "surface", "--no-mesh", "--projection", "perspective",
                "--elevation", "28", "--azimuth", "-65", "--zoom", "1.2", "-o", str(output),
                "--export-atlas", str(first),
            ]), 0)
            original = load_atlas(first)
            self.assertEqual(cli.main(["render", str(first), "-o", str(output), "--export-atlas", str(second)]), 0)
            restored = load_atlas(second)
        for field in ("view", "show_mesh", "surface_projection", "surface_elevation", "surface_azimuth",
                      "surface_zoom", "surface_camera"):
            self.assertEqual(restored.style[field], original.style[field])
        self.assertFalse(restored.style["show_mesh"])

    def test_explicit_camera_options_replace_saved_camera_and_reenable_mesh(self):
        patches, viz = make_atlas()
        captured = []
        plot = InterfaceVisualizer.plot_patches

        def capture(instance, *args, **kwargs):
            captured.append(dict(kwargs["style_config"]))
            return plot(instance, *args, **kwargs)

        with tempfile.TemporaryDirectory() as tmp:
            source, target = Path(tmp) / "source.npz", Path(tmp) / "target.npz"
            output = Path(tmp) / "surface.png"
            save_atlas(source, patches, viz)
            self.assertEqual(cli.main(["render", str(source), "--view", "surface", "--no-mesh",
                                       "--hide-seams", "--hide-residue-borders",
                                       "-o", str(output), "--export-atlas", str(source)]), 0)
            self.assertIn("surface_camera", load_atlas(source).style)
            with mock.patch.object(InterfaceVisualizer, "plot_patches", capture):
                self.assertEqual(cli.main(["render", str(source), "--elevation", "15", "--show-mesh",
                                           "--show-seams", "--show-residue-borders",
                                           "-o", str(output), "--export-atlas", str(target)]), 0)
            updated = load_atlas(target)
        self.assertNotIn("surface_camera", captured[0])
        self.assertTrue(updated.style["show_mesh"])
        self.assertTrue(updated.style["show_seams"])
        self.assertTrue(updated.style["show_residue_borders"])
        self.assertEqual(updated.style["surface_elevation"], 15.)


if __name__ == "__main__":
    unittest.main()
