"""Preserve interactive cameras and the physical arrangement of surface patches."""

import tempfile
import unittest
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from test_atlas_io import make_atlas

from topoppi.atlas.uv import as_corner_uv
from topoppi.visualization.atlas_io import load_atlas, save_atlas


class SurfaceCameraTests(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_arcball_roll_and_pan_survive_save_reopen_and_atlas_view(self):
        patches, viz = make_atlas()
        figure = viz.plot_patches(patches, style_config={"view": "surface"}, show=False)
        axis = figure._topoppi_surface["axis"]
        if not hasattr(axis, "roll"):
            self.skipTest("This Matplotlib version uses the elevation/azimuth camera.")
        axis.view_init(elev=35., azim=20., roll=40.)
        axis.set_xlim3d(-.3, .7)
        axis.set_ylim3d(-.6, .4)
        axis.set_zlim3d(-.45, .55)
        figure.canvas.draw()

        with tempfile.TemporaryDirectory() as tmp:
            saved = Path(tmp) / "rotated.npz"
            save_atlas(saved, patches, viz)
            restored = load_atlas(saved)
            figure = restored.visualizer.plot_patches(
                restored.patches, style_config=restored.style, show=False,
            )
            restored_axis = figure._topoppi_surface["axis"]
            self.assertEqual((restored_axis.elev, restored_axis.azim, restored_axis.roll), (35., 20., 40.))
            np.testing.assert_allclose(restored_axis.get_xlim3d(), [-.3, .7])
            np.testing.assert_allclose(restored_axis.get_ylim3d(), [-.6, .4])
            np.testing.assert_allclose(restored_axis.get_zlim3d(), [-.45, .55])

            # A saved 2D view carries its previous 3D camera for later reopening.
            restored.visualizer.plot_patches(
                restored.patches, style_config={**restored.visualizer.last_style, "view": "atlas"}, show=False,
            )
            save_atlas(saved, restored.patches, restored.visualizer)
            atlas_view = load_atlas(saved)
            figure = atlas_view.visualizer.plot_patches(
                atlas_view.patches, style_config={**atlas_view.style, "view": "surface"}, show=False,
            )
            self.assertEqual(figure._topoppi_surface["axis"].roll, 40.)

    def test_multiple_patches_keep_their_physical_separation_and_saved_geometry(self):
        patches, viz = make_atlas()
        first = patches[0]
        second = first.copy()
        separation = np.asarray([17., 4., 2.])
        second.vertices += separation
        original_uv = [as_corner_uv(first).copy(), as_corner_uv(second).copy()]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "multipatch.npz"
            save_atlas(path, [first, second], viz)
            restored = load_atlas(path)
            figure = restored.visualizer.plot_patches(
                restored.patches, style_config={**restored.style, "view": "surface"}, show=False,
            )
        triangles = figure._topoppi_surface["source_triangles"]
        face_count = len(first.faces)
        centre_a = triangles[:face_count].mean(axis=(0, 1))
        centre_b = triangles[face_count:].mean(axis=(0, 1))
        self.assertAlmostEqual(np.linalg.norm(centre_b - centre_a), np.linalg.norm(separation))
        np.testing.assert_array_equal(restored.patches[0].vertices, first.vertices)
        np.testing.assert_array_equal(restored.patches[1].vertices, second.vertices)
        for patch, uv in zip(restored.patches, original_uv, strict=True):
            np.testing.assert_array_equal(as_corner_uv(patch), uv)


if __name__ == "__main__":
    unittest.main()
