"""Exercise 2D/3D view changes and saved-camera behavior through the desktop."""

import os
import sys
import tempfile
import tkinter as tk
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
import trimesh
from matplotlib.colors import to_rgba

from topoppi.atlas.uv import set_uv_layout
from topoppi.config import VisualizationConfig
from topoppi.gui_app.app import ProtSurfApp
from topoppi.io.io_loader import PDBLoader
from topoppi.visualization.atlas_io import load_atlas
from topoppi.visualization.visualizer import InterfaceVisualizer

FIXTURE = Path(__file__).parent / "fixtures" / "tiny_complex.pdb"


@unittest.skipUnless(os.environ.get("DISPLAY") or sys.platform in {"win32", "darwin"}, "Tk needs a graphical display")
class SurfaceDesktopTests(unittest.TestCase):
    def setUp(self):
        self.root = tk.Tk()
        self.root.withdraw()
        with mock.patch.object(ProtSurfApp, "_load_recent_items"):
            self.app = ProtSurfApp(self.root)
        self.app._remember_recent_file = lambda _path: None
        self.app._remember_recent_output_dir = lambda _path: None
        self.app.combo_map_style.set("Residue footprints")
        self.app.combo_residue_scope.set("Full patch context")
        self.app.var_color_type.set(False)
        self.app.var_avoid_overlap.set(False)
        self.app.var_auto_save.set(False)
        self.error_patch = mock.patch("topoppi.gui_app.ui_mixin.messagebox.showerror")
        self.errors = self.error_patch.start()
        loader = PDBLoader(FIXTURE)
        coords_a, atoms_a = loader.get_chain_data("A")
        coords_b, atoms_b = loader.get_chain_data("B")
        self.viz = InterfaceVisualizer(
            atoms_a, coords_a, coords_b, atoms_b, chain_a_id="A", chain_b_id="B",
            structure_label="tiny_complex",
            config=VisualizationConfig(map_style="footprints", residue_scope="patch", min_points=1),
        )
        self.patch = self.make_patch()

    def tearDown(self):
        self.error_patch.stop()
        if not self.app._closed:
            self.root.update_idletasks()
            for callback in self.root.tk.call("after", "info"):
                self.root.after_cancel(callback)
            self.app.close()

    @staticmethod
    def make_patch(offset=0.0):
        patch = trimesh.Trimesh(
            vertices=[[0., 0., offset], [1., 0., offset], [1., 1., offset], [0., 1., offset]],
            faces=[[0, 1, 2], [0, 2, 3]], process=False,
        )
        patch.metadata["source_atom_indices"] = np.array([0, 1, 4, 5])
        uv = patch.vertices[:, :2] + [offset * 2, 0]
        for key in ("uv", "uv_optcuts", "uv_global"):
            set_uv_layout(patch, uv, key=key)
        return patch

    def render(self, *, view="surface", patches=None):
        self.app.combo_view.set("3D interface" if view == "surface" else "2D atlas")
        success = self.app.update_plot(
            self.viz, patches or [self.patch], self.app.get_style_config(), complete_task=True,
            run_params={"path": str(FIXTURE), "chain_a": "A", "chain_b": "B", "min_points": 1,
                        "auto_save": False, "cutoff": 4., "res": 2., "sigma": 1.},
            run_manifest={"config": {}, "prolif_source": "none", "run_id": "surface-test"},
        )
        self.assertTrue(success, self.errors.call_args_list)
        return self.app.current_fig

    def test_view_switch_restores_all_patches_and_retains_the_camera_without_optimization(self):
        hidden = self.make_patch(2.)
        self.app.combo_map_style.set("Residue markers")
        with (
            mock.patch.object(self.viz, "count_patch_interaction_residues", side_effect=lambda p: int(p is self.patch)),
            mock.patch("topoppi.gui_app.workflow_mixin.OptCutsUVOptimizer") as optimizer,
        ):
            self.render(view="atlas", patches=[self.patch, hidden])
            self.assertEqual(len(self.app._successful_single_run["patches"]), 1)
            self.app.combo_view.set("3D interface")
            self.app._view_changed()
            self.assertEqual(len(self.app._successful_single_run["patches"]), 2)
            surface = self.app.current_fig._topoppi_surface
            np.testing.assert_allclose(
                surface["source_triangles"] @ surface["rotation"] + surface["centre"],
                np.concatenate([p.triangles for p in [self.patch, hidden]]), atol=1e-12,
            )
            surface["axis"].view_init(elev=32., azim=18.)
            self.app.on_mouse_release(SimpleNamespace())
            self.app.combo_view.set("2D atlas")
            self.app._view_changed()
            self.assertEqual(self.viz.last_style["view"], "atlas")
            self.assertEqual(len(self.app._successful_single_run["patches"]), 1)
            self.app.combo_view.set("3D interface")
            self.app._view_changed()
            restored = self.app.current_fig._topoppi_surface["axis"]
            self.assertAlmostEqual(restored.elev, 32.)
            self.assertAlmostEqual(restored.azim, 18.)
            optimizer.assert_not_called()
        self.errors.assert_not_called()

    def test_save_reopen_and_reset_use_the_current_camera(self):
        figure = self.render()
        surface = figure._topoppi_surface
        default_camera = surface["capture_camera"]()
        axis = surface["axis"]
        roll = {"roll": 19.} if hasattr(axis, "roll") else {}
        axis.view_init(elev=24., azim=41., **roll)
        axis.set_xlim3d(-.9, .7)
        # Save captures edits even before a mouse-release callback is delivered.
        with tempfile.TemporaryDirectory() as directory:
            path = str(Path(directory) / "surface.atlas.npz")
            with mock.patch("topoppi.gui_app.ui_mixin.filedialog.asksaveasfilename", return_value=path):
                self.app.save_atlas()
            document = load_atlas(path)
            self.assertEqual(document.style["view"], "surface")
            self.assertAlmostEqual(document.style["surface_camera"]["elevation"], 24.)
            if roll:
                self.assertAlmostEqual(document.style["surface_camera"]["roll"], 19.)
            np.testing.assert_allclose(document.style["surface_camera"]["limits"][0], [-.9, .7])
            self.app.reset_surface_view()
            self.assertAlmostEqual(axis.elev, default_camera["elevation"])
            with (
                mock.patch("topoppi.gui_app.ui_mixin.filedialog.askopenfilename", return_value=path),
                mock.patch("topoppi.gui_app.workflow_mixin.PDBLoader") as loader,
                mock.patch("topoppi.gui_app.workflow_mixin.OptCutsUVOptimizer") as optimizer,
            ):
                self.app.open_atlas()
            loader.assert_not_called()
            optimizer.assert_not_called()
        self.assertEqual(self.app.combo_view.get(), "3D interface")
        restored = self.app.current_fig._topoppi_surface["axis"]
        self.assertAlmostEqual(restored.elev, 24.)
        self.assertAlmostEqual(restored.azim, 41.)
        if roll:
            self.assertAlmostEqual(restored.roll, 19.)
        np.testing.assert_allclose(restored.get_xlim3d(), [-.9, .7])
        self.app.reset_surface_view()
        reset = self.app._successful_single_run["style"]["surface_camera"]
        self.assertAlmostEqual(reset["elevation"], default_camera["elevation"])
        self.assertAlmostEqual(reset["azimuth"], default_camera["azimuth"])
        np.testing.assert_allclose(reset["limits"], default_camera["limits"])
        self.errors.assert_not_called()

    def test_view_change_displays_a_corrected_pending_result_without_recalculation(self):
        self.app.combo_view.set("3D interface")
        self.app.var_highlight_residues.set("A:99999")
        self.app._run_style = self.app.get_style_config()
        self.app.accept_pipeline_result(
            self.viz, [self.patch], {"config": {}, "prolif_source": "none", "run_id": "pending-test"},
            {"path": str(FIXTURE), "chain_a": "A", "chain_b": "B", "min_points": 1,
             "auto_save": False, "cutoff": 4., "res": 2., "sigma": 1.},
        )
        self.assertIsNotNone(self.app._pending_single_run)
        self.assertIsNone(self.app._successful_single_run)
        self.app.var_highlight_residues.set("")
        with mock.patch("topoppi.gui_app.workflow_mixin.OptCutsUVOptimizer") as optimizer:
            self.app.combo_view.set("2D atlas")
            self.app._view_changed()
        optimizer.assert_not_called()
        self.assertIsNone(self.app._pending_single_run)
        self.assertEqual(self.app._successful_single_run["style"]["view"], "atlas")

    def test_rotation_does_not_open_color_dialog_and_double_click_recolors_only_one_residue(self):
        figure = self.render()
        surface = figure._topoppi_surface
        gid, objects = next(iter(self.viz.artist_map.items()))
        event = SimpleNamespace(inaxes=surface["axis"], xdata=0., ydata=0., button=1, dblclick=False)
        with mock.patch("topoppi.gui_app.plot_mixin.colorchooser.askcolor", return_value=(None, "#336699")) as chooser:
            self.app.on_pick(SimpleNamespace(artist=objects["collection"]))
            self.app.on_mouse_press(event)
            self.app.on_mouse_move(event)
            self.app.on_mouse_release(event)
            chooser.assert_not_called()
            self.assertIsNone(self.app._drag_state)
            event.dblclick = True
            with (
                mock.patch.dict(surface, {"pick_residue": lambda _event: gid}),
                mock.patch.object(objects["collection"], "set_facecolor") as paint_whole_surface,
            ):
                self.app.on_mouse_press(event)
            paint_whole_surface.assert_not_called()
            chooser.assert_called_once()
        self.assertEqual(self.app.residue_color_overrides, {objects["residue_key"]: "#336699"})
        repainted = self.app.current_fig._topoppi_surface
        selected = repainted["cell_gids"] == gid
        self.assertTrue(selected.any())
        np.testing.assert_allclose(repainted["base_colors"][selected], np.tile(to_rgba("#336699"), (selected.sum(), 1)))
        np.testing.assert_allclose(
            repainted["base_colors"][~selected], np.tile(to_rgba("#DCE8EF"), ((~selected).sum(), 1))
        )

    def test_surface_markers_share_numeric_annotations_and_view_controls(self):
        self.app.combo_map_style.set("Residue markers")
        self.app.annotation_values = {"A:GLY:1": .4, "A:ALA:2": None}
        self.app.var_highlight_residues.set("A:GLY:1")
        self.app.var_footprint_labels.set("highlighted")
        self.app.var_show_mesh.set(False)
        self.app.combo_surface_projection.set("Perspective")
        figure = self.render()
        self.assertTrue(self.app.surface_controls.grid_info())
        self.assertTrue(self.app.footprint_controls.grid_info())
        self.assertEqual(self.viz.last_style["annotation_values"], self.app.annotation_values)
        self.assertEqual(self.viz.last_style["surface_projection"], "perspective")
        self.assertFalse(self.viz.last_style["show_mesh"])
        self.assertFalse(any(c.get_gid() == "mesh_segments" for c in figure._topoppi_surface["axis"].collections))
        with mock.patch("topoppi.gui_app.plot_mixin.colorchooser.askcolor") as chooser:
            self.app.on_mouse_press(SimpleNamespace(inaxes=figure._topoppi_surface["axis"], button=1, dblclick=True))
        chooser.assert_not_called()
        self.app.combo_map_style.set("Residue footprints")
        self.app._map_style_changed()
        self.assertEqual(self.viz.last_style["annotation_values"], self.app.annotation_values)
        self.assertEqual(self.viz.last_style["footprint_labels"], "highlighted")
        self.assertEqual(self.viz.last_style["view"], "surface")
        self.errors.assert_not_called()

    def test_surface_recolor_replaces_prior_piece_colors_in_the_2d_atlas(self):
        self.app.combo_map_style.set("Residue markers")
        corners = self.patch.triangles[:, :, :2].copy()
        corners[1, :, 0] += 2.
        for key in ("uv", "uv_optcuts", "uv_global"):
            set_uv_layout(self.patch, corners, key=key)
        with mock.patch.object(self.viz, "count_patch_interaction_residues", return_value=2):
            self.render(view="atlas")
            piece_ids = [uid for uid, objects in self.viz.artist_map.items() if objects["residue_key"] == "A:GLY:1"]
            other_id = next(uid for uid, objects in self.viz.artist_map.items() if objects["residue_key"] == "A:ALA:2")
            self.assertEqual(len(piece_ids), 2)
            self.app.marker_color_overrides = {uid: "#CC5500" for uid in piece_ids + [other_id]}
            self.app.redraw_plot()
            self.app.combo_view.set("3D interface")
            self.app._view_changed()
            surface = self.app.current_fig._topoppi_surface
            gid = next(uid for uid, objects in self.viz.artist_map.items() if objects["residue_key"] == "A:GLY:1")
            with (
                mock.patch.dict(surface, {"pick_residue": lambda _event: gid}),
                mock.patch("topoppi.gui_app.plot_mixin.colorchooser.askcolor", return_value=(None, "#338855")),
            ):
                self.app.on_mouse_press(SimpleNamespace(inaxes=surface["axis"], button=1, dblclick=True))
            self.assertEqual(self.app.marker_color_overrides, {other_id: "#CC5500"})
            self.app.combo_view.set("2D atlas")
            self.app._view_changed()
        for uid in piece_ids:
            np.testing.assert_allclose(self.viz.artist_map[uid]["scatter"].get_facecolor()[0], to_rgba("#338855"))
        np.testing.assert_allclose(self.viz.artist_map[other_id]["scatter"].get_facecolor()[0], to_rgba("#CC5500"))
        self.errors.assert_not_called()

    def test_new_run_receives_surface_settings_and_starts_with_a_fresh_camera(self):
        self.render()
        self.app.current_fig._topoppi_surface["axis"].view_init(elev=22., azim=37.)
        self.app.on_mouse_release(SimpleNamespace())
        self.app.var_input_path.set(str(FIXTURE))
        self.app.var_chain_a.set("A")
        self.app.var_chain_b.set("B")
        self.app.combo_surface_projection.set("Perspective")
        self.app.var_show_mesh.set(False)
        with mock.patch("topoppi.gui_app.ui_mixin.threading.Thread") as thread:
            thread.return_value.is_alive.return_value = False
            self.app.start_analysis()
        self.app._finish_task()
        _params, config = thread.call_args.kwargs["args"]
        self.assertEqual(config.visualization.view, "surface")
        self.assertEqual(config.visualization.surface_projection, "perspective")
        self.assertFalse(config.visualization.show_mesh)
        self.assertNotIn("surface_camera", self.app._run_style)
        self.errors.assert_not_called()
