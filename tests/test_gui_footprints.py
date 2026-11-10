"""GUI actions that must preserve editable geometry, styles, and annotations."""

import os
import subprocess
import sys
import tempfile
import threading
import time
import tkinter as tk
import unittest
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
import trimesh
from matplotlib.colors import to_rgba

from topoppi.atlas.uv import set_uv_layout
from topoppi.config import OptCutsConfig, VisualizationConfig
from topoppi.gui_app.app import InterfaceApp
from topoppi.gui_app.forms import parse_single_run_form
from topoppi.gui_app.workflow_mixin import WorkflowMixin
from topoppi.io.io_loader import PDBLoader
from topoppi.optimization.optcuts import OptCutsUVOptimizer
from topoppi.visualization.atlas_io import load_atlas, save_atlas
from topoppi.visualization.visualizer import InterfaceVisualizer

FIXTURE = Path(__file__).parent / "fixtures" / "tiny_complex.pdb"


def make_visualizer():
    loader = PDBLoader(FIXTURE)
    coords_a, atoms_a = loader.get_chain_data("A")
    coords_b, atoms_b = loader.get_chain_data("B")
    return InterfaceVisualizer(
        chain_A_atoms=atoms_a,
        chain_A_coords=coords_a,
        chain_B_coords=coords_b,
        chain_B_atoms=atoms_b,
        chain_a_id="A",
        chain_b_id="B",
        structure_label="tiny_complex",
        config=VisualizationConfig(map_style="footprints", min_points=1),
    )


def make_patch(indices=(0, 1, 4, 5)):
    patch = trimesh.Trimesh(
        vertices=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        faces=[[0, 1, 2], [0, 2, 3]],
        process=False,
    )
    patch.metadata["source_atom_indices"] = np.array(indices)
    set_uv_layout(patch, patch.vertices[:, :2])
    set_uv_layout(patch, patch.vertices[:, :2], key="uv_optcuts")
    set_uv_layout(patch, patch.vertices[:, :2], key="uv_global")
    return patch


class FootprintFormTests(unittest.TestCase):
    def test_form_keeps_geometric_weight_source_and_surface_cutoff_separate(self):
        form = parse_single_run_form(
            {
                "path": str(FIXTURE),
                "chain_a": "A",
                "chain_b": "B",
                "res": "2",
                "sigma": "1",
                "interaction_source": "geometric",
                "map_style": "footprints",
                "cutoff": "4",
                "contact_distance_angstrom": "6",
            }
        )
        config = form.to_config()
        self.assertEqual(config.interaction_source, "geometric")
        self.assertEqual(config.contact_distance_angstrom, 6.0)
        self.assertEqual(config.topology.distance_cutoff, 4.0)
        self.assertEqual(config.visualization.map_style, "footprints")

    def test_gui_geometric_run_bypasses_prolif_before_surface_generation(self):
        class Harness(WorkflowMixin):
            def __init__(self):
                self._cancel_event = threading.Event()
                self._successful_single_run = None
                self.messages = []
                self.errors = []

            def log(self, message):
                self.messages.append(message)

            def show_error(self, message):
                self.errors.append(message)

            def set_stage_progress(self, *args):
                pass

            @staticmethod
            def post_to_ui(callback, *args):
                callback(*args)

        form = parse_single_run_form(
            {
                "path": str(FIXTURE),
                "chain_a": "A",
                "chain_b": "B",
                "res": "2",
                "sigma": "1",
                "interaction_source": "geometric",
                "map_style": "footprints",
                "auto_save": False,
            }
        )
        harness = Harness()
        with (
            mock.patch("topoppi.gui_app.workflow_mixin.OptCutsUVOptimizer"),
            mock.patch("topoppi.gui_app.workflow_mixin.generate_prolif_interactions") as prolif,
            mock.patch("topoppi.gui_app.workflow_mixin.prepare_interface") as surface,
        ):
            surface.side_effect = RuntimeError("Reached surface generation")
            harness.run_pipeline(form.to_params(), form.to_config())
        prolif.assert_not_called()
        surface.assert_called_once()
        self.assertTrue(any("Geometric contacts:" in message for message in harness.messages))
        self.assertEqual(harness.errors, ["Reached surface generation"])


@unittest.skipUnless(os.environ.get("DISPLAY") or sys.platform in {"win32", "darwin"}, "Tk needs a graphical display")
class FootprintDesktopTests(unittest.TestCase):
    def setUp(self):
        error_dialog = mock.patch("topoppi.gui_app.ui_mixin.messagebox.showerror")
        self.errors = error_dialog.start()
        self.addCleanup(error_dialog.stop)
        self.root = tk.Tk()
        self.root.withdraw()
        self.app = InterfaceApp(self.root)
        self.app.combo_map_style.set("Residue footprints")
        self.app.combo_residue_scope.set("Full patch context")
        self.app.var_avoid_overlap.set(False)
        self.app.var_auto_save.set(False)
        self.app._remember_recent_file = lambda _path: None
        self.app._remember_recent_output_dir = lambda _path: None

    def tearDown(self):
        if not self.app._closed:
            for callback in self.root.tk.call("after", "info"):
                self.root.after_cancel(callback)
            self.app.close()

    def wait_until(self, condition, timeout=5):
        deadline = time.monotonic() + timeout
        while not condition() and time.monotonic() < deadline:
            self.root.update()
            time.sleep(0.01)
        self.assertTrue(condition(), "GUI task did not complete within the test deadline")

    def render(self, visualizer=None, patches=None, style=None):
        visualizer = visualizer or make_visualizer()
        patches = patches or [make_patch()]
        success = self.app.update_plot(
            visualizer,
            patches,
            style or self.app.get_style_config(),
            complete_task=True,
            run_params={
                "path": str(FIXTURE),
                "chain_a": "A",
                "chain_b": "B",
                "min_points": 1,
                "auto_save": False,
                "cutoff": 4.0,
                "res": 2.0,
                "sigma": 1.0,
            },
            run_manifest={"config": {}, "prolif_source": "none", "run_id": "gui-test"},
        )
        self.assertTrue(success, self.errors.call_args_list)
        return visualizer

    def test_invalid_highlight_preserves_the_previous_editable_map(self):
        visualizer = self.render()
        original_figure = self.app.current_fig
        original_artists = visualizer.artist_map
        original_style = dict(visualizer.last_style)
        self.app.var_highlight_residues.set("A:99999")
        self.app.redraw_plot()
        self.assertIs(self.app.current_fig, original_figure)
        self.assertIs(visualizer.artist_map, original_artists)
        self.assertEqual(visualizer.last_style, original_style)
        objects = next(iter(original_artists.values()))
        with mock.patch("topoppi.gui_app.plot_mixin.colorchooser.askcolor", return_value=(None, "#336699")):
            self.app.on_pick(SimpleNamespace(artist=objects["collection"]))
        np.testing.assert_allclose(objects["collection"].get_facecolor()[0], to_rgba("#336699"))

    def test_failed_auto_save_warns_and_keeps_the_map_available_for_manual_save(self):
        self.render()
        result = self.app._successful_single_run
        with tempfile.TemporaryDirectory() as folder:
            params = dict(result["params"], auto_save=True, output_dir=folder)
            with (
                mock.patch("topoppi.visualization.export.save_figure", side_effect=OSError("No space left on device")),
                mock.patch("topoppi.gui_app.plot_mixin.messagebox.showwarning") as warning,
            ):
                self.app.accept_pipeline_result(result["viz"], result["patches"], result["manifest"], params)
            warning.assert_called_once()
            self.assertIn("Save Figure", warning.call_args.args[1])
            self.assertIn("Save Atlas", warning.call_args.args[1])
            self.assertIn("No space left on device", warning.call_args.args[1])
            self.assertIsNotNone(self.app.current_fig)
            self.assertFalse(self.app._busy)
            self.assertNotIn("disabled", self.app.btn_save.state())
            destination = str(Path(folder) / "recovered.png")
            with (
                mock.patch("topoppi.gui_app.ui_mixin.filedialog.asksaveasfilename", return_value=destination),
                mock.patch("topoppi.gui_app.ui_mixin.messagebox.showinfo"),
            ):
                self.app.btn_save.invoke()
            self.assertTrue(Path(destination).is_file())
            self.assertTrue(Path(folder, "recovered.topoppi.json").is_file())

    def test_click_color_and_label_drag_survive_redraw_and_atlas_reopen(self):
        visualizer = self.render()
        uid, objects = next(iter(visualizer.artist_map.items()))
        key = objects["residue_key"]
        with mock.patch("topoppi.gui_app.plot_mixin.colorchooser.askcolor", return_value=(None, "#336699")):
            self.app.on_pick(SimpleNamespace(artist=objects["collection"]))
        self.assertEqual(self.app.residue_color_overrides[key], "#336699")
        anchor = objects["anchor"]
        self.app._drag_state = {"gid": uid}
        self.app.on_mouse_move(SimpleNamespace(inaxes=True, xdata=anchor[0] + 0.2, ydata=anchor[1] + 0.3))
        self.app.on_mouse_release(SimpleNamespace())
        self.app.redraw_plot()
        np.testing.assert_allclose(visualizer.artist_map[uid]["collection"].get_facecolor()[0], to_rgba("#336699"))
        np.testing.assert_allclose(self.app.label_offsets[uid], [0.2, 0.3])
        self.app.annotation_values = {"A:GLY:1": 0.4, "A:ALA:2": None}
        self.app.redraw_plot()
        with tempfile.TemporaryDirectory() as folder:
            path = str(Path(folder) / "map.npz")
            with mock.patch("topoppi.gui_app.ui_mixin.filedialog.asksaveasfilename", return_value=path):
                self.app.save_atlas()
            self.app.residue_color_overrides.clear()
            self.app.label_offsets.clear()
            self.app.annotation_values = None
            with (
                mock.patch("topoppi.gui_app.ui_mixin.filedialog.askopenfilename", return_value=path),
                mock.patch("topoppi.gui_app.ui_mixin.messagebox.showerror") as errors,
                mock.patch("topoppi.gui_app.workflow_mixin.OptCutsUVOptimizer") as optimizer,
                mock.patch("topoppi.gui_app.workflow_mixin.PDBLoader") as loader,
            ):
                self.app.open_atlas()
                self.app.redraw_plot()
            errors.assert_not_called()
            optimizer.assert_not_called()
            loader.assert_not_called()
        self.assertEqual(self.app.residue_color_overrides[key], "#336699")
        np.testing.assert_allclose(self.app.label_offsets[uid], [0.2, 0.3])
        self.assertEqual(self.app.annotation_values, {"A:GLY:1": 0.4, "A:ALA:2": None})
        reopened = self.app._successful_single_run["viz"]
        self.app.clear_annotations()
        np.testing.assert_allclose(reopened.artist_map[uid]["collection"].get_facecolor()[0], to_rgba("#336699"))
        self.assertEqual(self.app.var_annotation_file.get(), "")

    def test_numeric_coloring_keeps_the_colorbar_consistent_when_a_region_is_clicked(self):
        self.app.annotation_values = {"A:GLY:1": 0.4}
        visualizer = self.render()
        objects = next(iter(visualizer.artist_map.values()))
        with mock.patch("topoppi.gui_app.plot_mixin.colorchooser.askcolor") as chooser:
            self.app.on_pick(SimpleNamespace(artist=objects["collection"]))
        chooser.assert_not_called()
        self.assertFalse(self.app.residue_color_overrides)

    def test_new_csv_replaces_embedded_values_and_clear_returns_to_region_colors(self):
        self.app.annotation_values = {"A:GLY:1": 0.4}
        self.render()
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "new.csv"
            path.write_text("residue,value\nA:ALA:2,-0.7\n", encoding="utf-8")
            with mock.patch("topoppi.gui_app.ui_mixin.filedialog.askopenfilename", return_value=str(path)):
                self.app.browse_annotations()
        visualizer = self.app._successful_single_run["viz"]
        self.assertEqual(visualizer.last_style["annotation_values"], {"A:ALA:2": -0.7})
        self.assertEqual(visualizer.last_report["missing_value_residue_count"], 1)
        self.app.clear_annotations()
        self.assertIsNone(visualizer.last_style["annotation_values"])
        self.assertEqual(len(self.app.current_fig.axes), 1)

    def test_switching_to_footprints_restores_threshold_hidden_patch_without_optimization(self):
        visualizer = make_visualizer()
        selected, hidden = make_patch((0, 1, 2, 3)), make_patch((4, 5, 4, 5))

        def count(patch):
            return 1 if patch is selected else 0

        with (
            mock.patch.object(visualizer, "count_patch_interaction_residues", side_effect=count),
            mock.patch("topoppi.gui_app.workflow_mixin.OptCutsUVOptimizer") as optimizer,
        ):
            self.render(visualizer, [selected, hidden])
            self.app.var_color_type.set(False)
            self.app.combo_map_style.set("Residue markers")
            self.app._map_style_changed()
            self.assertEqual(len(self.app._successful_single_run["patches"]), 1)
            self.assertEqual(len(self.app._successful_single_run["all_patches"]), 2)
            self.app.combo_map_style.set("Residue footprints")
            self.app._map_style_changed()
        optimizer.assert_not_called()
        self.assertEqual(len(self.app._successful_single_run["patches"]), 2)
        self.assertEqual(visualizer.last_report["displayed_residue_count"], 2)

    def test_completed_geometry_survives_invalid_display_and_recovers_without_optimization(self):
        old_viz = self.render()
        self.app.var_value_min.set("unfinished")
        new_viz = make_visualizer()
        previous = self.app._successful_single_run
        self.app._begin_task("Rendering completed calculation", "determinate")
        self.app.accept_pipeline_result(new_viz, [make_patch()], previous["manifest"], previous["params"])
        self.assertFalse(self.app._busy)
        self.assertIs(self.app._successful_single_run["viz"], old_viz)
        self.assertIs(self.app._pending_single_run["viz"], new_viz)
        self.assertEqual(str(self.app.btn_save_atlas.cget("state")), "normal")
        self.app.var_value_min.set("")
        with mock.patch("topoppi.gui_app.workflow_mixin.OptCutsUVOptimizer") as optimizer:
            self.app.redraw_plot()
        optimizer.assert_not_called()
        self.assertIs(self.app._successful_single_run["viz"], new_viz)
        self.assertIsNone(self.app._pending_single_run)

    def test_invalid_appearance_has_visible_feedback_and_preserves_actual_view(self):
        self.render()
        previous = self.app.current_fig
        self.app.var_value_min.set("unfinished")
        self.app.redraw_plot()
        self.assertIs(self.app.current_fig, previous)
        self.assertIn("color limits", self.app.display_feedback_var.get())
        self.assertIn("Apply Style", self.app.display_feedback_var.get())
        self.assertIn("2D atlas", self.app.map_context_var.get())
        self.app.var_value_min.set("")
        self.app.redraw_plot()
        self.assertEqual(self.app.display_feedback_var.get(), "")
        self.assertIsNot(self.app.current_fig, previous)

    def test_failed_next_run_preserves_previous_colors_and_label_edits(self):
        self.app.var_highlight_residues.set("A:GLY:1")
        viz = self.render()
        uid, objects = next(iter(viz.artist_map.items()))
        with mock.patch("topoppi.gui_app.plot_mixin.colorchooser.askcolor", return_value=(None, "#336699")):
            self.app.on_pick(SimpleNamespace(artist=objects["collection"]))
        self.app.label_offsets[uid] = (0.2, 0.3)
        self.app._remember_interactive_style()
        self.app.var_input_path.set(str(FIXTURE))
        self.app.var_chain_b.set("Z")
        self.app.start_analysis()
        self.assertEqual(str(self.app.combo_map_style.cget("state")), "disabled")
        self.wait_until(lambda: not self.app._busy)
        self.app.redraw_plot()
        self.assertIs(self.app._successful_single_run["viz"], viz)
        np.testing.assert_allclose(viz.artist_map[uid]["collection"].get_facecolor()[0], to_rgba("#336699"))
        np.testing.assert_allclose(self.app.label_offsets[uid], (0.2, 0.3))
        self.assertEqual(viz.last_style["highlight_residues"], ["A:GLY:1"])
        self.assertEqual(str(self.app.combo_map_style.cget("state")), "readonly")
        self.assertIn("previous successful result", self.errors.call_args.args[1])

    def test_pending_atlas_can_be_saved_reopened_and_corrected_with_edits_intact(self):
        viz = self.render()
        uid, objects = next(iter(viz.artist_map.items()))
        key = objects["residue_key"]
        original = self.app._successful_single_run
        style = self.app.get_style_config()
        style.update({"highlight_residues": ("A:999",), "label_offsets": {uid: [0.2, 0.3]},
                      "residue_color_overrides": {key: "#336699"}})
        with tempfile.TemporaryDirectory() as directory:
            invalid_atlas = str(Path(directory) / "correctable.npz")
            saved_pending = str(Path(directory) / "pending.npz")
            save_atlas(invalid_atlas, [make_patch()], viz, style_config=style,
                       run_metadata={"params": original["params"], "manifest": original["manifest"]})
            with mock.patch("topoppi.gui_app.ui_mixin.filedialog.askopenfilename", return_value=invalid_atlas):
                self.app.open_atlas()
            self.assertIsNotNone(self.app._pending_single_run)
            with mock.patch("topoppi.gui_app.ui_mixin.filedialog.asksaveasfilename", return_value=saved_pending):
                self.app.save_atlas()
            self.assertTrue(Path(saved_pending).exists())
            with mock.patch("topoppi.gui_app.ui_mixin.filedialog.askopenfilename", return_value=saved_pending):
                self.app.open_atlas()
            self.assertIsNotNone(self.app._pending_single_run)
            self.app.var_highlight_residues.set("")
            self.app.redraw_plot()
        self.assertIsNone(self.app._pending_single_run)
        restored = self.app._successful_single_run["viz"]
        np.testing.assert_allclose(restored.artist_map[uid]["collection"].get_facecolor()[0], to_rgba("#336699"))
        np.testing.assert_allclose(self.app.label_offsets[uid], (0.2, 0.3))

    def test_cli_interaction_scope_survives_gui_open_redraw_and_selection(self):
        loader = PDBLoader(FIXTURE)
        coords_a, atoms_a = loader.get_chain_data("A")
        coords_b, atoms_b = loader.get_chain_data("B")
        viz = InterfaceVisualizer(atoms_a, coords_a, coords_b, atoms_b, chain_a_id="A", chain_b_id="B",
            contact_distance_angstrom=1.6,
            config=VisualizationConfig(map_style="footprints", min_points=1, residue_scope="patch",
                                       use_geometric_interaction_fallback=True))
        style = {**self.app.get_style_config(), "highlight_residues": ("A:GLY:1", "A:ALA:2")}
        with tempfile.TemporaryDirectory() as directory:
            source, target = Path(directory) / "source.npz", Path(directory) / "interaction.npz"
            save_atlas(source, [make_patch()], viz, style_config=style)
            command = [sys.executable, "-m", "topoppi.cli", "render", str(source), "--residue-scope", "interaction",
                       "-o", str(Path(directory) / "map.png"), "--export-atlas", str(target)]
            completed = subprocess.run(command, capture_output=True, text=True, timeout=30)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with mock.patch("topoppi.gui_app.ui_mixin.filedialog.askopenfilename", return_value=str(target)):
                self.app.open_atlas()
            self.app.redraw_plot()
            loaded = self.app._successful_single_run["viz"]
            self.assertEqual(loaded.last_report["displayed_label_count"], 1)
            self.assertEqual(loaded.last_style["residue_scope"], "interaction")
            self.app.combo_residue_scope.set("Full patch context")
            self.app.combo_residue_scope.event_generate("<<ComboboxSelected>>")
            self.assertEqual(loaded.last_report["displayed_label_count"], 2)
            self.app.combo_residue_scope.set("Interaction residues")
            self.app.combo_residue_scope.event_generate("<<ComboboxSelected>>")
            self.assertEqual(loaded.last_report["displayed_label_count"], 1)

    def test_next_run_csv_uses_selected_chain_and_leaves_current_map_values_intact(self):
        self.app.annotation_values = {"A:GLY:1": 0.4}
        original_viz = self.render()
        self.app.var_input_path.set(str(FIXTURE))
        self.app.var_chain_a.set("B")
        self.app.var_chain_b.set("A")
        self.app.var_annotation_target.set("Next run")
        self.app.var_highlight_residues.set("B:GLY:1")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "next.csv"
            path.write_text("residue,value\nB:1,-0.7\n", encoding="utf-8")
            with mock.patch("topoppi.gui_app.ui_mixin.filedialog.askopenfilename", return_value=str(path)):
                self.app.browse_annotations()
            self.assertEqual(self.app.annotation_values, {"A:GLY:1": 0.4})
            self.assertEqual(original_viz.last_style["annotation_values"], {"A:GLY:1": 0.4})
            with mock.patch("topoppi.gui_app.ui_mixin.threading.Thread") as thread:
                thread.return_value.is_alive.return_value = False
                self.app.start_analysis()
            self.app._finish_task()
        self.assertEqual(self.app._run_style["annotation_values"], {"B:GLY:1": -0.7})
        self.assertEqual(self.app._run_style["highlight_residues"], ("B:GLY:1",))
        self.assertIn("Displayed map: tiny_complex.pdb | A to B", self.app.map_context_var.get())
        self.assertIn("Next run: tiny_complex.pdb | B to A", self.app.map_context_var.get())
        self.app.clear_annotations()
        self.assertIsNone(self.app._next_run_annotations)
        self.assertEqual(self.app.annotation_values, {"A:GLY:1": 0.4})

    def test_next_run_units_and_scale_do_not_change_current_map_or_saved_atlas(self):
        self.app.annotation_values = {"A:GLY:1": 0.4}
        self.render()
        self.app.var_annotation_label.set("Current assay log2 enrichment")
        self.app.var_value_min.set("-1")
        self.app.var_value_max.set("1")
        self.app._annotation_style_changed()
        self.app.var_input_path.set(str(FIXTURE))
        self.app.var_chain_a.set("B")
        self.app.var_chain_b.set("A")
        self.app.var_annotation_target.set("Next run")
        self.app._annotation_target_changed()
        with tempfile.TemporaryDirectory() as directory:
            csv_path = Path(directory) / "next.csv"
            csv_path.write_text("residue,value\nB:1,-0.7\n", encoding="utf-8")
            with mock.patch("topoppi.gui_app.ui_mixin.filedialog.askopenfilename", return_value=str(csv_path)):
                self.app.browse_annotations()
            self.app.var_annotation_label.set("Next assay percent change")
            self.app.var_value_min.set("-4")
            self.app.var_value_max.set("4")
            self.app._annotation_style_changed()
            self.app.redraw_plot()
            current = self.app._successful_single_run["style"]
            self.assertEqual(current["annotation_label"], "Current assay log2 enrichment")
            self.assertEqual((current["value_min"], current["value_max"]), (-1, 1))
            self.assertEqual(current["annotation_values"], {"A:GLY:1": 0.4})
            saved = Path(directory) / "current.atlas.npz"
            with mock.patch("topoppi.gui_app.ui_mixin.filedialog.asksaveasfilename", return_value=str(saved)):
                self.app.save_atlas()
            self.assertEqual(load_atlas(saved).style["annotation_label"], "Current assay log2 enrichment")
            self.app.var_annotation_target.set("Current map")
            self.app._annotation_target_changed()
            self.assertEqual(self.app.var_annotation_label.get(), "Current assay log2 enrichment")
            with mock.patch("topoppi.gui_app.ui_mixin.threading.Thread") as thread:
                thread.return_value.is_alive.return_value = False
                self.app.start_analysis()
            self.app._finish_task()
            self.assertEqual(self.app._run_style["annotation_label"], "Next assay percent change")
            self.assertEqual((self.app._run_style["value_min"], self.app._run_style["value_max"]), (-4, 4))
            self.assertEqual(self.app._run_style["annotation_values"], {"B:GLY:1": -0.7})

    def test_successful_run_stores_effective_marker_threshold_for_atlas_reuse(self):
        self.render()
        previous = self.app._successful_single_run
        params = {**previous["params"], "min_points": 7}
        self.app._run_style = {**self.app.get_style_config(), "min_points": 99}
        self.app.accept_pipeline_result(make_visualizer(), [make_patch()], previous["manifest"], params)
        self.assertEqual(self.app._successful_single_run["style"]["min_points"], 7)
        with tempfile.TemporaryDirectory() as directory:
            path = str(Path(directory) / "threshold.npz")
            with mock.patch("topoppi.gui_app.ui_mixin.filedialog.asksaveasfilename", return_value=path):
                self.app.save_atlas()
            self.assertEqual(load_atlas(path).style["min_points"], 7)

    def test_close_waits_for_optimizer_process_cleanup(self):
        optimizer = OptCutsUVOptimizer(replace(OptCutsConfig(), residue_fragmentation_weight=0),
                                      cancel_event=self.app._cancel_event)
        launched = threading.Event()
        real_popen, processes = subprocess.Popen, []

        def launch_solver(_command, **kwargs):
            process = real_popen([sys.executable, "-c", "import time; time.sleep(60)"], **kwargs)
            processes.append(process)
            launched.set()
            return process

        def optimize():
            try:
                optimizer._run_optcuts_for_patch(make_patch(), None, 0, 60)
            except RuntimeError:
                pass

        self.app._begin_task("Testing optimizer cancellation", "determinate")
        try:
            with mock.patch.object(optimizer, "_resolved_binary", return_value=(sys.executable, "test")), \
                 mock.patch("topoppi.optimization.optcuts.joint_optimizer.subprocess.Popen", side_effect=launch_solver):
                worker = threading.Thread(target=optimize, daemon=True)
                self.app._worker_thread = worker
                worker.start()
                self.assertTrue(launched.wait(5))
                self.app.close()
                self.wait_until(lambda: self.app._closed)
                self.assertFalse(worker.is_alive())
                self.assertIsNotNone(processes[0].poll())
        finally:
            for process in processes:
                if process.poll() is None:
                    process.kill()
                process.wait(timeout=5)


if __name__ == "__main__":
    unittest.main()
