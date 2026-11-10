"""Desktop interaction scenarios for the main visualization workflow."""

import os
import sys
import tkinter as tk
from pathlib import Path
from unittest.mock import patch

import pytest

from topoppi.gui_app.app import InterfaceApp

pytestmark = pytest.mark.skipif(
    sys.platform not in {"win32", "darwin"} and not os.environ.get("DISPLAY"),
    reason="A graphical display is required",
)


@pytest.fixture
def desktop(tmp_path, monkeypatch):
    from topoppi import plugins
    monkeypatch.setattr(plugins, "_default_registry", plugins.PluginRegistry(entry_points=()))
    with patch.object(InterfaceApp, "_recent_store_path", return_value=tmp_path / "recent.json"):
        root = tk.Tk()
        root.tk.call("tk", "scaling", 1.2)
        root.geometry("1000x720")
        app = InterfaceApp(root)
        root.geometry("1000x720")
        root.update()
        try:
            yield app
        finally:
            app.close()
            try:
                root.destroy()
            except tk.TclError:
                pass


def test_first_use_prioritizes_visualization_and_explains_next_action(desktop):
    app = desktop
    assert app.btn_run.winfo_width() > 100
    assert not app._gui_extensions
    assert not app.btn_cancel.winfo_ismapped()
    assert app.basic_tab.winfo_ismapped()
    assert not app.appearance_tab.winfo_ismapped()
    assert "Choose a structure" in app.run_hint_var.get()
    assert "disabled" in app.btn_run.state()
    plugin_buttons = [child.cget("text") for child in app.plugin_actions_frame.winfo_children()
                      if child.winfo_class() == "TButton"]
    assert plugin_buttons == ["Binding analysis...", "Plugins..."]
    assert not app.combo_recent_file.winfo_ismapped()


def test_saved_atlas_reviewer_can_start_directly_from_empty_workspace(desktop):
    button = next(button for button in desktop.empty_state_actions if button.cget("text") == "Open Atlas...")
    with patch("topoppi.gui_app.ui_mixin.filedialog.askopenfilename", return_value="") as choose:
        button.invoke()
    choose.assert_called_once()
    assert "editable atlas" in choose.call_args.kwargs["title"]
    assert "editable map" in desktop.figure_help.cget("text")


def test_unreadable_structure_has_inline_recovery_guidance(desktop, tmp_path):
    path = tmp_path / "invalid.cif"
    path.write_text("invalid mmCIF content\n")
    desktop.var_input_path.set(str(path))
    desktop._commit_structure_path()
    assert "readable PDB or mmCIF" in desktop.run_hint_var.get()
    assert "readable PDB or mmCIF" in desktop.chain_preview_var.get()
    assert "disabled" in desktop.btn_run.state()
    desktop.var_input_path.set(str(Path(__file__).parent / "fixtures" / "synthetic_complex.pdb"))
    desktop._commit_structure_path()
    assert "A: 2 residues" in desktop.chain_preview_var.get()
    assert "disabled" not in desktop.btn_run.state()


def test_benchmark_missing_manifest_stays_actionable_and_recovers(desktop, tmp_path):
    desktop.var_run_mode.set("benchmark")
    desktop.var_input_path.set(str(tmp_path))
    desktop.var_chain_selection_mode.set("manifest")
    desktop._form_changed()
    assert "existing benchmark manifest CSV" in desktop.run_hint_var.get()
    assert "disabled" in desktop.btn_run.state()
    manifest = tmp_path / "manifest.csv"
    manifest.write_text("structure_path,chain_a,chain_b\nexample.pdb,A,B\n")
    desktop.entry_manifest_path.insert(0, str(manifest))
    desktop._form_changed()
    assert desktop.run_hint_var.get() == "Ready to run the selected benchmark."
    assert "disabled" not in desktop.btn_run.state()


def test_example_prepares_approved_input_and_writable_output(desktop, tmp_path):
    with patch("topoppi.gui_app.ui_mixin.filedialog.askdirectory", return_value=str(tmp_path)):
        desktop.btn_example.invoke()
    assert Path(desktop.var_input_path.get()).is_file()
    assert desktop.var_chain_a.get() == "B"
    assert desktop.var_chain_b.get() == "C"
    assert Path(desktop.var_output_dir.get()).is_relative_to(tmp_path)
    assert desktop.chain_residue_counts["C"] == 129
    assert "disabled" not in desktop.btn_run.state()
    assert not desktop._busy


def test_pasted_path_detects_chains_and_invalid_selection_has_inline_guidance(desktop):
    path = Path(__file__).parent / "fixtures" / "synthetic_complex.pdb"
    desktop.var_input_path.set(str(path))
    desktop.entry_file.focus_force()
    desktop.root.update()
    desktop.entry_file.event_generate("<Return>")
    desktop.root.update()
    assert "B" in desktop.chain_residue_counts
    desktop.var_chain_b.set(desktop.var_chain_a.get())
    desktop._form_changed()
    assert "different chains" in desktop.run_hint_var.get()
    assert "disabled" in desktop.btn_run.state()
    desktop.var_chain_b.set("Z")
    desktop._form_changed()
    assert "present in this structure" in desktop.run_hint_var.get()
    assert "disabled" in desktop.btn_run.state()
    desktop.var_chain_b.set("B")
    desktop._form_changed()
    assert "disabled" not in desktop.btn_run.state()


def test_busy_cancel_and_retry_restore_the_workflow(desktop):
    desktop.var_input_path.set(str(Path(__file__).parent / "fixtures" / "synthetic_complex.pdb"))
    desktop._commit_structure_path()
    desktop._begin_task("Creating interface map...", progress_mode="determinate")
    desktop.root.update()
    assert desktop.btn_cancel.winfo_ismapped()
    assert "disabled" in desktop.btn_run.state()
    assert all("disabled" in button.state() for button in desktop.empty_state_actions)
    assert "disabled" in desktop.combo_view.state()
    desktop.btn_cancel.invoke()
    assert desktop._cancel_event.is_set()
    assert "disabled" in desktop.btn_cancel.state()
    desktop.finish_cancelled("Run cancelled.")
    desktop.root.update()
    assert not desktop.btn_cancel.winfo_ismapped()
    assert "disabled" not in desktop.btn_run.state()
    assert not desktop._extensions_need_restart()
    assert str(desktop.combo_view.cget("state")) == "readonly"


def test_settings_navigation_keeps_selected_benchmark_mode(desktop, tmp_path):
    desktop.var_input_path.set(str(tmp_path))
    desktop.var_run_mode.set("benchmark")
    desktop._sync_mode_controls()
    for page in ("appearance", "basic", "advanced"):
        desktop.var_settings_page.set(page)
        desktop._sync_settings_page()
        assert desktop.var_run_mode.get() == "benchmark"
        assert desktop.btn_run.cget("text") == "Run Benchmark"


def test_entering_sidebar_keeps_keyboard_focus_in_input(desktop):
    desktop.entry_file.focus_force()
    desktop.root.update()
    desktop._bind_sidebar_wheel(None)
    assert desktop.root.focus_get() is desktop.entry_file


def test_benchmark_setup_uses_matching_browse_action_and_guidance(desktop, tmp_path):
    desktop.var_input_path.set(str(tmp_path))
    desktop.var_run_mode.set("benchmark")
    desktop._sync_mode_controls()
    desktop.var_settings_page.set("basic")
    desktop._sync_settings_page()
    assert "Run Benchmark" in desktop.empty_state_hint.get()
    assert desktop.btn_browse_input.cget("text") == "Folder..."
    with patch.object(desktop, "browse_folder") as folder, patch.object(desktop, "browse_file") as file:
        desktop.btn_browse_input.invoke()
    folder.assert_called_once()
    file.assert_not_called()


def test_single_run_shows_activity_without_inventing_elapsed_percentage(desktop):
    desktop._begin_task("Creating interface map...", progress_mode="determinate")
    desktop._set_stage_progress_ui("OptCuts", 64, "Optimizing UVs")
    assert str(desktop.progress.cget("mode")) == "indeterminate"
    assert desktop.stage_status_var.get() == "OptCuts: Optimizing UVs"
    desktop._append_log_message("Detailed optimizer report")
    assert desktop.status_var.get() == "OptCuts: Optimizing UVs"
    desktop._finish_task()


def test_activity_details_and_diagnostic_actions_expand_together(desktop):
    toggle = next(child for child in desktop.activity_frame.winfo_children()[0].winfo_children()
                  if child.winfo_class() == "TCheckbutton")
    actions = desktop.activity_details.winfo_children()[0].winfo_children()
    assert {button.cget("text") for button in actions} == {"Copy diagnostics", "Save log..."}
    assert not desktop.log_text.winfo_ismapped()
    assert all(not button.winfo_ismapped() for button in actions)
    desktop._append_log_message("Interface map created")
    toggle.invoke()
    desktop.root.update()
    assert desktop.log_text.winfo_ismapped()
    assert "Interface map created" in desktop.log_text.get("1.0", tk.END)
    assert all(button.winfo_ismapped() for button in actions)
    toggle.invoke()
    desktop.root.update()
    assert not desktop.log_text.winfo_ismapped()
    assert all(not button.winfo_ismapped() for button in actions)


def test_plot_navigation_controls_track_buttons_and_keyboard_modes(desktop):
    from matplotlib.backend_bases import KeyEvent, key_press_handler
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.figure import Figure

    from topoppi.gui_app.plot_mixin import _MapToolbar

    fig = Figure()
    fig.add_subplot()
    desktop.current_fig = fig
    desktop.current_canvas = FigureCanvasTkAgg(fig, master=desktop.canvas_frame)
    desktop.current_toolbar = _MapToolbar(desktop.current_canvas, desktop.canvas_frame, pack_toolbar=False,
                                          on_mode_changed=desktop._sync_plot_navigation)
    desktop._build_plot_toolbar(desktop.canvas_frame)
    pan = desktop.plot_navigation_buttons["Pan"]
    zoom = desktop.plot_navigation_buttons["Zoom"]
    pan.invoke()
    assert str(pan.cget("relief")) == "sunken"
    assert "Pan active" in desktop.plot_navigation_hint.get()
    zoom.invoke()
    assert str(pan.cget("relief")) == "flat"
    assert str(zoom.cget("relief")) == "sunken"
    key_press_handler(KeyEvent("key_press_event", desktop.current_canvas, key="o"),
                      desktop.current_canvas, desktop.current_toolbar)
    assert str(zoom.cget("relief")) == "flat"
    assert "Drag labels" in desktop.plot_navigation_hint.get()


def test_navigation_gestures_leave_residue_edits_unchanged(desktop):
    from types import SimpleNamespace

    desktop.current_toolbar = SimpleNamespace(mode="pan/zoom")
    desktop._drag_state = {"gid": "previous-label"}
    with patch("topoppi.gui_app.plot_mixin.colorchooser.askcolor") as color:
        desktop.on_mouse_press(SimpleNamespace())
        desktop.on_pick(SimpleNamespace())
        desktop.on_mouse_move(SimpleNamespace())
    assert desktop._drag_state is None
    color.assert_not_called()
    desktop.current_toolbar = None


def _descendants(widget):
    for child in widget.winfo_children():
        yield child
        yield from _descendants(child)


def test_binding_launcher_refreshes_context_and_displayed_result_availability(desktop):
    from types import SimpleNamespace
    from unittest.mock import Mock

    desktop.show_binding_analysis()
    window = desktop._binding_window
    displayed = next(child for child in _descendants(window)
                     if child.winfo_class() == "TRadiobutton" and child.cget("text") == "Displayed result")
    assert "disabled" in displayed.state()
    desktop.var_input_path.set("updated.pdb")
    desktop.var_chain_a.set("A")
    desktop.var_chain_b.set("B")
    desktop._successful_single_run = {"params": {"path": "shown.pdb", "chain_a": "C", "chain_b": "D"},
                                    "viz": SimpleNamespace(chain_a_id="C", chain_b_id="D")}
    desktop._sync_extensions()
    assert "disabled" not in displayed.state()
    labels = [child for child in _descendants(window) if child.winfo_class() == "TLabel" and child.cget("textvariable")]
    assert any("updated.pdb" in str(child.getvar(child.cget("textvariable"))) for child in labels)
    opened = Mock()
    desktop.register_binding_task("mutation", title="Mutation effects", description="", plugin="affinity", open_task=opened)
    displayed.invoke()
    next(child for child in _descendants(window)
         if child.winfo_class() == "TButton" and child.cget("text") == "Mutation effects").invoke()
    assert opened.call_args.kwargs["setup"] == {"structure": "shown.pdb", "chains_a": "C", "chains_b": "D", "prediction_task": "mutation"}
    desktop._successful_single_run = None


def test_compact_variant_composer_requires_real_changes_and_previews_joint_grouping(desktop):
    from unittest.mock import Mock

    opened = Mock()
    desktop.register_binding_task("mutation", title="Mutation effects", description="", plugin="affinity", open_task=opened)
    setup = {"structure": str(Path(__file__).parent / "fixtures" / "synthetic_complex.pdb"), "chains_a": "A", "chains_b": "B"}
    composer = desktop.create_variants_from_residues([
        {"chain": "A", "position": 1, "wt": "A", "insertion": ""},
        {"chain": "B", "position": 52, "wt": "Y", "insertion": "A"},
    ], setup)
    composer.window.geometry("620x480")
    desktop.root.update()
    assert "disabled" in composer.add_button.state()
    assert "already A" in composer.feedback.get()
    composer.table.selection_set("A:1")
    desktop.root.update()
    composer.target.set("G")
    next(child for child in _descendants(composer.window) if child.winfo_class() == "TCombobox").event_generate("<<ComboboxSelected>>")
    joint = next(child for child in _descendants(composer.window)
                 if child.winfo_class() == "TRadiobutton" and child.cget("text").startswith("One joint"))
    joint.invoke()
    desktop.root.update()
    assert len(composer.preview.get_children()) == 1
    assert "disabled" not in composer.add_button.state()
    assert composer.add_button.winfo_rooty() + composer.add_button.winfo_height() <= composer.window.winfo_rooty() + composer.window.winfo_height()
    composer.add_button.invoke()
    sent = opened.call_args.kwargs["setup"]
    assert sent["append_variants"] == [{"variant_id": "joint_variant", "mutations": "A:A1G,B:Y52AA"}]
    assert sent["structure"] == setup["structure"]
    assert not composer.window.winfo_exists()


def test_missing_source_has_an_explicit_checked_recovery_without_rewriting_provenance(desktop, tmp_path):
    from unittest.mock import Mock

    original = {"structure": str(tmp_path / "moved.pdb"), "chains_a": "A", "chains_b": "B"}
    opened = Mock()
    desktop.register_binding_task("mutation", title="Mutation effects", description="", plugin="affinity", open_task=opened)
    composer = desktop.create_variants_from_residues([{"chain": "A", "position": 1, "wt": "A"}], original)
    composer.targets["A:1"].set("G")
    assert "disabled" in composer.add_button.state()
    assert "Locate structure" in composer.feedback.get()
    assert len(composer.preview.get_children()) == 1
    wrong = Path(__file__).parent / "fixtures" / "tiny_complex.pdb"
    with patch("tkinter.filedialog.askopenfilename", return_value=str(wrong)):
        composer.locate_structure()
    assert "does not match selected residue" in composer.feedback.get()
    assert composer.setup == original
    matching = Path(__file__).parent / "fixtures" / "synthetic_complex.pdb"
    with patch("tkinter.filedialog.askopenfilename", return_value=str(matching)):
        composer.locate_structure()
    assert "disabled" not in composer.add_button.state()
    composer.add_button.invoke()
    assert opened.call_args.kwargs["setup"]["structure"] == str(matching.resolve())
    assert original["structure"] == str(tmp_path / "moved.pdb")


def test_long_joint_preview_can_be_scrolled_and_copied_without_losing_changes(desktop):
    setup = {"structure": str(Path(__file__).parent / "fixtures" / "synthetic_complex.pdb"), "chains_a": "A", "chains_b": "B"}
    composer = desktop.create_variants_from_residues([
        {"chain": "B", "position": number, "wt": "Y"} for number in range(100, 113)
    ], setup)
    composer.window.geometry("620x480")
    composer.mode.set("joint")
    composer.refresh_preview()
    desktop.root.update()
    row = composer.preview.get_children()[0]
    expected = ",".join(f"B:Y{number}A" for number in range(100, 113))
    assert composer.preview.xview()[1] < 1
    composer.preview.xview_moveto(1)
    assert composer.preview.xview()[1] == 1
    composer.preview.selection_set(row)
    composer.preview.focus_force()
    desktop.root.update()
    composer.preview.event_generate("<Control-c>")
    desktop.root.update()
    assert composer.window.clipboard_get() == f"variant_id\tmutations\njoint_variant\t{expected}"
    assert composer.copy_button.cget("text") == "Copied 1 variant"
    assert composer.add_button.winfo_rooty() + composer.add_button.winfo_height() <= composer.window.winfo_rooty() + composer.window.winfo_height()
