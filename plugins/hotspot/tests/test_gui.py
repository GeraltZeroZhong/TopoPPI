"""Actual Tk interactions with an independently generated structure and fixed model."""
import json
import time
import tkinter as tk
from pathlib import Path
from tkinter import ttk
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest
from topoppi_hotspot.gui import create_controller
from topoppi_hotspot.gui_workflow import RankingDocument
from topoppi_hotspot.plugin import plugin

import topoppi
from topoppi import plugins

FIXTURE = Path(__file__).parent / "fixtures/tiny_complex.pdb"


@pytest.fixture
def workspace(monkeypatch, tmp_path):
    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("A display is required for actual Tk interaction")
    root.withdraw()
    registry = plugins.PluginRegistry(entry_points=[])
    registry.register_plugin(plugin)
    monkeypatch.setattr(plugins, "_default_registry", registry)
    geometry = topoppi.compute_geometry(topoppi.get_interface(topoppi.load_structure(FIXTURE), "A", "B"))
    host = SimpleNamespace(root=root, analysis_menu=tk.Menu(root), plugin_actions_frame=ttk.Frame(root),
        var_input_path=tk.StringVar(master=root, value=str(FIXTURE)), var_chain_a=tk.StringVar(master=root, value="A"),
        var_chain_b=tk.StringVar(master=root, value="B"), _busy=False, _closing=False,
        get_current_geometry=lambda: geometry, update_plot=Mock(return_value=True))
    controller = create_controller(host)
    controller.initialize()
    panel = controller.open()
    panel.window.withdraw()
    yield controller, panel, tmp_path
    controller.close()
    root.destroy()


def test_prefill_missing_model_and_real_worker_results(workspace, monkeypatch):
    controller, panel, tmp_path = workspace
    assert panel.vars["structure"].get() == str(FIXTURE)
    assert panel.vars["chains_a"].get() == "A" and panel.vars["chains_b"].get() == "B"
    panel.vars["model_dir"].set(str(tmp_path / "missing-model"))
    panel.run_button.invoke()
    assert "Choose a model folder" in panel.error.get()
    assert panel.worker is None
    model = tmp_path / "model"
    model.mkdir()
    (model / "model.json").write_text(json.dumps({"kind": "linear", "columns": ["surface.area"], "seeds": [0]}))
    np.savez(model / "linear_model.npz", mean=[0.], scale=[1.], coefficient=[1.])
    panel.vars["model_dir"].set(str(model))
    panel.vars["output_dir"].set(str(tmp_path / "results"))
    panel.run_button.invoke()
    deadline = time.monotonic() + 15
    while panel._working and time.monotonic() < deadline:
        controller.host.root.update()
        time.sleep(.01)
    assert not panel._working, "The CPU ranking worker did not finish"
    assert panel.document is not None, panel.error.get()
    assert panel.document.metadata["geometry_reused"] is True
    assert panel.table.get_children() == ("A:GLY:1",)
    assert panel.show_surface()
    assert controller.host.update_plot.call_args.args[2]["view"] == "surface"
    assert controller.host._pending_single_run["geometry"] is panel.document.geometry
    target = tmp_path / "core-scores.csv"
    monkeypatch.setattr("topoppi_hotspot.gui.filedialog.asksaveasfilename", lambda **_: str(target))
    panel.export(annotations=True)
    assert target.read_text().startswith("residue,value\n")
    assert "CSV saved" in panel.status.get()


def test_worker_error_keeps_controls_ready_for_correction(workspace):
    controller, panel, tmp_path = workspace
    model = tmp_path / "invalid-model"
    model.mkdir()
    (model / "model.json").write_text('{"kind": "unsupported"}')
    panel.vars["model_dir"].set(str(model))
    panel.vars["output_dir"].set(str(tmp_path / "results"))
    panel.run_button.invoke()
    deadline = time.monotonic() + 10
    while panel._working and time.monotonic() < deadline:
        controller.host.root.update()
        time.sleep(.01)
    assert "Unsupported frozen model kind" in panel.error.get()
    assert str(panel.run_button.cget("state")) == "normal"
    assert str(panel.cancel_button.cget("state")) == "disabled"
    assert panel.document is None


def test_host_plot_images_are_collected_before_ranking_worker(workspace, monkeypatch):
    import gc
    import threading
    import weakref

    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure

    controller, panel, tmp_path = workspace
    host, root = controller.host, controller.host.root
    owner = threading.get_ident()
    events = []
    host.canvas_frame = ttk.Frame(root)

    def released(name):
        events.append((name, threading.get_ident()))

    def install_host_plot():
        host.figure = Figure(figsize=(2, 2))
        host.figure.subplots().plot([0, 1], [0, 1])
        host.current_canvas = FigureCanvasTkAgg(host.figure, master=host.canvas_frame)
        host.current_toolbar = NavigationToolbar2Tk(host.current_canvas, host.canvas_frame, pack_toolbar=False)
        host.current_canvas.draw()
        image = host.current_toolbar._buttons["Home"]._ntimage
        weakref.finalize(image, released, "toolbar_image")
        weakref.finalize(host.current_canvas._tkphoto, released, "canvas_image")
        return weakref.ref(image), weakref.ref(host.current_canvas._tkphoto)

    def clear_host_plot():
        for child in host.canvas_frame.winfo_children():
            child.destroy()
        host.figure = host.current_canvas = host.current_toolbar = None

    def computation(*args, **kwargs):
        events.append(("worker", threading.get_ident()))
        gc.collect()
        return ranking_document(tmp_path)

    monkeypatch.setattr("topoppi_hotspot.gui.run_ranking", computation)
    panel.vars["output_dir"].set(str(tmp_path / "results"))
    toolbar_ref, canvas_ref = install_host_plot()
    gc.collect()
    clear_host_plot()

    def finish():
        if panel.worker is not None and not panel.worker.is_alive():
            root.quit()
        else:
            root.after(10, finish)

    deadline = root.after(5000, root.quit)
    root.after(0, panel.start)
    root.after(10, finish)
    root.mainloop()
    root.after_cancel(deadline)
    assert panel.worker is not None and not panel.worker.is_alive()
    assert toolbar_ref() is None and canvas_ref() is None
    assert {name for name, _ in events} == {"toolbar_image", "canvas_image", "worker"}
    assert all(thread == owner for name, thread in events if name != "worker")
    assert next(i for i, event in enumerate(events) if event[0] == "worker") == 2


def ranking_document(tmp_path):
    return RankingDocument(tmp_path, {"structure": str(FIXTURE), "partner_groups": [["A"], ["B"]], "selected_partner": ["A"]}, [
        {"residue_key": "A:GLY:1", "chain": "A", "residue_number": 1, "insertion_code": "", "wt": "G",
         "rank_score": .8, "rank_min": 1, "rank_max": 1},
        {"residue_key": "A:TYR:2A", "chain": "A", "residue_number": 2, "insertion_code": "A", "wt": "Y",
         "rank_score": .2, "rank_min": 2, "rank_max": 2},
    ])


@pytest.mark.parametrize("temperature, label", [(None, "298.15 K"), (310.0, "310 K")])
def test_ddg_summary_uses_saved_reference_temperature(workspace, temperature, label):
    _controller, panel, tmp_path = workspace
    document = ranking_document(tmp_path)
    document.metadata["target"] = "ddg"
    if temperature is not None:
        document.metadata["temperature_K"] = temperature
    panel.display(document)
    assert f"Predicted Ala ddG (kcal/mol) at {label}." in panel.summary.get()
    assert panel.table.heading("score", "text") == "ddG (kcal/mol)"


def test_selected_residues_keep_source_and_author_identity(workspace):
    controller, panel, tmp_path = workspace
    controller.host.create_variants_from_residues = Mock(return_value=True)
    panel.display(ranking_document(tmp_path))
    panel.vars["structure"].set("a-different-structure.pdb")
    panel.vars["chains_a"].set("C")
    panel.vars["chains_b"].set("D")
    panel.table.selection_set(*panel.table.get_children())
    panel.selected()
    assert "2 residues selected" in panel.selection_status.get()
    assert str(panel.table.cget("selectmode")) == "extended"
    assert panel.create_variants()
    controller.host.create_variants_from_residues.assert_called_once_with([
        {"chain": "A", "position": 1, "insertion": "", "wt": "G", "residue_key": "A:GLY:1"},
        {"chain": "A", "position": 2, "insertion": "A", "wt": "Y", "residue_key": "A:TYR:2A"},
    ], setup={"structure": str(FIXTURE), "chains_a": "A", "chains_b": "B"})
    assert panel.vars["structure"].get() == "a-different-structure.pdb"


def test_task_buttons_use_current_inputs_while_variants_keep_result_source(workspace):
    from topoppi.gui_app.binding_analysis import BindingAnalysisMixin

    controller, panel, tmp_path = workspace
    panel.close()
    host = controller.host
    host._navigate_binding_task = Mock()
    host.add_binding_navigation = lambda parent, current_task, setup_provider: BindingAnalysisMixin.add_binding_navigation(
        host, parent, current_task, setup_provider)
    panel = controller.open()
    panel.display(ranking_document(tmp_path))
    next_input = {"structure": "next-complex.pdb", "chains_a": "C,D", "chains_b": "E"}
    controller.open(setup=next_input)
    navigation = panel.window.winfo_children()[0].winfo_children()[0]
    button = next(child for child in navigation.winfo_children()
                  if child.winfo_class() == "TButton" and child.cget("text") == "Binding affinity")
    button.invoke()
    assert host._navigate_binding_task.call_args.args[:2] == ("absolute", next_input)
    host.create_variants_from_residues = Mock(return_value=True)
    assert panel.create_variants()
    assert host.create_variants_from_residues.call_args.kwargs["setup"] == {
        "structure": str(FIXTURE), "chains_a": "A", "chains_b": "B"}


def test_variant_handoff_requires_selection_and_idle_workspace(workspace):
    controller, panel, tmp_path = workspace
    controller.host.create_variants_from_residues = Mock(return_value=True)
    panel.display(ranking_document(tmp_path))
    panel.table.selection_remove(*panel.table.selection())
    panel.selected()
    assert str(panel.variants_button.cget("state")) == "disabled"
    assert panel.create_variants() is False
    panel.table.selection_set(panel.table.get_children()[0])
    panel._set_working(True)
    assert str(panel.variants_button.cget("state")) == "disabled"
    assert panel.create_variants() is False
    assert "Wait for" in panel.error.get()
    with pytest.raises(ValueError, match="Wait for the Hotspot"):
        controller.open(setup={"structure": "other.pdb"})
    assert panel.vars["structure"].get() == str(FIXTURE)
    panel._set_working(False)
    controller.host._busy = True
    assert panel.create_variants() is False
    controller.host.create_variants_from_residues.assert_not_called()
    controller.host._busy = False


def test_unavailable_affinity_preserves_ranking(workspace):
    controller, panel, tmp_path = workspace
    document = ranking_document(tmp_path)
    panel.display(document)
    controller.host.create_variants_from_residues = Mock(side_effect=ValueError("Enable Affinity in Plugins to create variants."))
    assert panel.create_variants() is False
    assert panel.error.get() == "Enable Affinity in Plugins to create variants."
    assert panel.document is document
    assert panel.table.get_children() == ("A:GLY:1", "A:TYR:2A")


def test_old_table_without_identity_is_reviewable_but_not_guessed(workspace):
    controller, panel, tmp_path = workspace
    document = ranking_document(tmp_path)
    del document.rows[0]["wt"]
    controller.host.create_variants_from_residues = Mock()
    panel.display(document)
    assert panel.create_variants() is False
    assert "lacks explicit residue identities" in panel.error.get()
    assert str(panel.review_buttons[1].cget("state")) == "normal"
    controller.host.create_variants_from_residues.assert_not_called()


def test_export_retry_clears_failed_destination_feedback(workspace, monkeypatch):
    _controller, panel, tmp_path = workspace
    panel.display(ranking_document(tmp_path))
    destination = tmp_path / "missing" / "rankings.csv"
    monkeypatch.setattr("topoppi_hotspot.gui.filedialog.asksaveasfilename", lambda **_: str(destination))
    panel.export()
    assert "writable CSV destination" in panel.error.get()
    destination = tmp_path / "recovered.csv"
    panel.export()
    assert destination.is_file()
    assert "CSV saved" in panel.status.get()
    assert panel.error.get() == ""


def test_long_result_path_wraps_inside_compact_workspace(workspace):
    from topoppi.config import DEFAULT_GUI_CONFIG
    from topoppi.gui_app.app import InterfaceApp

    controller, panel, tmp_path = workspace
    controller.host.config = DEFAULT_GUI_CONFIG
    InterfaceApp._configure_style(controller.host)
    panel.display(ranking_document(tmp_path))
    panel.window.geometry("850x650")
    panel.window.deiconify()
    panel.status.set("CSV saved · /project/" + "long-experiment-name-" * 8 + "/exports/ranked-residue-values.csv")
    controller.host.root.update()
    outer = panel.window.winfo_children()[0]
    label = next(child for child in outer.winfo_children()
                 if child.winfo_class() == "TLabel" and str(child.cget("textvariable")) == str(panel.status))
    assert label.winfo_reqwidth() <= label.winfo_width()
    assert panel.table.winfo_height() >= 45


def test_surface_review_and_reopen_keep_all_selected_residues(workspace):
    controller, panel, tmp_path = workspace
    document = ranking_document(tmp_path)
    document.geometry = controller.host.get_current_geometry()
    for row in document.rows:
        row["footprint_area_A2"] = 1.0
    panel.display(document)
    selection = panel.table.get_children()
    panel.table.selection_set(selection)
    panel.selected()
    assert panel.show_surface()
    assert controller.host.update_plot.call_args.args[2]["highlight_residues"] == selection
    manifest = controller.host.update_plot.call_args.kwargs["run_manifest"]
    assert manifest["hotspot"]["selected_residues"] == list(selection)
    assert manifest["hotspot"]["selected"] == selection[0]
    controller.document_loaded(manifest)
    assert panel.table.selection() == selection
    panel.close()
    reopened = controller.open()
    assert reopened.table.selection() == selection
    assert "2 residues selected" in reopened.selection_status.get()
    # Existing saved figures recorded a single selected key.
    del manifest["hotspot"]["selected_residues"]
    controller.document_loaded(manifest)
    assert reopened.table.selection() == selection[:1]


def test_closing_and_reopening_retains_unfinished_ranking_settings(workspace):
    controller, panel, tmp_path = workspace
    panel.display(ranking_document(tmp_path))
    panel.vars["model_dir"].set(str(tmp_path / "chosen-model"))
    panel.vars["output_dir"].set(str(tmp_path / "next-experiment"))
    panel.vars["embedding_registry"].set(str(tmp_path / "prepared-embeddings"))
    panel.partner.set("Partner 2")
    before = {name: value.get() for name, value in panel.vars.items()}
    panel.close()
    reopened = controller.open(setup={"structure": before["structure"], "chains_a": "A", "chains_b": "B"})
    assert {name: value.get() for name, value in reopened.vars.items()} == before
    assert reopened.partner.get() == "Partner 2"
    assert reopened.document.metadata["selected_partner"] == ["A"]
    assert "Ranking setup restored" in reopened.status.get()
    assert reopened.worker is None


def test_saved_snapshot_recovers_partner_groups_without_replacing_interface_metadata(workspace):
    controller, _panel, tmp_path = workspace
    document = ranking_document(tmp_path)
    manifest = {"hotspot": {"directory": str(tmp_path), "metadata": document.metadata, "rows": document.rows}}
    run = controller.host._successful_single_run = {"manifest": dict(manifest)}
    controller.document_loaded(manifest)
    assert run["manifest"]["interface"] == {
        "source": {"path": str(FIXTURE)}, "groups": [["A"], ["B"]]}
    original = {"source": {"path": "original.pdb"}, "groups": [["C"], ["D"]], "analysis_id": "saved"}
    run["manifest"]["interface"] = original
    controller.document_loaded(manifest)
    assert run["manifest"]["interface"] is original


def test_shared_navigation_replaces_separate_plugin_button(workspace):
    controller, panel, _ = workspace
    panel.close()
    controller.host.register_binding_task = Mock()
    controller.host.add_binding_navigation = Mock()
    shared = create_controller(controller.host)
    shared.initialize()
    controller.host.register_binding_task.assert_called_once_with(
        "hotspot", title="Rank hotspots", description="Rank interface residues within one binding partner.",
        plugin="hotspot", open_task=shared.open)
    assert shared.button is None and shared.menu_index is None
    shared.sync_controls()
    new_panel = shared.open()
    controller.host.add_binding_navigation.assert_called_once()
    assert controller.host.add_binding_navigation.call_args.kwargs["current_task"] == "hotspot"
    new_panel.close()


@pytest.mark.parametrize("size", ["1040x760", "850x650"])
@pytest.mark.parametrize("sequence_open", [False, True])
def test_compact_workspace_keeps_result_actions_inside_window(workspace, size, sequence_open):
    from topoppi.config import DEFAULT_GUI_CONFIG
    from topoppi.gui_app.app import InterfaceApp
    from topoppi.gui_app.binding_analysis import BindingAnalysisMixin

    controller, panel, tmp_path = workspace
    panel.close()
    host = controller.host
    host.config = DEFAULT_GUI_CONFIG
    InterfaceApp._configure_style(host)
    host.add_binding_navigation = lambda parent, current_task, setup_provider: BindingAnalysisMixin.add_binding_navigation(
        host, parent, current_task, setup_provider)
    panel = controller.open()
    panel.display(ranking_document(tmp_path))
    panel.window.geometry(size)
    panel.window.deiconify()
    if sequence_open:
        panel.sequence_button.invoke()
    host.root.update()
    for widget in [panel.variants_button, *panel.review_buttons]:
        assert widget.winfo_viewable()
        assert widget.winfo_rootx() >= panel.window.winfo_rootx()
        assert widget.winfo_rooty() >= panel.window.winfo_rooty()
        assert widget.winfo_rootx() + widget.winfo_width() <= panel.window.winfo_rootx() + panel.window.winfo_width()
        assert widget.winfo_rooty() + widget.winfo_height() <= panel.window.winfo_rooty() + panel.window.winfo_height()
    assert panel.table.winfo_height() >= 45, "Keep the ranking header and at least one selectable row visible"


def test_default_model_runs_without_model_selection(workspace):
    controller, panel, tmp_path = workspace
    assert panel.vars["model_dir"].get() == ""
    panel.vars["output_dir"].set(str(tmp_path / "official-results"))
    panel.run_button.invoke()
    deadline = time.monotonic() + 20
    while panel._working and time.monotonic() < deadline:
        controller.host.root.update()
        time.sleep(.01)
    assert not panel._working
    assert panel.document is not None, panel.error.get()
    assert panel.document.metadata["model_id"] == "topoppi-hotspot-mutation-surface"
    assert panel.table.heading("score", "text") == "ddG (kcal/mol)"
    assert panel.table.get_children() == ("A:GLY:1",)
    assert panel.show_surface()
