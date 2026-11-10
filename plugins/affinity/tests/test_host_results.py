"""Host/plugin result journeys: retained canvases, partners and portable numbers."""

import builtins
import json
import os
import shutil
import sys
import tkinter as tk
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import trimesh
from topoppi.atlas.uv import set_uv_layout
from topoppi.config import SurfaceConfig, TopologyConfig
from topoppi.gui_app.app import InterfaceApp
from topoppi.interface import PartnerInterface, PreparedInterface
from topoppi.interface_io import save_interface
from topoppi.io.io_loader import PDBLoader
from topoppi.visualization.atlas_io import load_atlas, save_atlas

from topoppi_affinity.gui import AffinityPanel, variant_lines
from topoppi_affinity.workflow import RESIDUE_COLUMNS, load_analysis

pytestmark = pytest.mark.skipif(
    sys.platform not in {"win32", "darwin"} and not os.environ.get("DISPLAY"),
    reason="A graphical display is required",
)
FIXTURE = Path(__file__).parent / "fixtures" / "tiny_complex.pdb"


@pytest.fixture
def desktop(tmp_path, monkeypatch):
    monkeypatch.setenv("TOPOPPI_AFFINITY_CACHE", str(tmp_path / "resources"))
    monkeypatch.setattr(InterfaceApp, "_recent_store_path", lambda _: tmp_path / "recent.json")
    monkeypatch.setattr(AffinityPanel, "_missing_dependencies", staticmethod(lambda task: []))
    errors = []
    monkeypatch.setattr("tkinter.messagebox.showerror", lambda *args, **kwargs: errors.append(args))
    root = tk.Tk()
    root.withdraw()
    app = InterfaceApp(root)
    root.update()
    try:
        yield app, errors
    finally:
        app.close()
        if not app._closed:
            root.destroy()


@pytest.fixture
def analysis_factory(tmp_path):
    def create(name="analysis", *, failure=None):
        output = tmp_path / name
        output.mkdir()
        sites_a = [{"chain": "A", "position": 1, "insertion": "", "wt": "G", "mt": "W",
                    "residue_key": "A:GLY:1", "side": 0, "coverage": "footprint"}]
        sites_b = [{"chain": "B", "position": 1, "insertion": "", "wt": "G", "mt": "A",
                    "residue_key": "B:GLY:1", "side": 1, "coverage": "footprint"}]
        rows = pd.DataFrame([
            {"variant_id": "left", "mutations": "A:G1W", "predicted": 1.25, "sites": json.dumps(sites_a)},
            {"variant_id": "joint", "mutations": "A:G1W,B:G1A", "predicted": .75, "sites": json.dumps(sites_a + sites_b)},
            {"variant_id": "right", "mutations": "B:G1A", "predicted": -.6, "sites": json.dumps(sites_b)},
            {"variant_id": "invalid", "mutations": "A:G99W", "predicted": "", "sites": "[]"},
        ])
        rows["sample_id"] = rows.variant_id
        rows["status"], rows["error"] = "ok", ""
        rows.loc[rows.variant_id == "invalid", ["status", "error"]] = ["error", "Mutation residue absent: A:99"]
        if failure:
            rows["predicted"], rows["status"], rows["error"] = "", failure, "Interface preparation " + failure
        rows.to_csv(output / "predictions.csv", index=False)
        pd.DataFrame(columns=RESIDUE_COLUMNS).to_csv(output / "residues.csv", index=False)
        report = {"format": "topoppi-affinity-analysis", "version": 1,
                  "model": {"kind": "known", "target": "ddg", "training": {}},
                  "model_path": str(output / "model.joblib"),
                  "input": {"path": str(FIXTURE), "chains_a": "A", "chains_b": "B", "temperature": 298.15,
                            "wt_kd": 1e-8, "source": "user", "complex_id": "tiny"},
                  "files": {"predictions": "predictions.csv", "residues": "residues.csv"},
                  "prediction_status": failure or "partial", "layout_status": "unavailable" if failure else "ok",
                  "interface_status": failure or "ok", "layout": {}}
        if failure:
            report["interface_error"] = "Interface preparation " + failure
        else:
            loader = PDBLoader(FIXTURE)
            partners = []
            for chain in ("A", "B"):
                coords, atoms = loader.get_chain_data(chain)
                patch = trimesh.Trimesh(vertices=[[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]],
                                        faces=[[0, 1, 2], [0, 2, 3]], process=False)
                patch.metadata["source_atom_indices"] = np.array([0, 1, 4, 5])
                partners.append(PartnerInterface((chain,), coords, atoms, patches=[patch]))
            prepared = PreparedInterface(name, {"path": str(FIXTURE), "complex_id": "tiny"}, tuple(partners),
                                         SurfaceConfig(), TopologyConfig())
            save_interface(prepared, output / "interface.npz")
            report.update(analysis_id=name, interface=prepared.metadata())
            report["files"]["interface"] = "interface.npz"
            for side, partner in enumerate(partners):
                layout = partner.patches[0].copy()
                for key in ("uv", "uv_global"):
                    set_uv_layout(layout, layout.vertices[:, :2], key=key)
                name_atlas = f"partner_{side}.atlas.npz"
                save_atlas(output / name_atlas, [layout], prepared.visualizer(side),
                           run_metadata={"analysis_id": name, "side": side, "interface": prepared.metadata()})
                report["layout"][str(side)] = {"status": "ok", "atlas": name_atlas}
        (output / "analysis.json").write_text(json.dumps(report), encoding="utf-8")
        return load_analysis(output)
    return create


def _choose(app, variant, side=0):
    app.affinity_table.selection_set(variant)
    app.affinity_partner.current(side)
    app._affinity_row_selected()
    return app.show_affinity_interface()


def test_main_window_reveals_affinity_section_when_a_result_arrives(desktop, analysis_factory):
    app, errors = desktop
    assert not app.affinity_frame.winfo_manager()
    assert app.btn_binding_analysis.winfo_manager() == "pack"
    assert {"absolute", "mutation"}.issubset(app._binding_tasks)
    app.receive_affinity_result(analysis_factory())
    app.root.update()
    assert app.affinity_frame.winfo_manager() == "pack"
    children = app.right_frame.pack_slaves()
    assert children.index(app.affinity_frame) < children.index(app.canvas_frame)
    assert app.var_show_affinity.get()
    assert app.affinity_table.get_children()
    assert not errors


@pytest.mark.parametrize("has_atlas", [False, True])
def test_detailed_review_retains_prediction_and_uses_available_view(desktop, analysis_factory, monkeypatch, has_atlas):
    app, errors = desktop
    document = analysis_factory()
    if not has_atlas:
        document.report["layout"] = {}
    panel = app.open_extension("affinity")
    panel.document = document
    panel.table.insert("", "end", iid="joint", values=("joint", "A:G1W,B:G1A", ".75", "Complete"))
    panel.table.selection_set("joint")
    panel.side.set("1")
    sessions = []
    monkeypatch.setattr("topoppi_affinity.atlas_review.gui.ReviewWindow", lambda _parent, session: sessions.append(session))
    panel.open_linked_review()
    session, = sessions
    assert session.condition == ("A" if has_atlas else "G")
    assert session.variant_id == "joint" and session.side == 1
    assert session.selected == {"A:GLY:1", "B:GLY:1"}
    assert not errors


@pytest.mark.parametrize("failure", ["error", "cancelled"])
def test_failed_analysis_without_interface_id_retains_rows_and_previous_canvas(desktop, analysis_factory, monkeypatch, failure):
    app, errors = desktop
    app.receive_affinity_result(analysis_factory("previous"), show=True)
    previous, figure = app._successful_single_run, app.current_fig
    assert previous is not None and figure is not None
    failed = analysis_factory("failed", failure=failure)
    assert "analysis_id" not in failed.report
    monkeypatch.setattr("topoppi_affinity.integration.filedialog.askopenfilename", lambda **_: str(failed.path))
    app.open_affinity_result()
    app.root.update()
    assert app.affinity_table.get_children() == ("left", "joint", "right", "invalid")
    assert all(row["status"] == failure for row in app._affinity_snapshot["predictions"])
    assert app._successful_single_run is previous and app.current_fig is figure
    assert "disabled" in app.btn_show_affinity.state()
    assert "disabled" not in app.btn_export_affinity.state()
    assert "Interface preparation" in app.affinity_detail_var.get()
    assert "Prediction workspace" in app.affinity_detail_var.get()
    assert "Show interface" not in app.affinity_detail_var.get()
    if failure == "cancelled":
        assert app.affinity_table.item("left", "values")[-1] == "Cancelled"
    app.affinity_table.selection_set("joint")
    app.affinity_table.event_generate("<<TreeviewSelect>>")
    app.root.update()
    assert app._successful_single_run is previous and app.current_fig is figure
    assert "retry interface preparation" in app.affinity_detail_var.get()
    assert "generate 2D atlases" not in app.affinity_detail_var.get()
    assert not errors


def test_new_numeric_result_waits_for_explicit_canvas_handoff(desktop, analysis_factory):
    app, errors = desktop
    app.receive_affinity_result(analysis_factory("first"), show=True)
    previous, figure = app._successful_single_run, app.current_fig
    assert previous is not None and figure is not None
    app.receive_affinity_result(analysis_factory("next"))
    app.root.update()
    assert app._affinity_snapshot["analysis_id"] == "next"
    assert app._successful_single_run is previous and app.current_fig is figure
    assert "Show interface" in app.affinity_detail_var.get()
    app.btn_show_affinity.invoke()
    assert app._successful_single_run["manifest"]["affinity"]["analysis_id"] == "next"
    assert not errors


def test_joint_variant_and_opposite_partner_keep_prediction_identity(desktop, analysis_factory):
    app, errors = desktop
    document = analysis_factory()
    app.receive_affinity_result(document, show=True)
    assert _choose(app, "joint", 1)
    run = app._successful_single_run
    assert run["params"]["chain_a"] == "B" and run["params"]["chain_b"] == "A"
    assert run["style"]["highlight_residues"] == ["B:GLY:1"]
    assert run["manifest"]["affinity"]["variant_id"] == "joint"
    assert run["manifest"]["affinity"]["side"] == 1
    assert "Canvas: joint, partner B" in app.affinity_detail_var.get()
    assert app.affinity_table.item("joint", "values")[2] == "0.75"
    assert _choose(app, "left", 1)
    assert not app._successful_single_run["style"]["highlight_residues"]
    assert app.affinity_table.item("left", "values")[2] == "1.25"
    assert not errors


def test_unavailable_atlas_and_render_error_preserve_last_successful_map(desktop, analysis_factory, monkeypatch):
    app, errors = desktop
    document = analysis_factory()
    app.receive_affinity_result(document, show=True)
    previous, figure = app._successful_single_run, app.current_fig
    assert previous is not None and figure is not None
    document.report["layout"]["0"] = {"status": "error", "error": "Atlas solver stopped"}
    assert not app.show_affinity_interface(mode="atlas")
    assert app._successful_single_run is previous and app.current_fig is figure
    assert "Atlas solver stopped" in app.affinity_detail_var.get()
    with monkeypatch.context() as context:
        context.setattr("topoppi.visualization.visualizer.InterfaceVisualizer.plot_patches",
                        lambda *args, **kwargs: (_ for _ in ()).throw(ValueError("Style could not render")))
        assert not app.show_affinity_interface(mode="surface")
    assert app._successful_single_run is previous and app.current_fig is figure
    assert app._pending_single_run is not None
    assert any("Style could not render" in str(args) for args in errors)


def test_native_only_view_offers_figure_and_numbers_before_editable_atlas(desktop, analysis_factory):
    app, errors = desktop
    app.receive_affinity_result(analysis_factory(), show=True)
    assert not any("uv" in patch.metadata for patch in app._successful_single_run["patches"])
    assert "disabled" in app.btn_save_atlas.state()
    assert "disabled" not in app.btn_save.state()
    assert "disabled" not in app.btn_export_affinity.state()
    assert "Generate 2D" in app.affinity_detail_var.get()
    assert not errors


def test_portable_atlas_restores_embedded_predictions_without_plugin_or_analysis(desktop, analysis_factory, monkeypatch, tmp_path):
    app, errors = desktop
    document = analysis_factory()
    app.receive_affinity_result(document, variant_id="joint", side=1, show=True)
    assert app.show_affinity_interface(mode="atlas")
    before = app._successful_single_run["manifest"]["affinity"]
    original = tmp_path / "saved.atlas.npz"
    monkeypatch.setattr("topoppi.gui_app.ui_mixin.filedialog.asksaveasfilename", lambda **_: str(original))
    app.save_atlas()
    assert load_atlas(original).metadata["manifest"]["affinity"] == before
    moved = tmp_path / "relocated" / "saved.atlas.npz"
    moved.parent.mkdir()
    original.rename(moved)
    shutil.rmtree(document.path.parent)
    original_import = builtins.__import__
    def without_plugin(name, *args, **kwargs):
        if name.startswith("topoppi_affinity"):
            raise ImportError("Affinity plugin is unavailable in this viewing environment")
        return original_import(name, *args, **kwargs)
    with monkeypatch.context() as context:
        context.setattr(builtins, "__import__", without_plugin)
        context.setattr("topoppi.gui_app.ui_mixin.filedialog.askopenfilename", lambda **_: str(moved))
        app.open_atlas()
        app.root.update()
    assert app._affinity_document is None
    assert app._affinity_snapshot == before
    assert app._affinity_selected_id == "joint" and app.affinity_partner.current() == 1
    assert app.affinity_table.heading("value")["text"] == "ddG (kcal/mol)"
    assert app.affinity_table.item("joint", "values")[2] == "0.75"
    assert app.affinity_table.item("invalid", "values")[2:] == ("Unavailable", "Failed")
    assert app._successful_single_run["style"]["highlight_residues"] == ["B:GLY:1"]
    assert "disabled" in app.btn_show_affinity.state()
    assert "disabled" not in app.btn_export_affinity.state()
    target = tmp_path / "embedded.csv"
    monkeypatch.setattr("topoppi_affinity.integration.filedialog.asksaveasfilename", lambda **_: str(target))
    app.export_affinity_predictions()
    restored = pd.read_csv(target, keep_default_na=False)
    assert restored.variant_id.tolist() == ["left", "joint", "right", "invalid"]
    assert restored.error.iloc[-1] == "Mutation residue absent: A:99"
    assert not errors


def test_host_workspace_roundtrip_keeps_variant_partner_and_edited_next_inputs(desktop, analysis_factory):
    app, errors = desktop
    document = analysis_factory()
    app.receive_affinity_result(document, show=True)
    assert _choose(app, "joint", 1)
    app.open_affinity_workspace()
    panel = app._affinity_panel
    app.root.update()
    assert panel.table.selection() == ("joint",) and panel.side.get() == "1"
    panel.variants.delete("1.0", tk.END)
    panel.variants.insert("1.0", "next_run\tA:G1F")
    assert _choose(app, "right", 0)
    app.open_affinity_workspace()
    app.root.update()
    assert panel.table.selection() == ("right",) and panel.side.get() == "0"
    assert panel.variants.get("1.0", "end-1c") == "next_run\tA:G1F"
    panel.table.selection_set("joint")
    panel.side.set("1")
    panel.show_in_main()
    assert app._affinity_selected_id == "joint" and app.affinity_partner.current() == 1
    assert app._successful_single_run["manifest"]["affinity"]["variant_id"] == "joint"
    assert not errors


@pytest.mark.parametrize("saved_prior", [False, True])
def test_reuse_saved_inputs_discards_stale_overrides_and_restores_matching_prior(desktop, analysis_factory, saved_prior):
    app, errors = desktop
    document = analysis_factory()
    prior = document.path.parent / "prior.csv"
    if saved_prior:
        pd.DataFrame({"sample_id": ["left", "joint", "right"], "predicted": [1., 2., 3.],
                      "target": "ddg", "method": "RDE", "status": "ok", "error": ""}).to_csv(prior, index=False)
        document.report["files"]["prior"] = "prior.csv"
        document.path.write_text(json.dumps(document.report), encoding="utf-8")
    app._open_affinity()
    panel = app._affinity_panel
    panel.load(document.path)
    for key in ("prior_csv", "prior_repository", "prior_checkpoint", "sequence_dir"):
        panel.vars[key].set("stale/" + key)
    panel.reuse_inputs()
    assert panel.vars["prior_csv"].get() == (str(prior) if saved_prior else "")
    assert all(panel.vars[key].get() == "" for key in ("prior_repository", "prior_checkpoint", "sequence_dir"))
    assert variant_lines(panel.variants.get("1.0", "end")).variant_id.tolist() == ["left", "joint", "right", "invalid"]
    assert panel.document.predictions.equals(document.predictions)
    assert panel.vars["output"].get() == str(document.path.parent)
    assert panel.retry.get()
    assert not errors


def test_prediction_activity_copies_during_work_and_keeps_worker_traceback(desktop, monkeypatch, tmp_path):
    import threading
    import time
    from tkinter import ttk

    app, errors = desktop
    app._open_affinity()
    panel = app._affinity_panel
    panel.task.set("Mutation effects (ddG)")
    panel.mutation_context.set("known")
    panel.vars["structure"].set(str(FIXTURE))
    panel.vars["chains_a"].set("A")
    panel.vars["chains_b"].set("B")
    panel.vars["output"].set(str(tmp_path / "next-prediction"))
    release = threading.Event()
    def controlled_prediction_failure():
        if not release.wait(5):
            raise RuntimeError("Test worker was not released")
        raise ValueError("Controlled prediction failure for diagnostics")
    actions = next(child for child in panel.activity_tab.winfo_children() if isinstance(child, ttk.Frame))
    buttons = {button.cget("text"): button for button in actions.winfo_children()}
    target = tmp_path / "prediction-diagnostics.txt"
    monkeypatch.setattr("topoppi_affinity.gui.filedialog.asksaveasfilename", lambda **_: str(target))
    panel._launch(controlled_prediction_failure)
    try:
        assert panel._working and panel.worker.is_alive()
        for button in buttons.values():
            assert "disabled" not in button.state()
        buttons["Copy diagnostics"].invoke()
        busy_text = panel.window.clipboard_get()
        busy_data = json.loads(busy_text.split("\n\nActivity\n", 1)[0])
        assert busy_data["task"] == "Mutation effects (ddG)"
        assert busy_data["mutation_context"] == "known"
        assert busy_data["fields"]["structure"] == str(FIXTURE)
        assert busy_data["fields"]["output"] == str(tmp_path / "next-prediction")
        buttons["Save log..."].invoke()
        assert "Run started" in target.read_text()
    finally:
        release.set()
        deadline = time.monotonic() + 5
        while panel._working and time.monotonic() < deadline:
            app.root.update()
            time.sleep(.01)
    assert not panel._working and not panel.worker.is_alive()
    activity = panel.log_text.get("1.0", "end-1c")
    assert "Traceback (most recent call last)" in activity
    assert "controlled_prediction_failure" in activity
    assert "ValueError: Controlled prediction failure for diagnostics" in activity
    buttons["Copy diagnostics"].invoke()
    copied = panel.window.clipboard_get()
    buttons["Save log..."].invoke()
    assert target.read_text() == copied
    assert "Traceback (most recent call last)" in copied
    assert "controlled_prediction_failure" in copied
    assert not errors


def test_main_diagnostics_keeps_next_input_displayed_result_and_activity(desktop, analysis_factory, monkeypatch, tmp_path):
    app, errors = desktop
    document = analysis_factory()
    app.receive_affinity_result(document, variant_id="joint", side=1, show=True)
    app.var_input_path.set(str(tmp_path / "next-structure.pdb"))
    app.var_chain_a.set("X")
    app.var_chain_b.set("Y")
    app.log("Controlled main-window diagnostic message")
    app.copy_diagnostics()
    copied = app.root.clipboard_get()
    data = json.loads(copied.split("\n\nActivity\n", 1)[0])
    assert data["input"] == {"path": str(tmp_path / "next-structure.pdb"), "chain_a": "X", "chain_b": "Y"}
    shown = data["displayed_run"]["affinity"]
    assert shown["input"]["path"] == str(FIXTURE)
    assert shown["variant_id"] == "joint" and shown["side"] == 1
    assert data["extensions"]["affinity"]["analysis"] == str(document.path)
    assert "Controlled main-window diagnostic message" in copied
    assert data["python"] and data["executable"] and data["topoppi"]
    target = tmp_path / "main-diagnostics.txt"
    monkeypatch.setattr("topoppi_affinity.integration.filedialog.asksaveasfilename", lambda **_: str(target))
    app.save_diagnostics()
    assert target.read_text() == copied
    assert not errors


def test_variant_comparison_keeps_applied_style_and_partner_annotations(desktop, analysis_factory):
    app, errors = desktop
    app.receive_affinity_result(analysis_factory(), show=True)
    app.spin_size.set("18")
    app.var_show_mesh.set(False)
    app.annotation_values = {"A:GLY:1": 2.5}
    app.var_annotation_label.set("Measured residue value")
    app.btn_redraw.invoke()
    applied = app._successful_single_run["style"]
    assert applied["font_size"] == 18 and not applied["show_mesh"]
    app.current_fig.axes[0].view_init(elev=22, azim=43)
    app.affinity_table.selection_set("joint")
    app.affinity_table.event_generate("<<TreeviewSelect>>")
    app.root.update()
    current = app._successful_single_run["style"]
    assert current["font_size"] == 18 and not current["show_mesh"]
    assert current["surface_camera"]["elevation"] == pytest.approx(22)
    assert current["surface_camera"]["azimuth"] == pytest.approx(43)
    assert app.current_fig.axes[0].elev == pytest.approx(22)
    assert app.current_fig.axes[0].azim == pytest.approx(43)
    assert current["annotation_values"] == {"A:GLY:1": 2.5}
    assert current["annotation_label"] == "Measured residue value"
    assert current["highlight_residues"] == ["A:GLY:1"]
    assert app._successful_single_run["manifest"]["affinity"]["variant_id"] == "joint"
    app.affinity_partner.current(1)
    app.affinity_partner.event_generate("<<ComboboxSelected>>")
    app.root.update()
    assert app._successful_single_run["style"].get("annotation_values") is None
    assert app._successful_single_run["style"]["highlight_residues"] == ["B:GLY:1"]
    assert not errors


def test_affinity_view_switch_retains_selected_residues_until_prediction_identity_changes(desktop, analysis_factory):
    app, errors = desktop
    document = analysis_factory('view_switch')
    app.receive_affinity_result(document, show=True)
    app._toggle_binding_residue('A:GLY:1')
    selected = dict(app._binding_selection)
    assert selected
    assert app.show_affinity_interface(mode='atlas')
    assert app._binding_selection == selected
    assert app.show_affinity_interface(mode='surface')
    assert app._binding_selection == selected
    assert _choose(app, 'joint', 0)
    assert not app._binding_selection
    app._toggle_binding_residue('A:GLY:1')
    assert _choose(app, 'joint', 1)
    assert not app._binding_selection
    app._toggle_binding_residue('B:GLY:1')
    app.receive_affinity_result(analysis_factory('different_source'), show=True)
    assert not app._binding_selection
    assert not errors
