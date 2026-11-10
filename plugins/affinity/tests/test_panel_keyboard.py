"""Compact prediction windows keep keyboard editing and plot controls usable."""

import json
import os
import sys
import tkinter as tk
from pathlib import Path
from tkinter import ttk

import numpy as np
import pandas as pd
import pytest
import trimesh
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
from topoppi.config import SurfaceConfig, TopologyConfig
from topoppi.gui_app.app import InterfaceApp
from topoppi.interface import PartnerInterface, PreparedInterface
from topoppi.io.io_loader import PDBLoader

from topoppi_affinity.gui import AffinityPanel
from topoppi_affinity.workflow import AnalysisDocument

pytestmark = pytest.mark.skipif(
    sys.platform not in {"win32", "darwin"} and not os.environ.get("DISPLAY"),
    reason="A graphical display is required",
)


def _widgets(parent):
    for child in parent.winfo_children():
        yield child
        yield from _widgets(child)


def _entry(panel, name):
    return next(widget for widget in _widgets(panel.window)
                if isinstance(widget, ttk.Entry) and widget.cget("textvariable") == str(panel.vars[name]))


def _settings_canvas(panel, index):
    page = panel.window.nametowidget(panel.tabs.tabs()[index])
    return next(child for child in page.winfo_children() if isinstance(child, tk.Canvas))


@pytest.fixture
def panel(tmp_path, monkeypatch, request):
    monkeypatch.setenv("TOPOPPI_AFFINITY_CACHE", str(tmp_path / "resources"))
    root = tk.Tk()
    root.withdraw()
    root.tk.call("tk", "scaling", 1.2)
    callbacks = []
    root.report_callback_exception = lambda *exc: callbacks.append(exc)
    host = None
    if getattr(request, "param", False):
        monkeypatch.setattr(InterfaceApp, "_recent_store_path", lambda _: tmp_path / "recent.json")
        host = InterfaceApp(root)
        host._open_affinity()
        view = host._affinity_panel
    else:
        view = AffinityPanel(root)
    view.task_selector.current(1)
    view.mutation_context.set("unseen")
    view.model_source.set("Local model file")
    root.update()
    try:
        yield view
        assert not callbacks
    finally:
        if host is not None:
            host.close()
        else:
            view.close()
            root.destroy()


@pytest.mark.parametrize("geometry", ["900x620", "760x500"])
def test_focused_prediction_and_advanced_fields_are_revealed_in_compact_window(panel, geometry):
    panel.window.geometry(geometry)
    panel.window.update()
    for index, widget in ((0, _entry(panel, "structure")), (0, panel.task_selector),
                          (0, panel.mutation_chain), (0, _entry(panel, "output")),
                          (1, _entry(panel, "prior_checkpoint")), (1, _entry(panel, "layout_timeout"))):
        panel.tabs.select(index)
        canvas = _settings_canvas(panel, index)
        canvas.yview_moveto(0)
        panel.window.update()
        widget.event_generate("<FocusIn>")
        panel.window.update()
        top = widget.winfo_rooty() - canvas.winfo_rooty()
        assert 0 <= top <= canvas.winfo_height() - widget.winfo_height()


def test_paging_variant_editor_keeps_outer_settings_position(panel):
    panel.window.geometry("900x620")
    panel.variants.insert("1.0", "\n".join(f"record_{i}\tA:G1W" for i in range(40)))
    panel.variants.mark_set(tk.INSERT, "1.0")
    panel.window.update()
    panel.variants.focus_force()
    panel.variants.event_generate("<FocusIn>")
    panel.window.update()
    canvas = _settings_canvas(panel, 0)
    outer = canvas.yview()
    inner = panel.variants.yview()
    panel.variants.event_generate("<KeyPress-Next>")
    panel.window.update()
    assert panel.variants.yview()[0] > inner[0]
    assert canvas.yview() == outer
    assert panel.variants.get("1.0", "1.end") == "record_0\tA:G1W"


@pytest.fixture
def small_document(tmp_path):
    structure = Path(__file__).parent / "fixtures" / "tiny_complex.pdb"
    loader = PDBLoader(structure)
    partners = []
    for chain in ("A", "B"):
        coords, atoms = loader.get_chain_data(chain)
        mesh = trimesh.Trimesh(vertices=[[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]],
                               faces=[[0, 1, 2], [0, 2, 3]], process=False)
        mesh.metadata["source_atom_indices"] = np.array([0, 1, 4, 5])
        partners.append(PartnerInterface((chain,), coords, atoms, patches=[mesh]))
    interface = PreparedInterface("compact", {"path": str(structure)}, tuple(partners),
                                  SurfaceConfig(), TopologyConfig())
    sites = json.dumps([{"chain": "A", "position": 1, "insertion": "", "wt": "G", "mt": "W",
                         "residue_key": "A:GLY:1", "side": 0, "coverage": "footprint"}])
    predictions = pd.DataFrame([
        {"variant_id": f"record_{i}", "mutations": "A:G99W" if i == 3 else "A:G1W", "predicted": "" if i == 3 else 1.5,
         "status": "error" if i == 3 else "ok", "error": "Mutation residue absent: A:99" if i == 3 else "",
         "sites": "[]" if i == 3 else sites} for i in range(4)
    ])
    report = {"analysis_id": "compact", "input": {"path": str(structure)},
              "model": {"kind": "descriptor", "target": "ddg", "scope": "Conventional interface-descriptor regression"},
              "interface": interface.metadata(), "prediction_status": "ok", "layout_status": "not_requested", "layout": {}}
    return AnalysisDocument(tmp_path / "analysis.json", report, interface, predictions, pd.DataFrame())


@pytest.mark.parametrize("panel", [False, True], indirect=True, ids=["standalone", "host"])
@pytest.mark.parametrize("geometry", ["900x620", "760x500"])
def test_compact_result_reserves_canvas_and_toolbar_and_keeps_all_rows(panel, small_document, geometry, monkeypatch):
    monkeypatch.setattr("topoppi_affinity.workflow.load_analysis", lambda _: small_document)
    panel.window.geometry(geometry)
    panel.load(small_document.path)
    panel.window.update()
    toolbar = next(child for child in panel.plot.winfo_children() if isinstance(child, NavigationToolbar2Tk))
    canvas = panel.canvas.get_tk_widget()
    assert canvas.winfo_height() >= 70
    assert toolbar.winfo_ismapped() and toolbar.winfo_height() >= 24
    assert toolbar.winfo_rooty() + toolbar.winfo_height() <= panel.window.winfo_rooty() + panel.window.winfo_height()
    assert toolbar.winfo_rooty() >= canvas.winfo_rooty() + canvas.winfo_height()
    assert panel.table.get_children() == tuple(f"record_{i}" for i in range(4))
    panel.table.selection_set("record_3")
    panel.table.see("record_3")
    panel.table.event_generate("<<TreeviewSelect>>")
    panel.window.update()
    assert panel.table.bbox("record_3")
    assert panel.figure is not None
    if panel.host is not None:
        assert panel.show_main_button.winfo_ismapped()
    detail_toggle = next(widget for widget in _widgets(panel.window)
                         if isinstance(widget, ttk.Checkbutton) and widget.cget("text") == "Show result details")
    detail_toggle.invoke()
    panel.window.update()
    detail_label = next(widget for widget in _widgets(panel.window)
                        if isinstance(widget, ttk.Label) and widget.cget("textvariable") == str(panel.details))
    assert "Mutation residue absent: A:99" in panel.details.get()
    assert detail_label.winfo_height() >= detail_label.winfo_reqheight()
    assert detail_label.winfo_rooty() + detail_label.winfo_height() <= panel.window.winfo_rooty() + panel.window.winfo_height()
