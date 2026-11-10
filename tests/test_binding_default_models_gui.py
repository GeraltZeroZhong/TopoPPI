"""Real desktop journey through the three bundled CPU prediction models."""
from __future__ import annotations

import os
import socket
import time
import tkinter as tk
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("topoppi_affinity")
pytest.importorskip("topoppi_hotspot")

from topoppi_affinity.plugin import plugin as affinity_plugin
from topoppi_hotspot.plugin import plugin as hotspot_plugin

from topoppi import plugins
from topoppi.examples import copy_example
from topoppi.gui_app.app import InterfaceApp


def _button(parent, text):
    for child in parent.winfo_children():
        if child.winfo_class() in {"TButton", "Button"} and child.cget("text") == text:
            return child
        found = _button(child, text)
        if found is not None:
            return found
    return None


def _wait(root, panel):
    deadline = time.monotonic() + 60
    while panel._working and time.monotonic() < deadline:
        root.update()
        time.sleep(.02)
    root.update()
    assert not panel._working, "The local prediction worker did not complete"
    assert panel.document is not None, getattr(panel, "errors", getattr(panel, "error", None)).get()


def _capture(window, name):
    """Optionally retain screenshots for a manual desktop review."""
    destination = os.environ.get("TOPOPPI_GUI_CAPTURE_DIR")
    if not destination:
        return
    from PIL import ImageGrab

    folder = Path(destination)
    folder.mkdir(parents=True, exist_ok=True)
    window.geometry("+20+20")
    window.deiconify()
    window.lift()
    window.update()
    x, y = window.winfo_rootx(), window.winfo_rooty()
    ImageGrab.grab(bbox=(x, y, x + window.winfo_width(), y + window.winfo_height()),
                   xdisplay=os.environ.get("DISPLAY")).save(folder / f"{name}.png")


@pytest.fixture
def binding_desktop(tmp_path, monkeypatch):
    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("A graphical display is required")
    registry = plugins.PluginRegistry(entry_points=())
    registry.register_plugin(affinity_plugin)
    registry.register_plugin(hotspot_plugin)
    monkeypatch.setattr(plugins, "_default_registry", registry)
    monkeypatch.setattr(InterfaceApp, "_recent_store_path", lambda _self: tmp_path / "recent.json")
    monkeypatch.setenv("TOPOPPI_AFFINITY_CACHE", str(tmp_path / "empty-resources"))
    monkeypatch.setenv("TOPOPPI_HOTSPOT_CACHE", str(tmp_path / "empty-hotspot-resources"))

    def offline(*args, **kwargs):
        raise AssertionError("Bundled GUI prediction attempted a network connection")

    monkeypatch.setattr(socket.socket, "connect", offline)
    root.tk.call("tk", "scaling", 1.2)
    app = InterfaceApp(root)
    root.geometry("1100x760+10+10")
    root.update()
    assert set(app._binding_tasks) == {"absolute", "mutation", "hotspot"}, app._extension_errors
    try:
        yield app
    finally:
        app.close()
        try:
            root.destroy()
        except tk.TclError:
            pass


def test_defaults_from_structure_through_ranking_variants_and_saved_results(binding_desktop, tmp_path, monkeypatch):
    app = binding_desktop
    source = copy_example(tmp_path / "input", synthetic=True)
    app.entry_file.delete(0, tk.END)
    app.entry_file.insert(0, str(source))
    app.entry_file.event_generate("<Return>")
    app.var_chain_a.set("A")
    app.var_chain_b.set("B")
    app.btn_binding_analysis.invoke()
    app.root.update()
    _button(app._binding_window, "Binding affinity").invoke()
    affinity = app._affinity_panel
    assert affinity.vars["model"].get() == ""
    assert affinity.run_button.cget("text") == "Predict"
    affinity.vars["output"].set(str(tmp_path / "affinity-result"))
    affinity.run_button.invoke()
    _wait(app.root, affinity)
    assert affinity.document.predictions.status.eq("ok").all()
    assert np.isfinite(affinity.document.predictions.predicted_pkd).all()
    affinity_path = affinity.document.path
    _capture(affinity.window, "default-affinity-result")

    _button(affinity.window, "Rank hotspots").invoke()
    hotspot = app._gui_extensions["hotspot"].panel
    assert hotspot.vars["model_dir"].get() == ""
    hotspot.vars["output_dir"].set(str(tmp_path / "hotspot-result"))
    hotspot.run_button.invoke()
    _wait(app.root, hotspot)
    assert hotspot.document.metadata["model_id"] == "topoppi-hotspot-mutation-surface"
    assert hotspot.document.metadata["target"] == "ddg"
    assert hotspot.document.metadata["units"] == "kcal/mol"
    assert hotspot.document.metadata["temperature_K"] == 298.15
    assert hotspot.table.get_children() == ("A:PHE:2",)
    assert np.isfinite(float(hotspot.document.rows[0]["rank_score"]))
    assert float(hotspot.document.rows[0]["rank_score"]) == pytest.approx(
        float(hotspot.document.rows[0]["predicted_ddg_kcal_mol"]))
    hotspot_path = hotspot.document.directory / "metadata.json"
    _capture(hotspot.window, "default-hotspot-result")

    hotspot.table.selection_set("A:PHE:2")
    hotspot.table.event_generate("<<TreeviewSelect>>")
    app.root.update()
    hotspot.variants_button.invoke()
    composer_window = next(window for window in app.root.winfo_children()
                           if isinstance(window, tk.Toplevel) and window.title() == "Create variants - Binding analysis")
    _capture(composer_window, "hotspot-to-variant")
    _button(composer_window, "Add to mutation analysis").invoke()
    assert affinity.task.get() == "Mutation effects (ddG)"
    assert "A:F2A" in affinity.variants.get("1.0", "end-1c")
    assert affinity.mutation_context.get() == ""
    affinity.vars["output"].set(str(tmp_path / "mutation-result"))
    affinity.run_button.invoke()
    _wait(app.root, affinity)
    assert affinity.document.predictions.status.eq("ok").all()
    assert np.isfinite(affinity.document.predictions.predicted_ddg_kcal_mol).all()
    mutation_path = affinity.document.path
    _capture(affinity.window, "default-mutation-result")

    # A saved native interface and its numerical output remain usable after the
    # original input is removed. Only file dialogs are supplied by this test.
    source.unlink()
    for task, saved in (("affinity", affinity_path), ("mutation", mutation_path)):
        monkeypatch.setattr("topoppi_affinity.gui.filedialog.askopenfilename", lambda saved=saved, **_: str(saved))
        affinity.open_button.invoke()
        app.root.update()
        assert affinity.document.path == saved
        assert affinity.document.interface is not None
        exported = tmp_path / f"{task}-export.csv"
        monkeypatch.setattr("topoppi_affinity.gui.filedialog.asksaveasfilename", lambda exported=exported, **_: str(exported))
        affinity.export_predictions()
        assert exported.is_file() and "predicted" in exported.read_text()
    monkeypatch.setattr("topoppi_hotspot.gui.filedialog.askopenfilename", lambda **_: str(hotspot_path))
    hotspot.open_button.invoke()
    app.root.update()
    assert hotspot.document.geometry is not None
    exported = tmp_path / "hotspot-export.csv"
    monkeypatch.setattr("topoppi_hotspot.gui.filedialog.asksaveasfilename", lambda **_: str(exported))
    hotspot.export(annotations=True)
    assert exported.read_text().startswith("residue,value\n")
