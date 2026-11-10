"""Native interface handoffs reuse compatible geometry across installed tasks."""

import time
import tkinter as tk

import pytest

pytest.importorskip("topoppi_affinity")
pytest.importorskip("topoppi_hotspot")

from topoppi_affinity.workflow import analyze, load_analysis

from topoppi.examples import copy_example
from topoppi.gui_app.app import InterfaceApp


def test_affinity_surface_can_continue_to_hotspot_without_rebuilding_geometry(tmp_path, monkeypatch):
    source = copy_example(tmp_path / "input", synthetic=True)
    report = analyze(structure=source, chains_a="A", chains_b="B", task="absolute",
                     output=tmp_path / "affinity", device="cpu", log=lambda _: None)
    document = load_analysis(report["output_directory"])
    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("A graphical display is required")
    root.withdraw()
    monkeypatch.setattr(InterfaceApp, "_recent_store_path", lambda _: tmp_path / "recent.json")
    # Use the same plugin contract in installations that enable these providers
    # through source paths instead of distribution metadata.
    from topoppi_affinity.plugin import plugin as affinity_plugin
    from topoppi_hotspot.plugin import plugin as hotspot_plugin

    from topoppi import plugins

    registry = plugins.PluginRegistry(entry_points=())
    registry.register_plugin(affinity_plugin)
    registry.register_plugin(hotspot_plugin)
    monkeypatch.setattr(plugins, "_default_registry", registry)
    app = InterfaceApp(root)
    try:
        app.receive_affinity_result(document, show=True)
        panel = app.open_binding_task("hotspot", app.binding_setup(displayed=True))
        panel.vars["output_dir"].set(str(tmp_path / "hotspot"))

        def unexpected_rebuild(*_args, **_kwargs):
            raise AssertionError("The compatible displayed geometry should be reused")

        monkeypatch.setattr("topoppi.compute_geometry", unexpected_rebuild)
        panel.run_button.invoke()
        deadline = time.monotonic() + 20
        while panel._working and time.monotonic() < deadline:
            root.update()
            time.sleep(.01)
        assert not panel._working
        assert panel.document is not None, panel.error.get()
        assert panel.document.metadata["geometry_reused"] is True
        assert panel.document.geometry is document.interface
        assert panel.document.rows[0]["residue_key"] == "A:PHE:2"
    finally:
        app.close()
