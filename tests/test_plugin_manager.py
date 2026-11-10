"""Plugin discovery, local installation and Core GUI boundaries."""
import sys
import tkinter as tk
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from topoppi import __version__, plugins
from topoppi.gui_app.extensions import ExtensionsMixin
from topoppi.gui_app.plugin_manager import PluginManager, installation_command, plugin_catalog


def test_plugin_catalog_describes_default_and_optional_tools():
    catalog = plugin_catalog()
    assert catalog["chemistry"]["included"]
    assert "install" not in catalog["chemistry"]
    assert "topoppi-affinity" in catalog["affinity"]["install"]
    assert "topoppi-hotspot" in catalog["hotspot"]["install"]


def test_installation_uses_running_python_and_retains_core_version(tmp_path):
    package = str(tmp_path / "plugin with spaces.whl")
    command = installation_command(package)
    assert command[:4] == [sys.executable, "-m", "pip", "install"]
    assert package in command
    assert command[-1] == f"topoppi=={__version__}"


def test_refresh_preserves_existing_controller_without_duplicate_controls(monkeypatch):
    controller = SimpleNamespace(initialize=Mock())
    host = SimpleNamespace(_gui_extensions={"example": controller}, log_history=[])
    monkeypatch.setattr(plugins, "gui_actions", lambda: [SimpleNamespace(name="example")])
    opener = Mock(side_effect=AssertionError("Existing plugin must retain its working panel"))
    monkeypatch.setattr(plugins, "open_gui", opener)
    ExtensionsMixin._init_extensions(host)
    assert host._gui_extensions["example"] is controller
    opener.assert_not_called()


def test_current_geometry_tracks_displayed_result():
    geometry = object()
    host = SimpleNamespace(_successful_single_run={"geometry": geometry})
    assert ExtensionsMixin.get_current_geometry(host) is geometry
    host._successful_single_run = {"manifest": {}}
    assert ExtensionsMixin.get_current_geometry(host) is None


@pytest.fixture
def root():
    try:
        window = tk.Tk()
    except tk.TclError as exc:
        pytest.skip(str(exc))
    window.withdraw()
    yield window
    window.destroy()


def test_manager_discovers_third_party_gui_and_blocks_open_while_mapping(root, monkeypatch):
    info = plugins.PluginMetadata("toy", "0.1", "Synthetic interface statistics.")
    monkeypatch.setattr(plugins, "list", lambda: (info,))
    monkeypatch.setattr(plugins, "info", lambda name: info)
    host = SimpleNamespace(root=root, _gui_extensions={"toy": object()}, _busy=False,
                           open_extension=Mock(), _init_extensions=Mock())
    manager = PluginManager(host)
    try:
        manager.table.selection_set("toy")
        manager._select()
        assert manager.table.item("toy", "values") == ("Toy", "Installed", "0.1")
        assert str(manager.open_button.cget("state")) == tk.NORMAL
        host._busy = True
        manager._select()
        manager.open_selected()
        host.open_extension.assert_not_called()
        host._busy = False
        manager.open_selected()
        host.open_extension.assert_called_once_with("toy")
    finally:
        manager.close()


def test_failed_installation_recovers_controls_and_exposes_log(root, monkeypatch):
    monkeypatch.setattr(plugins, "list", lambda: ())
    monkeypatch.setattr(plugins, "refresh", lambda: None)
    host = SimpleNamespace(root=root, _gui_extensions={}, _busy=False, _init_extensions=Mock())
    manager = PluginManager(host)
    try:
        manager.window.after_cancel(manager.poll_id)
        manager.working = True
        manager.events.put(("line", "Choose a compatible local wheel.\n"))
        manager.events.put(("done", 1))
        manager._poll()
        assert not manager.working
        assert "Installation needs attention" in manager.status.get()
        assert "compatible local wheel" in manager.log.get("1.0", tk.END)
        assert str(manager.local_button.cget("state")) == tk.NORMAL
    finally:
        manager.close()
