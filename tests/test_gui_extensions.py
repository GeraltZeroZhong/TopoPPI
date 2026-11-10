"""Core GUI behavior depends on extension contracts rather than scientific tasks."""
import json
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from topoppi.gui_app.extensions import ExtensionsMixin


def test_missing_extension_has_actionable_error():
    host = SimpleNamespace(_gui_extensions={})
    with pytest.raises(RuntimeError, match="Install its corresponding TopoPPI plugin"):
        ExtensionsMixin.open_extension(host, "example")


def test_extension_open_redraw_and_shutdown_use_the_same_controller():
    controller = SimpleNamespace(open=Mock(return_value="panel"), redraw=Mock(return_value=True),
                                 close=Mock(side_effect=[False, True]), restart_required=True)
    host = SimpleNamespace(_gui_extensions={"example": controller})
    setup = {"input": "synthetic"}
    assert ExtensionsMixin.open_extension(host, "example", setup) == "panel"
    controller.open.assert_called_once_with(setup=setup)
    assert ExtensionsMixin._extensions_redraw(host, {"view": "surface"})
    assert ExtensionsMixin._extensions_need_restart(host)
    assert not ExtensionsMixin._close_extensions(host)
    assert ExtensionsMixin._close_extensions(host)


def test_unknown_saved_metadata_is_readable_without_loading_a_plugin():
    host = SimpleNamespace(_gui_extensions={}, _successful_single_run=None, _pending_single_run=None)
    data = {"uninstalled_provider": {"rows": [{"id": "toy", "value": 2.5}], "units": "synthetic"}}
    ExtensionsMixin._extensions_document_loaded(host, data)
    assert json.loads(ExtensionsMixin.saved_metadata_text(host)) == data


def test_legacy_installer_explains_optional_package_requirement(monkeypatch, capsys):
    from topoppi import install_affinity

    def unavailable(name):
        raise ModuleNotFoundError("optional package absent", name="topoppi_affinity")

    monkeypatch.setattr(install_affinity.importlib, "import_module", unavailable)
    with pytest.warns(DeprecationWarning), pytest.raises(RuntimeError, match="pip install topoppi-affinity"):
        _ = install_affinity.run_setup
    with pytest.warns(DeprecationWarning):
        assert install_affinity.main(["--status"]) == 1
    assert "pip install topoppi-affinity" in capsys.readouterr().err


def test_failed_extension_initialization_leaves_core_controls_usable(monkeypatch):
    from topoppi import plugins

    host = SimpleNamespace(_gui_extensions={}, log_history=[])
    monkeypatch.setattr(plugins, "gui_actions", lambda: [SimpleNamespace(name="broken")])
    controller = SimpleNamespace(initialize=Mock(side_effect=ImportError("optional renderer unavailable")))
    monkeypatch.setattr(plugins, "open_gui", lambda *args: controller)
    ExtensionsMixin._init_extensions(host)
    ExtensionsMixin._sync_extensions(host)
    assert host._gui_extensions == {}
    assert "optional renderer unavailable" in host.log_history[0]
