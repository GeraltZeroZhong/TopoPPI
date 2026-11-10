"""Task UI compatibility and lazy fallback remain owned by the plugin."""
import sys
from types import SimpleNamespace
from unittest import mock

from topoppi_affinity.integration import AffinityController, _LegacyHostMethods


def test_affinity_entry_opens_setup_when_plugin_is_missing():
    harness = SimpleNamespace(_busy=False, _closing=False)
    with mock.patch("topoppi_affinity.installer.open_installer") as installer, mock.patch.dict(
        sys.modules, {"topoppi_affinity.gui": None}
    ):
        _LegacyHostMethods._open_affinity(harness)
    installer.assert_called_once_with(harness)


def test_affinity_entry_reuses_the_installed_panel():
    harness = SimpleNamespace(_busy=False, _closing=False)
    panel = mock.Mock()
    with mock.patch.dict(sys.modules, {"topoppi_affinity.gui": SimpleNamespace(open_panel=panel)}):
        _LegacyHostMethods._open_affinity(harness)
    panel.assert_called_once_with(harness)


def test_affinity_entry_waits_for_active_core_work():
    harness = SimpleNamespace(_busy=True, _closing=False)
    with mock.patch("topoppi_affinity.installer.open_installer") as installer:
        _LegacyHostMethods._open_affinity(harness)
    installer.assert_not_called()


def test_dependency_restart_preserves_saved_result_review():
    host = SimpleNamespace(_busy=False, _closing=False, _affinity_restart_required=True,
                           btn_affinity=mock.Mock(), analysis_menu=mock.Mock(), _sync_affinity_controls=mock.Mock())
    controller = AffinityController(host)
    controller.menu_entries = [2, 3]
    controller.sync_controls()
    host.btn_affinity.configure.assert_called_once_with(state="disabled")
    assert host.analysis_menu.entryconfigure.call_args_list == [
        mock.call(2, state="disabled"), mock.call(3, state="normal")]



def test_controller_registers_two_scientific_tasks_with_shared_navigation(monkeypatch):
    monkeypatch.setattr(_LegacyHostMethods, '_init_affinity_results', lambda self: None)
    host = SimpleNamespace(register_binding_task=mock.Mock(), analysis_menu=mock.Mock())
    controller = AffinityController(host)
    controller.initialize()
    registrations = {call.args[0]: call.kwargs for call in host.register_binding_task.call_args_list}
    assert set(registrations) == {'absolute', 'mutation'}
    assert not hasattr(host, 'btn_affinity')
    controller.open = mock.Mock()
    registrations['mutation']['open_task'](setup={'structure': 'complex.pdb'})
    controller.open.assert_called_once_with(setup={'structure': 'complex.pdb', 'prediction_task': 'mutation'})
    host.analysis_menu.add_command.assert_called_once_with(label='Open affinity result...', command=host.open_affinity_result)
