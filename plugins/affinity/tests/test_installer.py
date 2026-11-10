import json
import os
import sys
import threading
import time
from pathlib import Path
from types import SimpleNamespace

import pytest
from topoppi import __version__

from topoppi_affinity import installer


def test_install_uses_current_interpreter_and_preserves_core(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(installer, "discover_plugin_source", lambda: None)
    monkeypatch.setattr(installer, "__version__", "3.0.0")

    def capture(command, progress, cancel):
        calls.append(command)
        if "--constraint" in command:
            path = Path(command[command.index("--constraint") + 1])
            assert path.read_text() == "topoppi==3.0.0\n"

    monkeypatch.setattr(installer, "_run", capture)
    installer.install_affinity(cache_dir=tmp_path, progress=lambda text: None)
    pip = next(command for command in calls if "pip" in command)
    assert pip[:4] == [sys.executable, "-m", "pip", "install"]
    assert pip[-1] == "topoppi-affinity==3.0.0"
    assert calls[-1] == [sys.executable, "-m", "topoppi_affinity", "models", "download", "--task", "absolute", "--cache-dir", str(tmp_path)]


def test_source_install_is_explicit_and_download_can_be_deferred(monkeypatch, tmp_path):
    folder = tmp_path / "plugins" / "affinity"
    folder.mkdir(parents=True)
    (folder / "pyproject.toml").write_text("[project]\nname='topoppi-affinity'\n")
    calls = []
    monkeypatch.setattr(installer, "_run", lambda command, *args: calls.append(command))
    installer.install_affinity(source=tmp_path, download=False, progress=lambda text: None)
    assert calls[-1][-1] == str(folder)
    assert not any("download" in command for command in calls)
    calls.clear()
    installer.install_affinity(source=tmp_path, download=False, tasks=("known",), progress=lambda text: None)
    assert calls[-1][-1] == f"{folder}[models]"
    assert not any("download" in command for command in calls)
    with pytest.raises(ValueError, match="No plugin pyproject"):
        installer.plugin_requirement(tmp_path / "missing")
    (folder / "pyproject.toml").write_text("[project]\nname='topoppi'\n")
    with pytest.raises(ValueError, match="not the topoppi-affinity project"):
        installer.plugin_requirement(folder)


def test_subprocess_failure_preserves_actionable_output():
    with pytest.raises(RuntimeError, match="No space left"):
        installer._run([sys.executable, "-c", "import sys;print('No space left');sys.exit(7)"],
                       lambda text: None, threading.Event())


def test_cancel_stops_running_process_before_returning(tmp_path):
    marker = tmp_path / "late.txt"
    cancelled = threading.Event()
    timer = threading.Timer(0.2, cancelled.set)
    timer.start()
    started = time.monotonic()
    try:
        with pytest.raises(installer.InstallationCancelled):
            installer._run([sys.executable, "-c", f"import time,pathlib;time.sleep(20);pathlib.Path({str(marker)!r}).touch()"],
                           lambda text: None, cancelled)
    finally:
        timer.cancel()
    assert time.monotonic() - started < 10
    assert not marker.exists()


def test_status_does_not_start_installer(monkeypatch, capsys):
    monkeypatch.setattr(installer, "install_affinity", lambda **kwargs: pytest.fail("status must be read-only"))
    assert installer.main(["--status"]) == 0
    assert sys.executable in capsys.readouterr().out


def test_restart_preserves_selected_inputs_using_supported_gui_arguments():
    from topoppi.gui import build_parser

    def value(text):
        return SimpleNamespace(get=lambda: text)

    host = SimpleNamespace(var_input_path=value("/tmp/structure with spaces.pdb"), var_chain_a=value("A,B"),
                           var_chain_b=value("C"), var_output_dir=value("/tmp/results"), entry_optcuts_bin=value("/opt/OptCuts"))
    command = installer.restart_command(host)
    assert command[:3] == [sys.executable, "-m", "topoppi.gui"]
    args = build_parser().parse_args(command[3:])
    assert (args.structure, args.chain_a, args.chain_b) == ("/tmp/structure with spaces.pdb", "A,B", "C")
    assert args.output_dir == "/tmp/results" and args.optcuts_bin == "/opt/OptCuts" and args.plugin == "affinity"


@pytest.fixture
def gui_package_requirement(monkeypatch):
    """Keep GUI worker tests independent of source-checkout discovery."""
    monkeypatch.setattr(installer, "plugin_requirement",
                        lambda source=None, **kwargs: f"topoppi-affinity=={__version__}")


@pytest.mark.skipif(not os.environ.get("DISPLAY") and sys.platform.startswith("linux"), reason="Tk requires a display")
def test_gui_install_retries_and_closes_without_background_work(monkeypatch, gui_package_requirement):
    import tkinter as tk

    calls = []

    def prepare(**kwargs):
        calls.append(kwargs)
        if len(calls) == 1:
            kwargs["on_install_started"]()
            kwargs["on_installed"]()
            raise RuntimeError("Download unavailable; retry")
        kwargs["progress"]("Models cached")

    monkeypatch.setattr(installer, "install_affinity", prepare)
    monkeypatch.setattr(installer, "installation_status", lambda: {"python": sys.executable, "core": __version__, "plugin": None})
    root = tk.Tk()
    root.withdraw()
    dialog = installer.open_installer(SimpleNamespace(root=root))
    try:
        for attempt in range(2):
            dialog.run()
            deadline = time.monotonic() + 3
            while not dialog.finished and time.monotonic() < deadline:
                root.update()
                time.sleep(0.01)
            assert dialog.finished
            assert dialog.ready == (attempt == 1)
        assert calls[0]["tasks"] == ("absolute",)
        assert calls[0]["install_packages"] and not calls[1]["install_packages"]
        assert "Models cached" in dialog.log.get("1.0", "end")
        assert dialog.close()
    finally:
        root.destroy()


@pytest.mark.skipif(not os.environ.get("DISPLAY") and sys.platform.startswith("linux"), reason="Tk requires a display")
def test_gui_restarts_owner_after_dependency_install_and_keeps_restart_on_reopen(monkeypatch, gui_package_requirement):
    import tkinter as tk
    from unittest.mock import Mock

    root = tk.Tk()
    root.withdraw()
    owner = SimpleNamespace(close=Mock(), _busy=False)
    panel = SimpleNamespace(root=root, host=owner, _busy=lambda active: None, vars={
        "structure": tk.StringVar(value="input.pdb"), "chains_a": tk.StringVar(value="A,B"), "chains_b": tk.StringVar(value="C")})
    launches = Mock()
    monkeypatch.setattr(installer.subprocess, "Popen", launches)

    def prepare(**kwargs):
        assert kwargs["cache_dir"] == "/tmp/models"
        assert not kwargs["download"]
        kwargs["on_installed"]()
        raise RuntimeError("An actionable preparation failure")

    monkeypatch.setattr(installer, "install_affinity", prepare)
    dialog = installer.open_installer(panel, download=False)
    try:
        dialog.cache.set("/tmp/models")
        dialog.run()
        deadline = time.monotonic() + 3
        while not dialog.finished and time.monotonic() < deadline:
            root.update()
            time.sleep(.01)
        assert dialog.restart_required and not dialog.ready
        launches.assert_not_called()
        assert dialog.close()
        dialog = installer.open_installer(panel, download=False)
        assert dialog.restart_required
        dialog.window.geometry("500x440")
        root.update()
        assert dialog.restart.winfo_rooty() + dialog.restart.winfo_height() <= dialog.window.winfo_rooty() + dialog.window.winfo_height()
        dialog.restart_app()
        launches.assert_called_once()
        command = launches.call_args.args[0]
        assert command[:3] == [sys.executable, "-m", "topoppi.gui"]
        assert command[command.index("--structure") + 1] == "input.pdb"
        snapshot = Path(command[command.index("--plugin-setup") + 1])
        assert json.loads(snapshot.read_text()) == {"task": "absolute"}
        snapshot.unlink()
        owner.close.assert_called_once()
    finally:
        root.destroy()


@pytest.mark.skipif(not os.environ.get("DISPLAY") and sys.platform.startswith("linux"), reason="Tk requires a display")
def test_restarted_gui_restores_input_before_opening_affinity(monkeypatch, tmp_path):
    import json

    from topoppi import gui

    captured = []
    snapshot = tmp_path / "restart.json"
    state = {"task": "known", "variants": "sample-1\tGB1A", "layout": False}
    snapshot.write_text(json.dumps(state))

    def opened(app):
        captured.append((app.var_chain_a.get(), app.var_chain_b.get(), app.entry_optcuts_bin.get()))
        app._affinity_panel = SimpleNamespace(restore_setup=captured.append)
        app.root.after(10, app.root.destroy)

    from topoppi_affinity.integration import _LegacyHostMethods
    monkeypatch.setattr(_LegacyHostMethods, "_open_affinity", opened)
    fixture = Path(__file__).parent / "fixtures" / "tiny_complex.pdb"
    assert gui.main(["--structure", str(fixture), "-A", "B", "-B", "A", "--optcuts-bin", "/tmp/solver",
                     "--affinity", "--plugin-setup", str(snapshot)]) == 0
    assert captured == [("B", "A", "/tmp/solver"), state]
    assert not snapshot.exists()


def test_local_checkout_is_discovered_and_explicit_source_wins(monkeypatch, tmp_path):
    checkout = tmp_path / "checkout"
    plugin = checkout / "plugins" / "affinity"
    plugin.mkdir(parents=True)
    (plugin / "pyproject.toml").write_text("[project]\nname = 'topoppi-affinity'\n")
    nested = checkout / "work" / "analysis"
    nested.mkdir(parents=True)
    monkeypatch.setattr(installer, "__file__", str(tmp_path / "installed" / "install_affinity.py"))
    monkeypatch.chdir(nested)
    assert installer.discover_plugin_source() == plugin
    assert installer.plugin_requirement() == str(plugin)
    other = tmp_path / "other-plugin"
    other.mkdir()
    (other / "pyproject.toml").write_text("[project]\nname='topoppi-affinity'\n")
    assert installer.plugin_requirement(other) == str(other)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(installer, "__file__", str(checkout / "src" / "topoppi" / "install_affinity.py"))
    assert installer.discover_plugin_source() == plugin


def test_unpublished_development_version_requires_source(monkeypatch):
    monkeypatch.setattr(installer, "discover_plugin_source", lambda: None)
    monkeypatch.setattr(installer, "__version__", "3.0.0.dev0")
    monkeypatch.setattr(installer, "_run", lambda *args: pytest.fail("No pip call without an available source"))
    with pytest.raises(ValueError, match="matching local plugin source"):
        installer.install_affinity(progress=lambda text: None)


def test_models_only_selects_tasks_without_pip_or_source_lookup(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(installer, "_run", lambda command, *args: calls.append(command))
    monkeypatch.setattr(installer, "plugin_requirement", lambda *args: pytest.fail("Resource retry must not inspect package source"))
    installer.install_affinity(install_packages=False, tasks=["known", "unseen", "known"], cache_dir=tmp_path,
                               progress=lambda text: None,
                               on_installed=lambda: pytest.fail("No package installation occurred"))
    assert calls == [[sys.executable, "-m", "topoppi_affinity", "models", "download",
                      "--task", "known", "--task", "unseen", "--cache-dir", str(tmp_path)]]
    with pytest.raises(ValueError, match="Choose a prediction task"):
        installer.install_affinity(install_packages=False, tasks=[])


@pytest.mark.parametrize("args,tasks,packages,download", [
    ([], None, True, True),
    (["--task", "known", "--task", "unseen", "--models-only"], ["known", "unseen"], False, True),
    (["--all"], ("absolute", "mutation"), True, True),
    (["--skip-download"], None, True, False),
])
def test_cli_task_and_preparation_modes(monkeypatch, args, tasks, packages, download):
    calls = []
    monkeypatch.setattr(installer, "install_affinity", lambda **kwargs: calls.append(kwargs))
    assert installer.main(args) == 0
    assert calls[0]["tasks"] == tasks
    assert calls[0]["install_packages"] is packages
    assert calls[0]["download"] is download


def test_package_success_remains_distinguishable_from_resource_failure(monkeypatch):
    events = []
    monkeypatch.setattr(installer, "plugin_requirement", lambda source: "local-plugin[models]")

    def run(command, *args):
        if "download" in command:
            raise RuntimeError("Requested model assets have not been published")

    monkeypatch.setattr(installer, "_run", run)
    with pytest.raises(RuntimeError, match="--models-only"):
        installer.install_affinity(progress=lambda text: None, on_install_started=lambda: events.append("started"),
                                   on_installed=lambda: events.append("installed"))
    assert events == ["started", "installed"]


def test_restart_retains_full_existing_prediction_setup():
    state = {"task": "known", "fields": {"structure": "sample.pdb", "chains_a": "A"},
             "variants": "sample\tAA10V", "layout": False, "retry": True, "analysis": None}
    command = installer.restart_command(SimpleNamespace(export_setup=lambda: state), task="unseen")
    path = Path(command[command.index("--plugin-setup") + 1])
    try:
        assert json.loads(path.read_text()) == state
    finally:
        path.unlink()


@pytest.fixture
def installer_gui(monkeypatch, gui_package_requirement):
    if not os.environ.get("DISPLAY") and sys.platform.startswith("linux"):
        pytest.skip("Tk requires a display")
    import tkinter as tk
    from unittest.mock import Mock

    root = tk.Tk()
    root.withdraw()
    monkeypatch.setattr(installer, "installation_status", lambda: {"python": sys.executable, "core": __version__, "plugin": None})
    owner = SimpleNamespace(root=root, close=Mock(), _busy=False)
    yield root, owner
    dialog = getattr(owner, "_affinity_installer", None)
    if dialog and dialog.worker and dialog.worker.is_alive():
        dialog.cancel.set()
        dialog.worker.join(timeout=2)
    if dialog and dialog.window.winfo_exists():
        dialog.close()
    root.destroy()


def _finish_gui(root, dialog):
    deadline = time.monotonic() + 3
    while (not dialog.finished or dialog.worker.is_alive()) and time.monotonic() < deadline:
        root.update()
        time.sleep(.01)
    assert dialog.finished and not dialog.worker.is_alive()
    root.update()


def test_gui_cancellation_retains_choices_and_can_retry(monkeypatch, installer_gui):
    root, host = installer_gui
    calls = []

    def prepare(**kwargs):
        calls.append(kwargs)
        if len(calls) == 1:
            kwargs["on_install_started"]()
            assert kwargs["cancel_event"].wait(timeout=2)
            raise installer.InstallationCancelled("Cancelled")
        kwargs["on_installed"]()

    monkeypatch.setattr(installer, "install_affinity", prepare)
    dialog = installer.open_installer(host, tasks=("known", "unseen"))
    dialog.cache.set("/tmp/chosen models")
    dialog.run()
    dialog.cancel_run()
    _finish_gui(root, dialog)
    assert "cancelled" in dialog.status.get()
    assert dialog.restart_required and not dialog.dependencies_installed
    assert dialog.selected_tasks() == ("known", "unseen")
    dialog.run()
    _finish_gui(root, dialog)
    assert dialog.ready and calls[1]["install_packages"]
    assert calls[1]["cache_dir"] == "/tmp/chosen models"
    assert dialog.close()
    reopened = installer.open_installer(host)
    assert reopened.selected_tasks() == ("known", "unseen")
    assert reopened.cache.get() == "/tmp/chosen models"
    assert reopened.restart_required


def test_gui_requires_a_task_only_when_resources_are_selected(monkeypatch, installer_gui):
    root, host = installer_gui
    calls = []
    monkeypatch.setattr(installer, "install_affinity", lambda **kwargs: calls.append(kwargs))
    dialog = installer.open_installer(host)
    for variable in dialog.task_vars.values():
        variable.set(False)
    dialog.run()
    assert dialog.worker is None and not calls
    assert "Choose at least one" in dialog.status.get()
    dialog.prepare_models.set(False)
    dialog.run()
    _finish_gui(root, dialog)
    assert not calls[0]["download"]


def test_gui_closing_waits_for_cancelled_worker(monkeypatch, installer_gui):
    root, host = installer_gui

    def prepare(**kwargs):
        assert kwargs["cancel_event"].wait(timeout=2)
        raise installer.InstallationCancelled("Cancelled")

    monkeypatch.setattr(installer, "install_affinity", prepare)
    dialog = installer.open_installer(host)
    dialog.run()
    assert not dialog.close()
    deadline = time.monotonic() + 3
    while dialog.window.winfo_exists() and time.monotonic() < deadline:
        root.update()
        time.sleep(.01)
    assert not dialog.worker.is_alive()
    assert not dialog.window.winfo_exists()


def test_prepared_existing_plugin_opens_without_restart(monkeypatch, installer_gui):
    root, host = installer_gui
    from unittest.mock import Mock

    monkeypatch.setattr(installer, "installation_status", lambda: {"python": sys.executable, "core": __version__, "plugin": __version__})
    monkeypatch.setattr(installer.importlib.util, "find_spec", lambda module: object())
    restore = Mock()
    host._affinity_panel = SimpleNamespace(restore_setup=restore)
    host._open_affinity = Mock()
    calls = []
    monkeypatch.setattr(installer, "install_affinity", lambda **kwargs: calls.append(kwargs))
    monkeypatch.setattr(installer.subprocess, "Popen", lambda *args: pytest.fail("Resource preparation does not require restarting"))
    dialog = installer.open_installer(host, tasks=("unseen",))
    dialog.run()
    _finish_gui(root, dialog)
    assert not calls[0]["install_packages"] and not dialog.restart_required
    dialog.continue_app()
    host._open_affinity.assert_called_once()
    restore.assert_called_once_with({"task": "unseen"})
    host.close.assert_not_called()


def test_bundled_tasks_prepare_without_neural_dependencies(monkeypatch, installer_gui):
    root, host = installer_gui
    monkeypatch.setattr(installer, "installation_status", lambda: {
        "python": sys.executable, "core": __version__, "plugin": __version__})
    monkeypatch.setattr(installer.importlib.util, "find_spec", lambda module:
                        None if module in {"torch", "esm", "transformers", "torch_geometric"} else object())
    calls = []
    monkeypatch.setattr(installer, "install_affinity", lambda **kwargs: calls.append(kwargs))
    dialog = installer.open_installer(host, tasks=("absolute", "mutation"))
    assert dialog.dependencies_installed and not dialog.needs_installation()
    assert dialog.start.cget("text") == "Prepare selected tasks"
    dialog.run()
    _finish_gui(root, dialog)
    assert not calls[0]["install_packages"] and not dialog.restart_required
    dialog.task_vars["known"].set(True)
    assert dialog.needs_installation()
    assert dialog.start.cget("text") == "Install and prepare"


@pytest.mark.parametrize("missing_module", ["torch", "esm", "transformers", "torch_geometric"])
def test_base_only_plugin_installs_prediction_dependencies(monkeypatch, installer_gui, missing_module):
    root, host = installer_gui
    monkeypatch.setattr(installer, "installation_status", lambda: {"python": sys.executable, "core": __version__, "plugin": __version__})
    checked = []

    def find_spec(module):
        checked.append(module)
        return None if module == missing_module else object()

    monkeypatch.setattr(installer.importlib.util, "find_spec", find_spec)
    calls = []
    monkeypatch.setattr(installer, "_run", lambda command, *args: calls.append(command))
    dialog = installer.open_installer(host, download=True, tasks=("known",))
    assert dialog.dependencies_installed and dialog.needs_installation()
    assert missing_module in checked
    dialog.run()
    _finish_gui(root, dialog)
    assert any(command[:4] == [sys.executable, "-m", "pip", "install"] for command in calls)
    assert calls[-1][-2:] == ["--task", "known"]
    assert dialog.restart_required and dialog.ready
