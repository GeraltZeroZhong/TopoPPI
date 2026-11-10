"""Independent plugin providers exercise discovery, lazy loading, and boundaries."""
from __future__ import annotations

import sys
from types import SimpleNamespace

import pytest

from topoppi import plugins
from topoppi.errors import (
    ExecutionCancelled,
    PluginDependencyError,
    PluginExecutionError,
    PluginLoadError,
    PluginUnavailableError,
)
from topoppi.execution import ExecutionContext
from topoppi.plugins import PluginMetadata, PluginRegistry


class EntryPoint:
    def __init__(self, plugin, *, name=None):
        self.plugin = plugin
        self.name = name or plugin.metadata.name
        self.dist = SimpleNamespace(version="1.0", metadata={"Name": "toy-plugin", "Summary": "Toy installed plugin"})
        self.loads = 0

    def load(self):
        self.loads += 1
        if isinstance(self.plugin, Exception):
            raise self.plugin
        return self.plugin


class ToyPlugin:
    metadata = PluginMetadata("toy", "1.0", "Independent toy plugin", ("toy.double",))

    def register(self, registry):
        registry.register_capability("toy.double", "toy_plugin_implementation:double")
        registry.register_command("toy", "toy_plugin_implementation:main")
        registry.register_gui("toy", "toy_plugin_implementation:panel", label="Toy analysis")


@pytest.fixture
def toy(monkeypatch):
    module = SimpleNamespace(
        double=lambda value, *, context: value * context.config_for("toy").get("factor", 2),
        main=lambda argv: 7 if argv == ["check"] else 0,
        panel=lambda owner: {"owner": owner},
    )
    monkeypatch.setitem(sys.modules, "toy_plugin_implementation", module)
    entry = EntryPoint(ToyPlugin())
    return entry, PluginRegistry([entry]), module


def test_discovery_and_available_do_not_load_provider_or_implementation():
    entry = EntryPoint(ToyPlugin())
    registry = PluginRegistry([entry])
    assert registry.available("toy")
    assert not registry.available("absent")
    assert registry.list() == (PluginMetadata("toy", "1.0", "Toy installed plugin"),)
    assert entry.loads == 0
    assert "toy_plugin_implementation" not in sys.modules
    assert registry.info("toy").capabilities == ("toy.double",)
    assert entry.loads == 1
    assert "toy_plugin_implementation" not in sys.modules


def test_installed_provider_runs_with_context_lifecycle_and_stays_registered(toy):
    entry, registry, _ = toy
    events = []
    context = ExecutionContext(plugin_config={"toy": {"factor": 3}, "absent": "ignored"}, progress=events.append)
    assert registry.run("toy.double", 4, context=context) == 12
    assert registry.run("toy.double", 5) == 10
    assert entry.loads == 1
    assert [event["stage"] for event in events] == ["started", "completed"]


def test_local_registration_uses_same_capability_command_and_gui_contract(toy):
    _, _, module = toy
    registry = PluginRegistry([])
    registry.register_plugin(ToyPlugin())
    assert registry.run("toy.double", 3) == 6
    assert registry.commands() == ("toy",)
    assert registry.run_command("toy", ["check"]) == 7
    assert registry.gui_actions()[0].label == "Toy analysis"
    owner = object()
    assert registry.open_gui("toy", owner) == module.panel(owner)


@pytest.mark.parametrize("status", [0, 2])
def test_command_parser_exit_preserves_status_and_completes_lifecycle(toy, status):
    _, registry, module = toy
    events = []

    def command(argv):
        raise SystemExit(status)

    module.main = command
    assert registry.run_command("toy", ["--help"], context=ExecutionContext(progress=events.append)) == status
    assert [event["stage"] for event in events] == ["started", "completed"]


def test_missing_plugin_and_missing_capability_report_install_guidance(toy):
    _, registry, _ = toy
    with pytest.raises(PluginUnavailableError, match="topoppi-affinity"):
        registry.run("affinity.predict", object())
    with pytest.raises(PluginUnavailableError, match="toy.unknown"):
        registry.run("toy.unknown", object())


@pytest.mark.parametrize("name", ["chemistry", "chemistry.detect"])
def test_missing_bundled_chemistry_points_to_the_topoppi_distribution(name):
    registry = PluginRegistry([])
    with pytest.raises(PluginUnavailableError) as result:
        if "." in name:
            registry.run(name, object())
        else:
            registry.info(name)
    message = str(result.value)
    assert "Chemistry is included with TopoPPI" in message
    assert "Reinstall TopoPPI" in message
    assert "python -m pip install ." in message
    assert "topoppi-chemistry" not in message
    assert "restart TopoPPI" in message


def test_missing_dependency_is_converted_at_execution_boundary():
    registry = PluginRegistry([EntryPoint(ToyPlugin())])
    with pytest.raises(PluginDependencyError, match="toy_plugin_implementation") as result:
        registry.run("toy.double", 2)
    assert isinstance(result.value.__cause__, ModuleNotFoundError)


def test_execution_failure_and_cancellation_emit_terminal_lifecycle(toy):
    _, registry, module = toy
    events = []
    context = ExecutionContext(progress=events.append)

    def fail(value, *, context):
        raise ValueError("toy invalid input")

    module.double = fail
    with pytest.raises(PluginExecutionError, match="toy invalid input"):
        registry.run("toy.double", 2, context=context)
    assert [event["stage"] for event in events] == ["started", "failed"]
    events.clear()
    context.cancel_event.set()
    with pytest.raises(ExecutionCancelled):
        registry.run("toy.double", 2, context=context)
    assert [event["stage"] for event in events] == ["started", "cancelled"]


def test_failed_registration_rolls_back_partial_handlers_and_does_not_break_other_plugins(toy):
    _, _, _ = toy

    class Broken:
        metadata = PluginMetadata("broken", "1", "Broken toy", ("broken.first",))

        def register(self, registry):
            registry.register_capability("broken.first", lambda *, context: None)
            registry.register_command("broken", lambda argv: 0)
            raise RuntimeError("registration failure")

    registry = PluginRegistry([EntryPoint(Broken()), EntryPoint(ToyPlugin())])
    assert registry.commands() == ("toy",)
    assert registry.run("toy.double", 3) == 6
    with pytest.raises(PluginLoadError, match="registration failure"):
        registry.run("broken.first")


def test_metadata_mismatch_and_duplicate_installed_names_are_clear_errors():
    duplicate = PluginRegistry([EntryPoint(ToyPlugin()), EntryPoint(ToyPlugin())])
    assert not duplicate.available("toy")
    with pytest.raises(PluginLoadError, match="Multiple installed"):
        duplicate.info("toy")
    mismatch = PluginRegistry([EntryPoint(ToyPlugin(), name="other")])
    with pytest.raises(PluginLoadError, match="that name"):
        mismatch.info("other")


def test_public_facade_uses_standard_entrypoint_group(monkeypatch, toy):
    entry, _, _ = toy
    groups = []

    def discover(*, group):
        groups.append(group)
        return [entry]

    monkeypatch.setattr(plugins.metadata, "entry_points", discover)
    monkeypatch.setattr(plugins, "_default_registry", None)
    assert plugins.available("toy")
    assert plugins.list()[0].name == "toy"
    assert plugins.info("toy").capabilities == ("toy.double",)
    assert plugins.run("toy.double", 9) == 18
    assert groups == ["topoppi.plugins"]
