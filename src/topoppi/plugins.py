"""Discover optional installed plugins and load their capabilities on demand.

Entry points in ``topoppi.plugins`` refer to a lightweight plugin object or a
zero-argument factory. Registration records callable import paths; scientific
implementations are imported only when their operation is selected.
"""
from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass
from importlib import import_module, metadata
from typing import Protocol

from topoppi.errors import (
    ExecutionCancelled,
    PluginDependencyError,
    PluginError,
    PluginExecutionError,
    PluginLoadError,
    PluginUnavailableError,
    TopoPPIError,
)
from topoppi.execution import ExecutionContext

ENTRY_POINT_GROUP = "topoppi.plugins"


@dataclass(frozen=True)
class PluginMetadata:
    name: str
    version: str
    description: str
    capabilities: tuple[str, ...] = ()


class TopoPPIPlugin(Protocol):
    metadata: PluginMetadata

    def register(self, registry: PluginRegistry) -> None: ...


@dataclass(frozen=True)
class GuiAction:
    name: str
    label: str
    plugin: str


@dataclass(frozen=True)
class _Handler:
    plugin: str
    target: str | Callable


def _unavailable(name: str) -> PluginUnavailableError:
    plugin = name.split(".", 1)[0]
    installation = (
        "Chemistry is included with TopoPPI. Reinstall TopoPPI "
        "(from a source checkout: python -m pip install .), "
        if plugin == "chemistry" else
        f"Install the corresponding plugin (usually: python -m pip install topoppi-{plugin}), "
    )
    return PluginUnavailableError(
        f"TopoPPI capability or plugin '{name}' is not available. "
        f"{installation}"
        "then restart TopoPPI or refresh the plugin registry."
    )


class PluginRegistry:
    """A registry with explicit local registration and lazy package discovery."""

    def __init__(self, entry_points: Iterable | None = None):
        self._provided_entry_points = entry_points
        self._discovered = False
        self._entries = {}
        self._plugins = {}
        self._failures = {}
        self._capabilities = {}
        self._commands = {}
        self._guis = {}
        self._gui_actions = {}
        self._registering = None

    def _discover(self):
        if self._discovered:
            return
        entries = (metadata.entry_points(group=ENTRY_POINT_GROUP)
                   if self._provided_entry_points is None else self._provided_entry_points)
        for entry in entries:
            if entry.name in self._entries:
                self._failures[entry.name] = PluginLoadError(
                    f"Multiple installed distributions register plugin '{entry.name}'. "
                    "Keep one provider installed for that plugin name."
                )
            else:
                self._entries[entry.name] = entry
        self._discovered = True

    def list(self) -> tuple[PluginMetadata, ...]:
        """List installed/local providers without importing registration modules.

        Before registration, version/description come from distribution metadata
        and capabilities are empty. ``info(name)`` reads the full declaration.
        """
        self._discover()
        result = []
        for name in sorted(set(self._entries) | set(self._plugins)):
            if name in self._plugins:
                result.append(self._plugins[name].metadata)
                continue
            distribution = getattr(self._entries[name], "dist", None)
            version = getattr(distribution, "version", "unknown")
            description = distribution.metadata.get("Summary", "") if distribution else ""
            result.append(PluginMetadata(name, version, description))
        return tuple(result)

    def available(self, name: str) -> bool:
        """Return whether the plugin provider is discovered, without loading it."""
        self._discover()
        return name in self._plugins or (name in self._entries and name not in self._failures)

    def _dependency_error(self, name, exc):
        entry = self._entries.get(name)
        distribution = getattr(entry, "dist", None)
        package = distribution.metadata.get("Name", name) if distribution else name
        missing = getattr(exc, "name", None)
        detail = f" '{missing}'" if missing else ""
        return PluginDependencyError(
            f"Plugin '{name}' requires an unavailable dependency{detail}. "
            f"Install the optional dependencies documented by '{package}'."
        )

    def register_plugin(self, plugin: TopoPPIPlugin) -> PluginMetadata:
        """Register a local plugin through the same contract as installed ones."""
        info = getattr(plugin, "metadata", None)
        if not isinstance(info, PluginMetadata) or not info.name:
            raise PluginLoadError("A plugin must declare PluginMetadata with a nonempty name.")
        if info.name in self._plugins:
            raise PluginLoadError(f"Plugin '{info.name}' is already registered.")
        if self._registering is not None:
            raise PluginLoadError("Plugin registration cannot start another plugin registration.")
        snapshots = tuple(mapping.copy() for mapping in (
            self._capabilities, self._commands, self._guis, self._gui_actions))
        self._registering = info.name
        try:
            plugin.register(self)
            actual = {name for name, handler in self._capabilities.items() if handler.plugin == info.name}
            if actual != set(info.capabilities):
                raise PluginLoadError(f"Plugin '{info.name}' registered capabilities that differ from its metadata.")
        except Exception:
            self._capabilities, self._commands, self._guis, self._gui_actions = snapshots
            raise
        finally:
            self._registering = None
        self._plugins[info.name] = plugin
        return info

    def _register(self, registry, name, target):
        plugin = self._registering
        if plugin is None:
            raise PluginLoadError("Register handlers inside a plugin's register(registry) method.")
        if name != plugin and not name.startswith(plugin + "."):
            raise PluginLoadError(f"Handler '{name}' must use its plugin namespace '{plugin}'.")
        if name in registry:
            raise PluginLoadError(f"Handler '{name}' is already registered.")
        if not callable(target) and not (isinstance(target, str) and ":" in target):
            raise PluginLoadError(f"Handler '{name}' must be a callable or a 'module:function' import path.")
        registry[name] = _Handler(plugin, target)

    def register_capability(self, name: str, target: str | Callable) -> None:
        self._register(self._capabilities, name, target)

    def register_command(self, name: str, target: str | Callable) -> None:
        self._register(self._commands, name, target)

    def register_gui(self, name: str, target: str | Callable, *, label: str | None = None) -> None:
        self._register(self._guis, name, target)
        self._gui_actions[name] = GuiAction(name, label or name, self._registering)

    def _load(self, name):
        self._discover()
        if name in self._failures:
            raise self._failures[name]
        if name in self._plugins:
            return self._plugins[name]
        if name not in self._entries:
            raise _unavailable(name)
        try:
            provider = self._entries[name].load()
            plugin = provider if hasattr(provider, "metadata") and not isinstance(provider, type) else provider()
            if getattr(getattr(plugin, "metadata", None), "name", None) != name:
                raise PluginLoadError(f"Entry point '{name}' must provide a plugin with that name.")
            self.register_plugin(plugin)
        except ImportError as exc:
            error = self._dependency_error(name, exc)
            self._failures[name] = error
            raise error from exc
        except PluginError as exc:
            self._failures[name] = exc
            raise
        except Exception as exc:
            error = PluginLoadError(f"Plugin '{name}' could not register: {exc}")
            self._failures[name] = error
            raise error from exc
        return plugin

    def info(self, name: str) -> PluginMetadata:
        """Load only the selected lightweight registration module."""
        return self._load(name).metadata

    def _load_registrations(self):
        for info in self.list():
            try:
                self._load(info.name)
            except PluginError:
                # A broken optional plugin must not prevent other commands or
                # GUI contributions; selecting it still reports its failure.
                continue

    def commands(self) -> tuple[str, ...]:
        self._load_registrations()
        return tuple(sorted(self._commands))

    def gui_actions(self) -> tuple[GuiAction, ...]:
        self._load_registrations()
        return tuple(self._gui_actions[name] for name in sorted(self._gui_actions))

    def _resolve(self, handler):
        if callable(handler.target):
            return handler.target
        module, attribute = handler.target.split(":", 1)
        try:
            value = import_module(module)
            for part in attribute.split("."):
                value = getattr(value, part)
        except ImportError as exc:
            raise self._dependency_error(handler.plugin, exc) from exc
        except (AttributeError, ValueError) as exc:
            raise PluginLoadError(f"Plugin '{handler.plugin}' has an invalid callable '{handler.target}'.") from exc
        if not callable(value):
            raise PluginLoadError(f"Plugin handler '{handler.target}' is not callable.")
        return value

    def _invoke(self, registry, name, args, kwargs, context, *, inject_context, command_exit=False):
        self._load(name.split(".", 1)[0])
        if name not in registry:
            raise _unavailable(name)
        handler = registry[name]
        context = context or ExecutionContext()
        try:
            with context.operation(name):
                function = self._resolve(handler)
                if inject_context:
                    kwargs = {**kwargs, "context": context}
                try:
                    return function(*args, **kwargs)
                except SystemExit as exc:
                    # argparse exits for help and usage errors. Keep those CLI
                    # statuses while allowing the execution lifecycle to finish.
                    if not command_exit:
                        raise
                    if exc.code is None or isinstance(exc.code, int):
                        return int(exc.code or 0)
                    raise PluginExecutionError(str(exc.code)) from exc
        except (TopoPPIError, ExecutionCancelled):
            raise
        except ImportError as exc:
            raise self._dependency_error(handler.plugin, exc) from exc
        except Exception as exc:
            raise PluginExecutionError(f"Plugin operation '{name}' failed: {exc}") from exc

    def run(self, capability: str, *args, context: ExecutionContext | None = None, **kwargs):
        return self._invoke(self._capabilities, capability, args, kwargs, context, inject_context=True)

    def run_command(self, name: str, argv, *, context: ExecutionContext | None = None):
        return self._invoke(self._commands, name, (argv,), {}, context, inject_context=False, command_exit=True)

    def open_gui(self, name: str, owner, *, context: ExecutionContext | None = None):
        return self._invoke(self._guis, name, (owner,), {}, context, inject_context=False)


_default_registry = None


def get_registry() -> PluginRegistry:
    global _default_registry
    if _default_registry is None:
        _default_registry = PluginRegistry()
    return _default_registry


def refresh() -> PluginRegistry:
    """Rediscover installed distributions after an explicit installation change."""
    global _default_registry
    _default_registry = PluginRegistry()
    return _default_registry


def list() -> tuple[PluginMetadata, ...]:
    return get_registry().list()


def available(name: str) -> bool:
    return get_registry().available(name)


def info(name: str) -> PluginMetadata:
    return get_registry().info(name)


def run(capability: str, *args, context: ExecutionContext | None = None, **kwargs):
    return get_registry().run(capability, *args, context=context, **kwargs)


def gui_actions() -> tuple[GuiAction, ...]:
    return get_registry().gui_actions()


def open_gui(name: str, owner, *, context: ExecutionContext | None = None):
    return get_registry().open_gui(name, owner, context=context)
