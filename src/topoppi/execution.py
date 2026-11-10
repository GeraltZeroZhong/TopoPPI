"""Shared, explicit configuration and lifecycle for platform operations."""
from __future__ import annotations

import logging
from collections.abc import Callable, Mapping
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from threading import Event

from topoppi.errors import ConfigurationError, ExecutionCancelled


@dataclass
class ExecutionContext:
    """Configuration, logging, progress, and cooperative cancellation for a run.

    Plugin namespaces are opaque until the selected plugin requests its own
    configuration. Constructing a context creates no files or worker resources.
    """

    core_config: Mapping = field(default_factory=dict)
    plugin_config: Mapping = field(default_factory=dict)
    logger: logging.Logger = field(default_factory=lambda: logging.getLogger("topoppi"))
    cancel_event: Event = field(default_factory=Event)
    progress: Callable[[dict], None] | None = None
    output_dir: Path | None = None

    def config_for(self, plugin: str) -> Mapping:
        value = self.plugin_config.get(plugin, {})
        if not isinstance(value, Mapping):
            raise ConfigurationError(f"Configuration for plugin '{plugin}' must be a mapping.")
        return value

    def log(self, message: str) -> None:
        self.logger.info(message)

    def check_cancelled(self) -> None:
        if self.cancel_event is not None and self.cancel_event.is_set():
            raise ExecutionCancelled("TopoPPI execution was cancelled.")

    def emit(self, stage: str, **details) -> None:
        if self.progress is not None:
            self.progress({"stage": stage, **details})

    @contextmanager
    def operation(self, name: str):
        """Report one operation's lifecycle without a global event dispatcher."""
        self.emit("started", operation=name)
        try:
            self.check_cancelled()
            yield self
            self.check_cancelled()
        except (ExecutionCancelled, InterruptedError) as exc:
            self.emit("cancelled", operation=name)
            if isinstance(exc, ExecutionCancelled):
                raise
            raise ExecutionCancelled(str(exc) or "TopoPPI execution was cancelled.") from exc
        except Exception as exc:
            self.emit("failed", operation=name, error=str(exc))
            raise
        else:
            self.emit("completed", operation=name)
