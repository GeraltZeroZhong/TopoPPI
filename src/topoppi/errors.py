"""Project-specific exceptions for TopoPPI."""


class TopoPPIError(Exception):
    """Base exception for expected TopoPPI runtime failures."""


class ConfigurationError(TopoPPIError):
    """Raised when user-provided configuration is invalid."""


class InputDataError(TopoPPIError):
    """Raised when input structures or chain selections are invalid."""


class PipelineError(TopoPPIError):
    """Raised when a pipeline stage fails."""


class ExternalToolError(TopoPPIError):
    """Raised when a required external tool is missing or fails."""


class PluginError(TopoPPIError):
    """Base exception for an optional plugin boundary."""


class PluginUnavailableError(PluginError):
    """Raised when an installed plugin does not provide the requested operation."""


class PluginLoadError(PluginError):
    """Raised when a plugin cannot be registered or its implementation loaded."""


class PluginDependencyError(PluginLoadError):
    """Raised when an optional plugin dependency cannot be imported."""


class PluginExecutionError(PluginError):
    """Raised when a selected plugin operation fails."""


class ExecutionCancelled(TopoPPIError):
    """Raised when an execution context receives a cancellation request."""
