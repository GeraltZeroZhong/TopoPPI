"""Select an installed chemistry backend without loading unrelated backends."""
from __future__ import annotations

from importlib import metadata

BACKEND_ENTRY_POINT_GROUP = "topoppi.chemistry.backends"


def _backend_entry_points():
    return metadata.entry_points(group=BACKEND_ENTRY_POINT_GROUP)


def backends(*, context=None):
    """Return installed backend names without importing chemical dependencies."""
    if context is not None:
        context.check_cancelled()
    return tuple(sorted({entry.name for entry in _backend_entry_points()}))


def detect(interface, *, backend="prolif", context=None):
    """Compute canonical records from selected Core atoms using one backend."""
    from topoppi.execution import ExecutionContext

    context = context or ExecutionContext()
    context.check_cancelled()
    entries = [entry for entry in _backend_entry_points() if entry.name == backend]
    if not entries:
        raise ValueError(f"Chemistry backend {backend!r} is not installed. Available: {', '.join(backends()) or 'none'}.")
    if len(entries) != 1:
        raise ValueError(f"Multiple installed providers register chemistry backend {backend!r}.")
    operation = entries[0].load()
    result = operation(interface, context=context)
    context.check_cancelled()
    return result
