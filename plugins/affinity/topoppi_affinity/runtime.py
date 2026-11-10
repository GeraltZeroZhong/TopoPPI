"""Small shared helpers for interactive analysis execution."""
from __future__ import annotations

import importlib.util


def check_cancel(cancel_event):
    if cancel_event is not None and cancel_event.is_set():
        raise InterruptedError("Analysis cancelled")


def resolve_device(device="auto", *, neural=True):
    requested = str(device).strip().lower()
    if requested == "cpu":
        return "cpu"
    if requested == "auto" and not neural:
        return "cpu"
    if importlib.util.find_spec("torch") is None:
        raise ImportError("Neural prediction requires the plugin's [models] dependencies; run topoppi-install-affinity")
    import torch

    if requested == "auto":
        return "cuda" if torch.cuda.is_available() else "cpu"
    if requested == "cuda" or requested.startswith("cuda:"):
        if not torch.cuda.is_available():
            raise ValueError("CUDA is unavailable in this environment. Select auto or cpu, or install a compatible CUDA PyTorch build.")
        resolved = torch.device(requested)
        if resolved.index is not None and resolved.index >= torch.cuda.device_count():
            raise ValueError(f"CUDA device {resolved.index} does not exist; {torch.cuda.device_count()} device(s) detected")
        return str(resolved)
    raise ValueError("Device must be auto, cpu, cuda or cuda:N; auto uses available CUDA, otherwise CPU")
