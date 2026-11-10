"""Deprecated affinity installer forwarding API.

Install ``topoppi-affinity`` first. Task setup now belongs to that package.
This module is retained solely for explicitly requested legacy imports.
"""
from __future__ import annotations

import importlib
import sys
import warnings


def _implementation():
    warnings.warn("topoppi.install_affinity is deprecated; use topoppi_affinity.installer",
                  DeprecationWarning, stacklevel=3)
    try:
        return importlib.import_module("topoppi_affinity.installer")
    except ModuleNotFoundError as exc:
        if exc.name == "topoppi_affinity":
            raise RuntimeError("Install the optional plugin first: python -m pip install topoppi-affinity") from exc
        raise


def __getattr__(name):
    return getattr(_implementation(), name)


def main(argv=None):
    try:
        implementation = _implementation()
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return implementation.main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
