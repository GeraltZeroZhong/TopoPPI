"""Binding conventions: molar Kd, kcal/mol, and positive ddG for weaker binding."""

from __future__ import annotations

import numpy as np

R_KCAL = 0.00198720425864083
REFERENCE_TEMPERATURE_K = 298.15


def _finite(value, name, *, positive=False):
    array = np.asarray(value, dtype=float)
    if not np.isfinite(array).all() or (positive and np.any(array <= 0)):
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    return array


def kd_to_pkd(kd_M):
    """pKd = -log10(Kd / 1 M)."""
    return -np.log10(_finite(kd_M, "Kd (M)", positive=True))


def pkd_to_kd(pkd):
    with np.errstate(over="raise", under="raise"):
        try:
            return np.power(10.0, -_finite(pkd, "pKd"))
        except FloatingPointError as exc:
            raise ValueError("Predicted pKd exceeds representable positive molar Kd") from exc


def pkd_to_dg(pkd, temperature_K=REFERENCE_TEMPERATURE_K):
    """Convert at the supplied temperature; this is not temperature extrapolation."""
    return -R_KCAL * _finite(temperature_K, "Temperature (K)", positive=True) * np.log(10) * _finite(pkd, "pKd")


def kd_to_dg(kd_M, temperature_K=REFERENCE_TEMPERATURE_K):
    return pkd_to_dg(kd_to_pkd(kd_M), temperature_K)


def ddg_from_kd(wt_kd_M, mutant_kd_M, temperature_K):
    return pkd_to_dg(kd_to_pkd(mutant_kd_M) - kd_to_pkd(wt_kd_M), temperature_K)

