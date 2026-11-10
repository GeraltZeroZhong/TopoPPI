"""Existing POT transport solvers on native interface area quadrature."""
from __future__ import annotations

import time

import numpy as np
import ot
from scipy.spatial.distance import cdist

from topoppi_affinity.graph_features import _mass_sample


def interface_measure(graph, side, maximum):
    """Preserve native footprint area under the existing equal-mass sampler."""
    mask = ((graph["side"][graph["point_node_index"]] == side) &
            (graph["point_interface_weight"] > 0))
    indices = np.flatnonzero(mask)
    selected, mass = _mass_sample(graph["point_interface_weight"][indices], maximum)
    indices = indices[selected]
    if not len(indices):
        raise ValueError(f"Partner {side} has no native positive-area footprint quadrature")
    normal = np.asarray(graph["point_normal"][indices], dtype=float)
    lengths = np.linalg.norm(normal, axis=1)
    if np.any(lengths == 0):
        raise ValueError(f"Partner {side} has a zero footprint normal")
    return {"pos": np.asarray(graph["point_pos"][indices], dtype=float),
            "normal": normal / lengths[:, None], "mass": mass / mass.sum(),
            "area": float(mass.sum()), "n_points": len(indices)}


def symmetric_sinkhorn(a, b, cost, regularization):
    """Use POT log-domain Sinkhorn in both directions to remove finite-stop ordering."""
    started = time.monotonic()
    forward, forward_log = ot.bregman.sinkhorn_log(a, b, cost, regularization,
        numItermax=5000, stopThr=1e-9, warn=False, log=True)
    reverse, reverse_log = ot.bregman.sinkhorn_log(b, a, cost.T, regularization,
        numItermax=5000, stopThr=1e-9, warn=False, log=True)
    plan = .5 * (forward + reverse.T)
    residual = max(float(np.max(np.abs(plan.sum(1) - a))), float(np.max(np.abs(plan.sum(0) - b))))
    if not np.isfinite(plan).all() or residual > 1e-5:
        raise ValueError(f"POT Sinkhorn failed area marginals: maximum residual={residual}")
    return plan, {"seconds": time.monotonic() - started, "marginal_max_abs_error": residual,
                  "forward_iterations": forward_log["niter"], "reverse_iterations": reverse_log["niter"]}


def transport_features(graph, maximum=128):
    """Separate geometry and normal matching; no labels or mixed-cost score.

    Geometry uses squared distances in native Angstrom units, epsilon=4 A^2.
    Normal mismatch is ||n_A+n_B||^2/4, epsilon=.05 (dimensionless).
    The linear transport costs do not include the entropy regularizer.
    """
    started = time.monotonic()
    left, right = (interface_measure(graph, side, maximum) for side in (0, 1))
    a, b = left["mass"], right["mass"]
    geometry = cdist(left["pos"], right["pos"], metric="sqeuclidean")
    normals = .25 * cdist(left["normal"], -right["normal"], metric="sqeuclidean")
    geometry_plan, geometry_log = symmetric_sinkhorn(a, b, geometry, 4.)
    normal_plan, normal_log = symmetric_sinkhorn(a, b, normals, .05)
    exact_started = time.monotonic()
    # At <=128 points the author's network-simplex solver is inexpensive and
    # reports whether the entropic feature differs from exact area-constrained OT.
    exact_geometry = float(ot.emd2(a, b, geometry, numThreads=1))
    exact_normal = float(ot.emd2(a, b, normals, numThreads=1))
    exact_seconds = time.monotonic() - exact_started
    rows, columns = geometry.argmin(axis=1), geometry.argmin(axis=0)
    features = {
        "geometry.sinkhorn_linear_cost": float(np.sum(geometry_plan * geometry)),
        "geometry.exact_ot_cost": exact_geometry,
        "geometry.area_chamfer": float(.5 * (a @ geometry.min(axis=1) + b @ geometry.min(axis=0))),
        "geometry.independent_mean_cost": float(a @ geometry @ b),
        "normal.sinkhorn_linear_cost": float(np.sum(normal_plan * normals)),
        "normal.exact_ot_cost": exact_normal,
        "normal.under_geometry_plan": float(np.sum(geometry_plan * normals)),
        "normal.under_geometry_nearest": float(.5 * (a @ normals[np.arange(len(a)), rows] +
                                                       b @ normals[columns, np.arange(len(b))])),
        "area.sum": left["area"] + right["area"],
        "area.absolute_difference": abs(left["area"] - right["area"]),
    }
    metadata = {"maximum_points_per_side": maximum, "n_left": left["n_points"], "n_right": right["n_points"],
                "area_left": left["area"], "area_right": right["area"], "geometry_epsilon_A2": 4.,
                "normal_epsilon": .05, "geometry_solver": geometry_log, "normal_solver": normal_log,
                "exact_solver_seconds": exact_seconds, "total_seconds": time.monotonic() - started}
    return features, metadata
