"""Physical queries on a decoded, self-contained surface carrier.

Masses are transported source areas (Å²), not areas recomputed from XYZ.
Only evaluation uses a reference carrier; queries themselves use their input.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.spatial import cKDTree

from ..chemistry import HYDROPATHY, classify_aa, is_hydrophobic
from ..schema import Surface

CHEMICAL_CLASSES = ("P", "A", "C", "X")
PROXIMITY_STATES = (
    "PP",
    "AA",
    "CC",
    "PA",
    "PC",
    "AC",
    "PX",
    "AX",
    "CX",
    "XX",
    "NO_PROXIMITY",
)


def _symmetric_composition(table: np.ndarray) -> np.ndarray:
    return np.asarray(
        [
            table[0, 0],
            table[1, 1],
            table[2, 2],
            table[0, 1] + table[1, 0],
            table[0, 2] + table[2, 0],
            table[1, 2] + table[2, 1],
            table[0, 3] + table[3, 0],
            table[1, 3] + table[3, 1],
            table[2, 3] + table[3, 2],
            table[3, 3],
            0.0,
        ]
    )


def _mass_fields(surface: Surface) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Total, H, and P/A/C/X source mass at each carrier cell."""
    mass = surface.mass
    if mass.shape != (len(surface.xyz), len(surface.amino_acids)):
        raise ValueError("Surface mass columns must match the residue amino-acid table")
    if not np.all(np.isfinite(mass.data)) or np.any(mass.data < 0):
        raise ValueError("Physical queries require finite, nonnegative source mass")
    one_hot = np.zeros((mass.shape[1], 4), dtype=float)
    h = np.zeros(mass.shape[1], dtype=float)
    for i, aa in enumerate(surface.amino_acids):
        one_hot[i, CHEMICAL_CLASSES.index(classify_aa(aa))] = 1.0
        h[i] = is_hydrophobic(aa)
    return surface.weights, np.asarray(mass @ h).ravel(), np.asarray(mass @ one_hot)


def _edges(surface: Surface, weights: np.ndarray) -> np.ndarray:
    """Physical edges, without accidental across-patch/partner connections."""
    edges = np.asarray(surface.edges, dtype=np.int64).reshape(-1, 2)
    if not len(edges):
        return edges
    if np.min(edges) < 0 or np.max(edges) >= len(weights):
        raise ValueError("Surface edge index is outside the carrier")
    a, b = edges.T
    keep = (
        (a != b)
        & (weights[a] > 0)
        & (weights[b] > 0)
        & (np.asarray(surface.side)[a] == np.asarray(surface.side)[b])
        & (np.asarray(surface.patch)[a] == np.asarray(surface.patch)[b])
    )
    return np.unique(np.sort(edges[keep], axis=1), axis=0)


def hydrophobic_organization(surface: Surface) -> dict:
    """Largest active H component / all H mass, separately for each partner.

    The denominator includes H in cells below the 0.5 occupancy threshold.
    A missing active domain gives L_H=0; reference applicability is separate.
    """
    weights, h_mass, _ = _mass_fields(surface)
    occupancy = np.divide(h_mass, weights, out=np.zeros_like(weights), where=weights > 0)
    active = (weights > 0) & (occupancy >= 0.5)
    edges = _edges(surface, weights)
    result = {}
    for side in (0, 1):
        side_mask = np.asarray(surface.side) == side
        total = float(h_mass[side_mask].sum())
        ids = np.flatnonzero(side_mask & active)
        largest = 0.0
        component_count = 0
        if len(ids):
            inverse = np.full(len(weights), -1, dtype=np.int64)
            inverse[ids] = np.arange(len(ids))
            pair = edges[(inverse[edges[:, 0]] >= 0) & (inverse[edges[:, 1]] >= 0)]
            local = inverse[pair]
            graph = csr_matrix(
                (np.ones(len(local)), (local[:, 0], local[:, 1])),
                shape=(len(ids), len(ids)),
            )
            component_count, labels = connected_components(graph, directed=False)
            largest = float(np.bincount(labels, weights=h_mass[ids]).max())
        result[str(side)] = {
            "H_total": total,
            "H_active": float(h_mass[side_mask & active].sum()),
            "H_inactive": float(h_mass[side_mask & ~active].sum()),
            "largest_H_component": largest,
            "L_H": largest / total if total > 0 else 0.0,
            "applicable": total > 0,
            "active_cells": int(len(ids)),
            "components": int(component_count),
            "no_active_H_domain": not len(ids),
        }
    return result


def proximity_composition(surface: Surface, block_size: int = 512) -> dict:
    """Area-product Gaussian proximity, kernel length 2 Å and cutoff 6 Å.

    Blocked sparse KD-tree pair enumeration avoids a dense N_A by N_B matrix.
    The eleven-state composition is partner-swap invariant; W has units Å⁴.
    """
    weights, _, chemistry_mass = _mass_fields(surface)
    sides = np.asarray(surface.side)
    left = np.flatnonzero((sides == 0) & (weights > 0))
    right = np.flatnonzero((sides == 1) & (weights > 0))
    q = np.divide(
        chemistry_mass,
        weights[:, None],
        out=np.zeros_like(chemistry_mass),
        where=weights[:, None] > 0,
    )
    table = np.zeros((4, 4), dtype=float)
    total = 0.0
    pairs = 0
    if block_size < 1:
        raise ValueError("KD-tree block size must be positive")
    if len(left) and len(right):
        # Reference dual cells (and uncollapsed native cells) carry one chemical
        # class. Accumulating their sixteen class pairs avoids a second sparse
        # matrix construction without changing the quadrature or kernel.
        pure = np.all(np.count_nonzero(chemistry_mass[np.r_[left, right]], axis=1) == 1)
        classes = np.argmax(chemistry_mass, axis=1) if pure else None
        # Even a fully dense local neighbourhood stays bounded in memory.
        block_size = min(block_size, max(1, 1_000_000 // len(right)))
        right_tree = cKDTree(np.asarray(surface.xyz)[right])
        for start in range(0, len(left), block_size):
            block = left[start : start + block_size]
            distances = cKDTree(np.asarray(surface.xyz)[block]).sparse_distance_matrix(
                right_tree,
                max_distance=6.0,
                output_type="coo_matrix",
            )
            contribution = (
                weights[block[distances.row]] * weights[right[distances.col]] * np.exp(-np.square(distances.data / 2.0))
            )
            if pure:
                a, b = block[distances.row], right[distances.col]
                values = contribution * q[a, classes[a]] * q[b, classes[b]]
                table += np.bincount(4 * classes[a] + classes[b], weights=values, minlength=16).reshape(4, 4)
            else:
                pair_weights = coo_matrix(
                    (contribution, (distances.row, distances.col)),
                    shape=(len(block), len(right)),
                ).tocsr()
                table += q[block].T @ (pair_weights @ q[right])
            total += float(contribution.sum())
            pairs += int(len(contribution))
    return _proximity_result(weights, chemistry_mass, left, right, table, total, pairs)


def _proximity_result(weights, chemistry_mass, left, right, table, total, pairs):
    composition = np.zeros(11, dtype=float)
    if total > 0:
        table /= total
        composition = _symmetric_composition(table)
    else:
        composition[-1] = 1.0
    # This control uses only whole-side mass composition, no pair distances.
    # It cannot detect NO_PROXIMITY when both sides exist but are far apart.
    product = np.zeros(11, dtype=float)
    if len(left) and len(right):
        marginal_left = chemistry_mass[left].sum(axis=0) / weights[left].sum()
        marginal_right = chemistry_mass[right].sum(axis=0) / weights[right].sum()
        product = _symmetric_composition(np.outer(marginal_left, marginal_right))
    else:
        product[-1] = 1.0
    return {
        "states": list(PROXIMITY_STATES),
        "composition": composition.tolist(),
        "W": total,
        "W_units": "angstrom^4",
        "pair_count": pairs,
        "applicable": bool(len(left) and len(right)),
        "global_composition_product": product.tolist(),
        "global_composition_product_TV": float(0.5 * np.abs(composition - product).sum()),
    }


def proximity_composition_cuda(surface: Surface, left_block=512, right_block=16384):
    """Same frozen integral, tiled float64 CUDA arithmetic; no distance sampling.

    The largest distance tile is bounded independently of
    the number of quadrature points, avoiding a full all-pairs allocation.
    """
    import torch

    weights, _, chemistry_mass = _mass_fields(surface)
    side = np.asarray(surface.side)
    left = np.flatnonzero((side == 0) & (weights > 0))
    right = np.flatnonzero((side == 1) & (weights > 0))
    if not len(left) or not len(right):
        return _proximity_result(weights, chemistry_mass, left, right, np.zeros((4, 4)), 0.0, 0)
    # Shared recentering reduces cancellation in ||a||²+||b||²−2a.b.
    xyz = np.asarray(surface.xyz, dtype=np.float64)
    xyz = xyz - xyz.mean(axis=0)
    a = torch.as_tensor(xyz[left], device="cuda")
    b = torch.as_tensor(xyz[right], device="cuda")
    ca = torch.as_tensor(chemistry_mass[left], device="cuda")
    cb = torch.as_tensor(chemistry_mass[right], device="cuda")
    wa = torch.as_tensor(weights[left], device="cuda")
    wb = torch.as_tensor(weights[right], device="cuda")
    nb = (b * b).sum(1)
    table = torch.zeros((4, 4), dtype=torch.float64, device="cuda")
    total = torch.zeros((), dtype=torch.float64, device="cuda")
    pairs = torch.zeros((), dtype=torch.int64, device="cuda")
    with torch.no_grad():
        for i in range(0, len(a), left_block):
            aa = a[i : i + left_block]
            na = (aa * aa).sum(1)[:, None]
            for j in range(0, len(b), right_block):
                bb = b[j : j + right_block]
                d2 = (na + nb[j : j + right_block][None, :] - 2 * (aa @ bb.T)).clamp_min_(0)
                near = d2 <= 36.0
                pairs += near.sum()
                kernel = d2.mul_(-0.25).exp_().mul_(near)
                table += ca[i : i + left_block].T @ (kernel @ cb[j : j + right_block])
                total += wa[i : i + left_block] @ (kernel @ wb[j : j + right_block])
    return _proximity_result(weights, chemistry_mass, left, right, table.cpu().numpy(), float(total), int(pairs))


def mass_anchors(surface: Surface) -> dict:
    """Carried area and chemical composition, never labelled geometric area."""
    weights, h_mass, chemistry_mass = _mass_fields(surface)
    result = {"mass_units": "angstrom^2", "mass_total": float(weights.sum()), "sides": {}}
    for side in (0, 1):
        mask = np.asarray(surface.side) == side
        total = float(weights[mask].sum())
        chemical = chemistry_mass[mask].sum(axis=0)
        result["sides"][str(side)] = {
            "mass": total,
            "H_mass": float(h_mass[mask].sum()),
            "chemical_mass": dict(zip(CHEMICAL_CLASSES, chemical.tolist(), strict=True)),
            "chemical_fraction": {
                label: float(value / total) if total > 0 else None
                for label, value in zip(CHEMICAL_CLASSES, chemical, strict=True)
            },
        }
    return result


def physical_targets(
    surface: Surface,
    proximity_surface: Surface | None = None,
    *,
    proximity_result: dict | None = None,
    proximity_device: str = "cpu",
) -> dict:
    """Compute targets without access to an original/reference geometry."""
    organization = hydrophobic_organization(surface)
    operator = proximity_composition_cuda if proximity_device == "cuda" else proximity_composition
    proximity = (
        proximity_result
        if proximity_result is not None
        else operator(
            surface if proximity_surface is None else proximity_surface,
        )
    )
    return {
        "O_by_side": organization,
        "O_valid": any(row["applicable"] for row in organization.values()),
        "X": proximity["composition"],
        "X_valid": proximity["applicable"],
        "X_states": proximity["states"],
        "W": proximity["W"],
        "W_units": proximity["W_units"],
        "X_pair_count": proximity["pair_count"],
        "X_global_composition_product": proximity["global_composition_product"],
        "X_global_composition_product_TV": proximity["global_composition_product_TV"],
        "anchors": mass_anchors(surface),
    }


def evaluate(
    reference: Surface,
    decoded: Surface,
    reference_x: Surface | None = None,
    decoded_x: Surface | None = None,
    reference_targets: dict | None = None,
    proximity_device: str = "cpu",
) -> dict:
    """Compare declared-domain targets; optional X quadrature is per carrier.

    Reference validity decides applicability. Decoded missing H is L_H=0;
    it is not omitted from O. A decoded partner loss gives NO_PROXIMITY.
    """
    ref = (
        physical_targets(reference, reference_x, proximity_device=proximity_device)
        if reference_targets is None
        else reference_targets
    )
    return compare_targets(ref, physical_targets(decoded, decoded_x, proximity_device=proximity_device))


def compare_targets(ref: dict, dec: dict) -> dict:
    """Compare precomputed targets; permits one reference integration per case."""
    side_errors = {
        side: abs(ref["O_by_side"][side]["L_H"] - dec["O_by_side"][side]["L_H"])
        if ref["O_by_side"][side]["applicable"]
        else None
        for side in ("0", "1")
    }
    observed = [value for value in side_errors.values() if value is not None]
    x_error = 0.5 * np.abs(np.asarray(ref["X"]) - np.asarray(dec["X"])).sum()
    result = {
        "O_error": float(np.mean(observed)) if observed else None,
        "O_valid": ref["O_valid"],
        "O_side_errors": side_errors,
        "X_error": float(np.clip(x_error, 0.0, 1.0)) if ref["X_valid"] else None,
        "X_valid": ref["X_valid"],
        "H_total_reference": {s: ref["O_by_side"][s]["H_total"] for s in ("0", "1")},
        "H_total_decoded": {s: dec["O_by_side"][s]["H_total"] for s in ("0", "1")},
        "W_reference": ref["W"],
        "W_decoded": dec["W"],
        "W_abs_error": abs(dec["W"] - ref["W"]),
        "W_rel_error": abs(dec["W"] - ref["W"]) / ref["W"] if ref["W"] > 0 else None,
        "mass_abs_error": abs(dec["anchors"]["mass_total"] - ref["anchors"]["mass_total"]),
        "reference_targets": ref,
        "decoded_targets": dec,
    }
    return result


def local_queries(
    surface: Surface,
    residues: list[str] | None = None,
    scales: tuple[float, ...] = (4, 12),
) -> dict:
    """Hydropathy about each residue's own positive-mass support.

    Distance is the minimum physical graph distance to that support; Gaussian
    length a is truncated at 3a. Missing support is uncovered, never value 0.
    Unknown-AA mass is reported and excluded from the hydropathy denominator.
    """
    scales = tuple(float(scale) for scale in scales)
    if not scales or any(not np.isfinite(scale) or scale <= 0 for scale in scales):
        raise ValueError("Query scales must be positive finite Å values")
    weights, _, _ = _mass_fields(surface)
    edges = _edges(surface, weights)
    lengths = np.linalg.norm(np.asarray(surface.xyz)[edges[:, 0]] - np.asarray(surface.xyz)[edges[:, 1]], axis=1)
    # Coincident true copies can have distance zero. A noncoincident seam edge
    # always retains its positive Euclidean length; no seam is a free shortcut.
    graph = csr_matrix(
        (
            np.concatenate((lengths, lengths)),
            (np.concatenate((edges[:, 0], edges[:, 1])), np.concatenate((edges[:, 1], edges[:, 0]))),
        ),
        shape=(len(weights), len(weights)),
    )
    known = np.asarray([aa in HYDROPATHY for aa in surface.amino_acids], dtype=float)
    values = np.asarray([HYDROPATHY.get(aa, 0.0) for aa in surface.amino_acids], dtype=float)
    known_mass = np.asarray(surface.mass @ known).ravel()
    hydro_mass = np.asarray(surface.mass @ values).ravel()
    lookup = {key: index for index, key in enumerate(surface.residue_keys)}
    requested = list(surface.residue_keys) if residues is None else list(residues)
    by_residue = {}
    covered_counts = {str(scale): 0 for scale in scales}
    columns = surface.mass.tocsc()
    for key in requested:
        index = lookup.get(key)
        if index is None:
            support = np.empty(0, dtype=int)
            support_mass = 0.0
        else:
            column = columns.getcol(index)
            support = column.indices[column.data > 0]
            support_mass = float(column.data[column.data > 0].sum())
        distances = None
        if len(support):
            distances = dijkstra(
                graph,
                directed=False,
                indices=support,
                min_only=True,
                limit=3.0 * max(scales),
            )
        entries = {}
        for scale in scales:
            label = str(scale)
            if distances is None:
                entries[label] = {
                    "value": None,
                    "covered": False,
                    "reason": "missing_support",
                    "support_mass": support_mass,
                }
                continue
            valid = np.isfinite(distances) & (distances <= 3.0 * scale)
            kernel = np.exp(-np.square(distances[valid] / scale) / 2.0)
            total = float(kernel @ weights[valid])
            known_total = float(kernel @ known_mass[valid])
            covered = known_total > 0
            entries[label] = {
                "value": float(kernel @ hydro_mass[valid] / known_total) if covered else None,
                "covered": covered,
                "reason": None if covered else "unknown_hydropathy",
                "support_mass": support_mass,
                "kernel_mass": total,
                "known_mass_fraction": known_total / total if total > 0 else None,
            }
            covered_counts[label] += int(covered)
        by_residue[key] = entries
    return {
        "scales": list(scales),
        "queries": by_residue,
        "n_requested": len(requested),
        "covered_counts": covered_counts,
        "coverage": {s: n / len(requested) if requested else None for s, n in covered_counts.items()},
    }
