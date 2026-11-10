"""Exploratory organization of experimentally measured interface hotspots.

Nulls shuffle paired hotspot/energy labels within a complex and partner side,
conditional on the experimentally measured sites eligible for each metric.
Unmeasured residues are never assigned negative labels. Surface islands are
therefore islands of the *observed* hotspots, not a complete interface map.
Spatial adjacency alone does not establish energetic cooperativity.
"""
from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import numpy as np
import pandas as pd

METRIC_DEFINITIONS = {
    "hotspot_nearest_ca_distance_angstrom": "Mean distance to the nearest other measured hotspot with an observed CA; lower means clustering.",
    "hotspot_ca10_edge_count": "Undirected same-side hotspot pairs within 10 angstrom CA distance; higher means clustering.",
    "hotspot_footprint_island_count": "Components of measured hotspot residues joined by positive native footprint boundary length; lower means clustering.",
    "hotspot_footprint_largest_island_fraction": "Fraction of measured footprint hotspots in the largest induced native-footprint component.",
    "hotspot_contact_degree_difference": "Mean cross-interface contact degree of measured hotspots minus that of measured non-hotspots.",
    "hotspot_contact_betweenness_difference": "Mean native contact-network betweenness of measured hotspots minus that of measured non-hotspots.",
    "hotspot_footprint_area_share": "Area assigned to measured hotspots divided by area assigned to all measured footprint residues.",
    "positive_er_hhi": "Sum of squared positive measured Er shares; descriptive because relabeling leaves it unchanged.",
    "positive_er_footprint_area_concentration": "Sum of squared positive measured Er shares divided by measured footprint-area shares; 1 means energy proportional to area.",
}


def _bh(p_values):
    """Benjamini-Hochberg across the finite p-values in one analysis call."""
    p = np.asarray(p_values, dtype=float)
    result = np.full(len(p), np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    order = valid[np.argsort(p[valid], kind="stable")]
    if len(order):
        adjusted = p[order] * len(order) / np.arange(1, len(order) + 1)
        result[order] = np.minimum(1., np.minimum.accumulate(adjusted[::-1])[::-1])
    return result


def _components(adjacency, hotspot):
    remaining = set(np.flatnonzero(hotspot).tolist())
    sizes = []
    while remaining:
        todo = [remaining.pop()]
        size = 0
        while todo:
            node = todo.pop()
            size += 1
            neighbors = remaining.intersection(np.flatnonzero(adjacency[node]).tolist())
            remaining.difference_update(neighbors)
            todo.extend(neighbors)
        sizes.append(size)
    return len(sizes), max(sizes, default=0)


def _graph_arrays(graph):
    if isinstance(graph, (str, Path)):
        with np.load(graph, allow_pickle=False) as source:
            return {name: source[name] for name in source.files}
    return graph


def _footprint_adjacency(graph):
    n = len(graph["residue_keys"])
    result = np.zeros((n, n), dtype=bool)
    if "footprint_edges" not in graph:
        return None
    edges = np.asarray(graph["footprint_edges"], dtype=int).reshape(-1, 2)
    if "footprint_weights" in graph:
        edges = edges[np.asarray(graph["footprint_weights"]) > 0]
    if edges.size:
        result[edges[:, 0], edges[:, 1]] = True
        result[edges[:, 1], edges[:, 0]] = True
    np.fill_diagonal(result, False)
    return result


def _test_reason(hotspot, minimum):
    n, n_hot = len(hotspot), int(np.sum(hotspot))
    if n == 0:
        return "no_eligible_measured_interface_residues"
    if n_hot == 0:
        return "no_measured_hotspots"
    if n_hot == n:
        return "all_measured_sites_are_hotspots"
    if n_hot < 2:
        return "fewer_than_two_measured_hotspots"
    if n - n_hot < 2:
        return "fewer_than_two_measured_nonhotspots"
    if n < minimum:
        return "too_few_eligible_measured_sites"
    return ""


def analyze_organization(labels: pd.DataFrame, structures: Mapping, *,
                         n_permutations: int = 999, seed: int = 42,
                         min_measured: int = 6):
    """Return tidy metric observations and a JSON-serializable coverage summary.

    ``labels`` must have one row per (complex_id, residue_key), finite ``Er``
    in kcal/mol, and binary ``hotspot``. Replicate/condition aggregation happens
    before this function. ``structures`` maps complex IDs to native residue
    graph dictionaries or NPZ paths. Graph fields are ``residue_keys``, ``side``,
    ``positions``, ``interface_mask``, ``footprint_area``, ``contact_degree``,
    ``contact_betweenness``, and ``footprint_edges`` (E by 2 node indices).
    ``interface_mask`` is the native feature extraction's <=6 angstrom minimum
    partner-atom-distance criterion. ``ca_resolved`` marks observed CA atoms,
    excluding centroid fallback from spatial statistics.

    Each permutation retains the hotspot count and Er multiset on the metric's
    measured eligible sites. P-values are directional, with the +1 Monte Carlo
    correction and inclusive ties. BH includes every testable complex-side
    metric in this call. Null quantiles are permutation reference intervals,
    not confidence intervals on the observed effect. No pooled residue-level
    or cross-complex significance test is performed.
    """
    required = {"complex_id", "residue_key", "Er", "hotspot"}
    missing = required - set(labels.columns)
    if missing:
        raise ValueError(f"Missing organization label columns: {sorted(missing)}")
    if labels.duplicated(["complex_id", "residue_key"]).any():
        raise ValueError("Aggregate assay replicates before organization analysis: duplicate residue labels")
    if not np.isfinite(labels["Er"].to_numpy(dtype=float)).all():
        raise ValueError("Organization requires finite measured Er labels")
    if not labels["hotspot"].isin([0, 1, False, True]).all():
        raise ValueError("Organization requires binary measured hotspot labels")
    if n_permutations < 1 or min_measured < 4:
        raise ValueError("Use at least one permutation and min_measured >= 4")

    rng = np.random.default_rng(seed)
    rows, coverage, missing_structures, unmatched = [], [], [], []
    for complex_id, group in labels.groupby("complex_id", sort=True):
        if complex_id not in structures:
            missing_structures.append(str(complex_id))
            continue
        graph = _graph_arrays(structures[complex_id])
        keys = list(map(str, graph["residue_keys"]))
        lookup = {key: i for i, key in enumerate(keys)}
        if len(keys) != len(lookup):
            raise ValueError(f"Duplicate graph residue keys for {complex_id}")
        degree = np.asarray(graph["contact_degree"], dtype=float)
        betweenness = np.asarray(graph["contact_betweenness"], dtype=float)
        area = np.asarray(graph["footprint_area"], dtype=float)
        interface = np.asarray(graph["interface_mask"], dtype=bool)
        side = np.asarray(graph["side"])
        pos = np.asarray(graph["positions"], dtype=float)
        ca_valid = np.isfinite(pos).all(axis=1) & np.asarray(graph["ca_resolved"], dtype=bool)
        footprint = _footprint_adjacency(graph)
        mapped = group[group["residue_key"].astype(str).isin(lookup)].copy()
        missing_keys = group.loc[~group["residue_key"].astype(str).isin(lookup), "residue_key"]
        unmatched.extend({"complex_id": str(complex_id), "residue_key": str(key)} for key in missing_keys)
        mapped["node"] = mapped["residue_key"].astype(str).map(lookup)

        for partner_side in np.unique(side):
            side_mask = side == partner_side
            all_side = mapped[side_mask[mapped["node"].to_numpy(dtype=int)]]
            measured = all_side[interface[all_side["node"].to_numpy(dtype=int)]].sort_values("node")
            indices = measured["node"].to_numpy(dtype=int)
            hotspot = measured["hotspot"].to_numpy(dtype=bool)
            energy = measured["Er"].to_numpy(dtype=float)
            n_interface = int(np.sum(side_mask & interface))
            interface_area = float(area[side_mask & interface].sum())
            measured_area = float(area[indices].sum())
            base = {"complex_id": str(complex_id), "side": str(partner_side),
                    "n_interface": n_interface, "n_measured": len(indices),
                    "n_hotspot": int(hotspot.sum()),
                    "measured_site_fraction": len(indices) / n_interface if n_interface else np.nan,
                    "measured_footprint_area_fraction": measured_area / interface_area if interface_area > 0 else np.nan}
            coverage.append({**base,
                             "n_measured_noninterface_excluded": len(all_side) - len(measured),
                             "n_measured_missing_ca": int(np.sum(~ca_valid[indices])),
                             "n_measured_zero_footprint": int(np.sum(area[indices] <= 0))})

            def add_metric(name, eligible, statistic, *, direction="greater", unavailable="", descriptive=False,
                           _hotspot=hotspot, _energy=energy, _base=base):
                eligible = np.asarray(eligible, dtype=bool)
                h, e = _hotspot[eligible], _energy[eligible]
                value = float(statistic(h, e)) if not unavailable else np.nan
                reason = unavailable or ("relabeling_invariant_descriptive_metric" if descriptive else _test_reason(h, min_measured))
                if not np.isfinite(value) and not reason:
                    reason = "undefined_observed_statistic"
                row = {**_base, "metric": name, "n_eligible": len(h),
                       "n_eligible_hotspot": int(h.sum()), "value": value,
                       "alternative": "descriptive" if descriptive else direction,
                       "null_mean": np.nan, "observed_minus_null_mean": np.nan,
                       "null_q025": np.nan, "null_q975": np.nan,
                       "p_value": np.nan, "q_value_bh": np.nan,
                       "n_permutations": 0, "test_status": "not_tested", "reason": reason}
                if not reason:
                    null = np.empty(n_permutations)
                    for iteration in range(n_permutations):
                        order = rng.permutation(len(h))
                        null[iteration] = statistic(h[order], e[order])
                    if not np.isfinite(null).all():
                        row["reason"] = "undefined_permutation_statistic"
                    else:
                        row.update(null_mean=float(null.mean()), observed_minus_null_mean=float(value - null.mean()),
                                   null_q025=float(np.quantile(null, .025)), null_q975=float(np.quantile(null, .975)),
                                   n_permutations=n_permutations)
                        if np.allclose(null, value, rtol=1e-12, atol=1e-12):
                            row["reason"] = "constant_statistic_on_sampled_null"
                        else:
                            extreme = null <= value + 1e-12 if direction == "less" else null >= value - 1e-12
                            row.update(p_value=float((1 + extreme.sum()) / (n_permutations + 1)),
                                       test_status="tested", reason="")
                rows.append(row)

            ca = ca_valid[indices]
            ca_positions = pos[indices[ca]]
            distances = np.linalg.norm(ca_positions[:, None] - ca_positions[None, :], axis=-1)
            np.fill_diagonal(distances, np.inf)
            ca_edges = np.triu(distances <= 10., k=1)

            def nearest(h, _e, _distances=distances):
                if h.sum() < 2:
                    return np.nan
                return np.min(_distances[np.ix_(h, h)], axis=1).mean()

            add_metric("hotspot_nearest_ca_distance_angstrom", ca, nearest, direction="less")
            add_metric("hotspot_ca10_edge_count", ca, lambda h, e, _edges=ca_edges: np.sum(_edges & h[:, None] & h[None, :]))
            on_footprint = (area[indices] > 0) & np.isfinite(area[indices])
            footprint_indices = indices[on_footprint]
            adjacency = footprint[np.ix_(footprint_indices, footprint_indices)] if footprint is not None else None
            missing_footprint = "missing_native_footprint_boundary_edges" if footprint is None else ""
            add_metric("hotspot_footprint_island_count", on_footprint,
                       lambda h, e, _adj=adjacency: _components(_adj, h)[0], direction="less", unavailable=missing_footprint)
            add_metric("hotspot_footprint_largest_island_fraction", on_footprint,
                       lambda h, e, _adj=adjacency: _components(_adj, h)[1] / h.sum() if h.sum() else np.nan,
                       unavailable=missing_footprint)
            finite_degree = np.isfinite(degree[indices])
            measured_degree = degree[indices[finite_degree]]
            add_metric("hotspot_contact_degree_difference", finite_degree,
                       lambda h, e, _values=measured_degree: _values[h].mean() - _values[~h].mean() if h.any() and (~h).any() else np.nan)
            finite_betweenness = np.isfinite(betweenness[indices])
            measured_betweenness = betweenness[indices[finite_betweenness]]
            add_metric("hotspot_contact_betweenness_difference", finite_betweenness,
                       lambda h, e, _values=measured_betweenness: _values[h].mean() - _values[~h].mean() if h.any() and (~h).any() else np.nan)
            measured_areas = area[footprint_indices]
            area_total = measured_areas.sum()
            add_metric("hotspot_footprint_area_share", on_footprint,
                       lambda h, e, _areas=measured_areas, _total=area_total: _areas[h].sum() / _total if _total > 0 else np.nan)

            def energy_hhi(_h, e):
                positive = np.maximum(e, 0.)
                return np.sum((positive / positive.sum()) ** 2) if positive.sum() > 0 else np.nan

            def area_concentration(_h, e, _areas=measured_areas, _total=area_total):
                positive = np.maximum(e, 0.)
                if positive.sum() <= 0 or _total <= 0:
                    return np.nan
                return np.sum((positive / positive.sum()) ** 2 / (_areas / _total))

            add_metric("positive_er_hhi", np.ones(len(indices), dtype=bool), energy_hhi, descriptive=True)
            add_metric("positive_er_footprint_area_concentration", on_footprint, area_concentration)

    columns = ["complex_id", "side", "n_interface", "n_measured", "n_hotspot", "measured_site_fraction",
               "measured_footprint_area_fraction", "metric", "n_eligible", "n_eligible_hotspot", "value",
               "alternative", "null_mean", "observed_minus_null_mean", "null_q025", "null_q975",
               "p_value", "q_value_bh", "n_permutations", "test_status", "reason"]
    observations = pd.DataFrame(rows, columns=columns)
    observations["q_value_bh"] = _bh(observations["p_value"])
    n_tested = int(observations["p_value"].notna().sum())
    summary = {"protocol": "measured_interface_hotspot_organization_v1",
               "seed": seed, "requested_permutations": n_permutations, "min_measured_for_test": min_measured,
               "n_label_complexes": int(labels["complex_id"].nunique()), "n_complex_sides": len(coverage),
               "n_tests": n_tested, "n_bh_q_lt_005": int((observations["q_value_bh"] < .05).sum()),
               "multiplicity_family": "all finite complex-side-metric p-values in this analysis call",
               "null": "Paired hotspot/Er labels shuffled across measured eligible residues within complex and side; missing and unmeasured sites remain excluded.",
               "inference_scope": "Exploratory within measured interface sites; assay selection can confound apparent organization. Similar complexes and metrics are not independent replications.",
               "cooperativity": "Spatial neighbors are not evidence of cooperative energy; this requires matched experimental double-mutant cycles.",
               "footprint_scope": "Islands exclude unmeasured sites; incomplete sampling can split a true island. Area concentration uses positive measured Er and is not a full-interface free-energy decomposition.",
               "metric_definitions": METRIC_DEFINITIONS, "coverage": coverage,
               "missing_structure_complexes": missing_structures, "unmapped_measured_residues": unmatched}
    return observations, summary
