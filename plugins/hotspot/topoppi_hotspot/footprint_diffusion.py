"""Fixed area-mass graph heat and wavelet-modulus residue features.

These are residue graph operators, not mesh Laplace--Beltrami discretizations.
All WT <=6 A interface nodes participate, including unmeasured/native alanine
residues. Each complex/partner has an independent operator. No labels are read.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .features import NEIGHBOR_FEATURES

KEYS = ["complex_id", "residue_key"]
TIMES = (.25, 1., 4.)
KINDS = ("native", "spatial", "native_unweighted")


def _interface_rows(features):
    if features.duplicated(KEYS).any():
        raise ValueError("WT feature rows must have unique complex/residue keys.")
    original = features.reset_index(drop=True)
    frame = original.loc[original["contact.partner_distance"] <= 6.].sort_values(KEYS)
    return frame.copy().reset_index(drop=True), frame.index.to_numpy()


@dataclass
class QueryHeatOperator:
    """One partner's eigensystem; row_indices index the collection row order."""

    row_indices: np.ndarray
    mass: np.ndarray
    active: np.ndarray
    isolated: np.ndarray
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    rate: float
    metadata: dict

    def apply(self, values, time):
        result = values.copy()
        if time == 0 or not len(self.active) or self.rate == 0:
            return result
        sqrt_mass = np.sqrt(self.mass[self.active])
        weighted = values[self.active] * sqrt_mass[:, None]
        coefficients = self.eigenvectors.T @ weighted
        coefficients *= np.exp(-time * self.eigenvalues)[:, None]
        result[self.active] = (self.eigenvectors @ coefficients) / sqrt_mass[:, None]
        # Disconnected zero-degree rows are mathematically exact identities.
        result[self.isolated] = values[self.isolated]
        return result


@dataclass
class DiffusionOperators:
    """An aligned set of independent partner operators.

    ``row_keys`` identifies the canonical interface row order. ``source_rows``
    contains positional/iloc indices into the supplied all-residue table.
    ``apply(values, time)`` expects values already in this interface order; for
    an all-residue embedding array, call ``apply(embeddings[source_rows], t)``.
    Feature channels are arbitrary: geometry, scalar chemistry, or embeddings.
    """

    kind: str
    row_keys: pd.DataFrame
    source_rows: np.ndarray
    queries: tuple[QueryHeatOperator, ...]
    metadata: dict

    def apply(self, values, time):
        values = np.asarray(values, dtype=float)
        if values.ndim not in (1, 2) or len(values) != len(self.row_keys):
            raise ValueError("Diffusion values must have one row per row_keys entry.")
        if not np.isfinite(values).all() or not np.isfinite(time) or time < 0:
            raise ValueError("Diffusion requires finite features and a finite nonnegative time.")
        vector = values.ndim == 1
        matrix = values[:, None] if vector else values
        output = matrix.copy()
        for query in self.queries:
            output[query.row_indices] = query.apply(matrix[query.row_indices], time)
        return output[:, 0] if vector else output


def _edge_arrays(graph, kind):
    key = "spatial_edges" if kind == "spatial" else "footprint_edges"
    edges = np.asarray(graph[key], dtype=int).reshape(-1, 2)
    if kind == "native":
        weights = np.asarray(graph["footprint_weights"], dtype=float)
        if weights.shape != (len(edges),) or not np.isfinite(weights).all() or (weights < 0).any():
            raise ValueError("Native edges require aligned finite nonnegative boundary lengths.")
    else:
        weights = np.ones(len(edges))
    n = len(graph["residue_keys"])
    if len(edges) and (edges.min() < 0 or edges.max() >= n):
        raise ValueError("Graph edge indices must refer to graph residue_keys.")
    return edges, weights


def build_diffusion_operators(features, graphs, kind="native"):
    """Build dimensionless heat operators from native or spatial graph edges.

    Native W is physical shared footprint boundary length. Spatial W is one on
    the existing <=10 A same-partner CA edges. The optional native_unweighted
    control uses native adjacency with unit weights. M is footprint area, with
    full residue surface area used only for zero-footprint rows. Rows with both
    areas zero remain identity rows and lose incident diffusion edges.

    B=M^-1/2(D-W)M^-1/2 is divided by its mean positive diagonal separately for
    each partner. Its heat time is therefore dimensionless, not physical A^2.
    """
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}.")
    frame, source_rows = _interface_rows(features)
    records = []
    for (complex_id, side), positions in frame.groupby(["complex_id", "side"], sort=True).indices.items():
        positions = np.asarray(positions, dtype=int)
        part = frame.iloc[positions]
        graph = graphs[str(complex_id)]
        graph_keys = pd.Index(np.asarray(graph["residue_keys"]).astype(str))
        if not graph_keys.is_unique:
            raise ValueError("Graph residue_keys must be unique within a complex.")
        graph_rows = graph_keys.get_indexer(part.residue_key.astype(str))
        if (graph_rows < 0).any():
            raise ValueError("Every WT interface feature row must map to a graph residue key.")
        lookup = np.full(len(graph_keys), -1, dtype=int)
        lookup[graph_rows] = np.arange(len(part))
        edges, weights = _edge_arrays(graph, kind)
        local = lookup[edges]
        eligible = (local >= 0).all(axis=1) & (local[:, 0] != local[:, 1]) & (weights > 0)
        local, weights = local[eligible], weights[eligible]

        area = part[["surface.footprint_area", "surface.area"]].to_numpy(float)
        if not np.isfinite(area).all() or (area < 0).any():
            raise ValueError("Physical footprint and surface areas must be finite and nonnegative.")
        fallback = area[:, 0] == 0
        mass = np.where(fallback, area[:, 1], area[:, 0])
        zero_mass = mass == 0
        retained = ~zero_mass[local].any(axis=1)
        removed_edges = int((~retained).sum())
        local, weights = local[retained], weights[retained]
        adjacency = np.zeros((len(part), len(part)), dtype=float)
        if len(local):
            left, right = local.T
            if kind == "native":
                np.add.at(adjacency, (left, right), weights)
                np.add.at(adjacency, (right, left), weights)
            else:
                adjacency[left, right] = adjacency[right, left] = 1.
        degree = adjacency.sum(axis=1)
        isolated = degree == 0
        active = np.flatnonzero(mass > 0)
        laplacian = np.diag(degree) - adjacency
        if len(active):
            inv_sqrt_mass = 1. / np.sqrt(mass[active])
            symmetric = laplacian[np.ix_(active, active)] * inv_sqrt_mass[:, None] * inv_sqrt_mass[None, :]
            diagonal = np.diag(symmetric)
            rate = float(diagonal[diagonal > 0].mean()) if (diagonal > 0).any() else 0.
            if rate:
                eigenvalues, eigenvectors = np.linalg.eigh(symmetric / rate)
                eigenvalues = np.maximum(eigenvalues, 0.)
            else:
                eigenvalues, eigenvectors = np.zeros(len(active)), np.eye(len(active))
        else:
            rate, eigenvalues, eigenvectors = 0., np.empty(0), np.empty((0, 0))
        metadata = {
            "complex_id": str(complex_id), "side": side.item() if isinstance(side, np.generic) else side,
            "n_nodes": len(part), "n_edges": int(np.count_nonzero(np.triu(adjacency, 1))),
            "n_isolates": int(isolated.sum()), "n_zero_footprint": int(fallback.sum()),
            "n_full_surface_mass_fallback": int((fallback & ~zero_mass).sum()),
            "n_zero_area_identity": int(zero_mass.sum()),
            "zero_area_identity_residue_keys": part.loc[zero_mass, "residue_key"].astype(str).tolist(),
            "n_incident_edges_removed_for_zero_area": removed_edges,
            "mean_positive_diagonal_rate": rate,
            "mass_sum_A2": float(mass.sum()),
        }
        records.append(QueryHeatOperator(positions, mass, active, isolated, eigenvalues, eigenvectors, rate, metadata))
    metadata = {
        "kind": kind, "operator": "area_mass_residue_graph_heat_exponential",
        "mesh_laplace_beltrami": False, "time_units": "dimensionless",
        "rate_normalization": "divide_B_by_mean_positive_diagonal_per_partner; edgeless_query_is_identity",
        "edge_source": "spatial_edges" if kind == "spatial" else "footprint_edges",
        "edge_weight": "physical_shared_footprint_boundary_length_A" if kind == "native" else "unit",
        "spatial_cutoff_A": 10. if kind == "spatial" else None,
        "interface_cutoff_A": 6., "query_columns": ["complex_id", "side"],
        "mass": "footprint_area_A2; full_residue_surface_area_only_when_footprint_is_zero",
        "zero_total_area": "identity_row; incident_edges_removed; no_invented_area",
        "node_population": "all_WT_interface_residues_including_unmeasured_and_native_alanine",
        "source_row_mapping": "source_rows_are_iloc_positions_in_input; row_keys_is_canonical_complex_residue_order",
        "n_source_rows": len(features), "n_interface_rows": len(frame), "n_partner_queries": len(records),
        "queries": [record.metadata for record in records],
    }
    return DiffusionOperators(kind, frame[KEYS + ["side"]].copy(), source_rows, tuple(records), metadata)


def _time_name(time):
    return f"{time:g}".replace(".", "p")


def add_graph_diffusion(features, graphs, *, kinds=("native", "spatial")):
    """Return (interface_frame, feature_blocks, metadata) without fitting heads.

    Input columns remain unchanged. Each kind contributes 48 heat columns and
    32 wavelet-modulus columns over the fixed 16 NEIGHBOR_FEATURES channels.
    The optional native_unweighted kind is available as an adjacency control.
    """
    frame, _ = _interface_rows(features)
    source = frame[NEIGHBOR_FEATURES].to_numpy(float)
    extras, blocks, descriptions = {}, {}, {}
    for kind in kinds:
        if kind in descriptions:
            raise ValueError("Diffusion kinds must not repeat.")
        operators = build_diffusion_operators(features, graphs, kind=kind)
        descriptions[kind] = operators.metadata
        heat = [operators.apply(source, time) for time in TIMES]
        heat_columns, wavelet_columns = [], []
        for time, values in zip(TIMES, heat, strict=True):
            for number, channel in enumerate(NEIGHBOR_FEATURES):
                name = f"diffusion.{kind}.t{_time_name(time)}.{channel}"
                extras[name] = values[:, number]
                heat_columns.append(name)
        for lower, upper, first, second in zip(TIMES[:-1], TIMES[1:], heat[:-1], heat[1:], strict=True):
            contrast = np.abs(first - second)
            for number, channel in enumerate(NEIGHBOR_FEATURES):
                name = f"scattering.{kind}.t{_time_name(lower)}_t{_time_name(upper)}.{channel}"
                extras[name] = contrast[:, number]
                wavelet_columns.append(name)
        blocks[f"{kind}_heat"] = heat_columns
        blocks[f"{kind}_scattering"] = wavelet_columns
    overlap = set(extras).intersection(frame.columns)
    if overlap:
        raise ValueError("Input already contains the requested graph diffusion columns.")
    metadata = {
        "protocol": "fixed_footprint_graph_heat_and_wavelet_modulus_v1",
        "source_channels": list(NEIGHBOR_FEATURES), "times": list(TIMES),
        "wavelet_bands": [[TIMES[0], TIMES[1]], [TIMES[1], TIMES[2]]],
        "wavelet_definition": "abs(H_lower_time_X_minus_H_upper_time_X)",
        "learned_operator": False, "labels_read": False, "original_columns_unchanged": True,
        "scattering_claim": "fixed_graph_wavelet_modulus_adaptation_not_published_transform_reproduction",
        "operators": descriptions,
    }
    return pd.concat([frame, pd.DataFrame(extras, index=frame.index)], axis=1), blocks, metadata
