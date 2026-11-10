"""Construct relation-specific normalized operators on residue graphs."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from .partner_features import KEYS


def build_graph_operators(frame, graphs):
    """Map full-graph indices to the supplied interface frame's positional rows.

    Native and spatial relations stay within a partner. Context crosses the
    two partners of the same complex. Non-interface endpoints are discarded,
    then symmetric weighted edges are row-normalized. Empty rows remain zero.
    No label column is read; native alanines and all other frame nodes remain.
    """
    if frame.duplicated(KEYS).any():
        raise ValueError("Interface nodes must have unique complex/residue keys.")
    buckets = {kind: {"row": [], "col": [], "weight": []} for kind in ("native", "spatial", "context")}
    records = []
    for complex_id, positions in frame.groupby("complex_id", sort=True).indices.items():
        part = frame.iloc[positions]
        if "family_group" in frame and part.family_group.nunique(dropna=False) != 1:
            raise ValueError("All nodes in a complex must belong to one externally assigned family.")
        graph = graphs[str(complex_id)]
        keys = pd.Index(np.asarray(graph["residue_keys"]).astype(str))
        if not keys.is_unique:
            raise ValueError("Graph residue_keys must be unique within a complex.")
        graph_rows = keys.get_indexer(part.residue_key.astype(str))
        if np.any(graph_rows < 0):
            raise ValueError("Every interface row must map to a graph residue key.")
        node_map = np.full(len(keys), -1, dtype=int)
        node_map[graph_rows] = positions
        side = frame.side.to_numpy()
        record = {"complex_id": str(complex_id), "n_nodes": len(part), "n_full_graph_nodes": len(keys)}
        for kind, edge_key, weight_key in (
            ("native", "footprint_edges", "footprint_weights"),
            ("spatial", "spatial_edges", None),
            ("context", "cross_edges", "cross_weights"),
        ):
            edges = np.asarray(graph[edge_key], dtype=int).reshape(-1, 2)
            weights = np.ones(len(edges)) if weight_key is None else np.asarray(graph[weight_key], dtype=float)
            if weights.shape != (len(edges),) or not np.isfinite(weights).all() or np.any(weights < 0):
                raise ValueError(f"{kind} edges require aligned finite nonnegative weights.")
            if len(edges) and (edges.min() < 0 or edges.max() >= len(keys)):
                raise ValueError("Graph edge endpoints must index the full residue-key registry.")
            mapped = node_map[edges]
            keep = (mapped >= 0).all(axis=1) & (mapped[:, 0] != mapped[:, 1]) & (weights > 0)
            mapped, weights = mapped[keep], weights[keep]
            same_side = side[mapped[:, 0]] == side[mapped[:, 1]]
            if np.any(same_side if kind == "context" else ~same_side):
                raise ValueError(f"{kind} edges violate their declared partner relation.")
            left, right = mapped.T
            bucket = buckets[kind]
            bucket["row"].extend(np.concatenate([left, right]).tolist())
            bucket["col"].extend(np.concatenate([right, left]).tolist())
            bucket["weight"].extend(np.tile(weights, 2).tolist())
            record[f"n_{kind}_edges"] = len(mapped)
        records.append(record)
    operators, details = {}, {}
    for kind, bucket in buckets.items():
        adjacency = sparse.coo_matrix((bucket["weight"], (bucket["row"], bucket["col"])),
                                      shape=(len(frame), len(frame)), dtype=np.float32).tocsr()
        degree = np.asarray(adjacency.sum(axis=1)).ravel()
        inverse = np.divide(1., degree, out=np.zeros_like(degree), where=degree > 0)
        operators[kind] = (sparse.diags(inverse) @ adjacency).tocsr()
        details[kind] = {"n_directed_entries": int(adjacency.nnz), "n_isolated_nodes": int((degree == 0).sum())}
    return operators, {
        "n_nodes": len(frame), "n_complexes": len(records),
        "n_native_alanine_nodes": int((frame.wt == "A").sum()) if "wt" in frame else None,
        "mapping": "full_graph_residue_keys_to_frame_iloc_positions_within_each_complex",
        "normalization": "symmetric_edges_then_row_division_by_weighted_degree; isolates_zero",
        "native_weights": "footprint_weights: physical shared footprint boundary lengths",
        "spatial_weights": "unit weights on exported same-partner spatial_edges",
        "context_weights": "cross_weights: exported inter-partner atom contact counts",
        "node_population": "all prepared WT interface nodes including unmeasured and native alanine",
        "complexes": records, "operators": details,
    }


def load_graph_operators(source_dir, frame):
    root = Path(source_dir)
    registry = json.loads((root / "graphs.json").read_text())
    graphs = {}
    for complex_id in frame.complex_id.unique():
        with np.load(root / registry[str(complex_id)], allow_pickle=False) as data:
            graphs[str(complex_id)] = {key: data[key] for key in (
                "residue_keys", "footprint_edges", "footprint_weights",
                "spatial_edges", "cross_edges", "cross_weights",
            )}
    return build_graph_operators(frame, graphs)
