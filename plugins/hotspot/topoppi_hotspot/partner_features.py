"""Partner-relative feature construction over the full structural interface."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from .features import FEATURE_SETS, NEIGHBOR_FEATURES, neighbor_means

KEYS = ["complex_id", "residue_key"]


QUERY_COLUMNS = ["complex_id", "side"]


SELF = tuple(FEATURE_SETS["self"])


SELF_CONTINUOUS = tuple(column for column in SELF if not column.startswith("aa."))


def _relative_columns(columns):
    return [f"relative.partner.{role}.{column}"
            for column in columns for role in ("zscore", "percentile")]


def feature_specifications():
    """Return the complete, fixed predictor lists for the four linear models."""
    relative_self = list(SELF) + _relative_columns(SELF_CONTINUOUS)
    specifications = {"absolute_self": list(SELF), "relative_self": relative_self}
    for name, kind in (("relative_spatial", "spatial"), ("relative_topoppi", "topology")):
        neighbors = [f"{kind}.neighbor_count"] + [f"{kind}.mean.{source}" for source in NEIGHBOR_FEATURES]
        specifications[name] = relative_self + neighbors + _relative_columns(neighbors)
    return specifications


def _query_ids(frame):
    return np.asarray([json.dumps([str(complex_id), str(side)], separators=(",", ":"))
                       for complex_id, side in frame[QUERY_COLUMNS].itertuples(index=False, name=None)])


def _shuffle_neighbor_sources(all_features, graphs, seed):
    """Recompute both fixed graph readouts with paired source permutations."""
    result = all_features.copy().reset_index(drop=True)
    rng = np.random.default_rng(seed)
    for complex_id, part in result.groupby("complex_id", sort=True):
        graph = graphs[str(complex_id)]
        keys = np.asarray(graph["residue_keys"]).astype(str)
        graph_frame = part.set_index("residue_key").loc[keys]
        positions = pd.Series(part.index.to_numpy(), index=part.residue_key).loc[keys].to_numpy()
        values = graph_frame[NEIGHBOR_FEATURES].to_numpy(float)
        source = values.copy()
        side = graph_frame.side.to_numpy()
        interface = graph_frame["contact.partner_distance"].to_numpy(float) <= 6.
        for partner in sorted(np.unique(side)):
            selected = np.flatnonzero((side == partner) & interface)
            selected = selected[np.argsort(keys[selected])]
            source[selected] = values[rng.permutation(selected)]
        for kind, edge_name in (("spatial", "spatial_edges"), ("topology", "native_neighbor_edges")):
            edges = np.asarray(graph[edge_name], dtype=int).reshape(-1, 2)
            if len(edges) and not ((side[edges[:, 0]] == side[edges[:, 1]]) & interface[edges].all(axis=1)).all():
                raise ValueError("Neighbor controls require same-partner interface edges.")
            means, counts = neighbor_means(source, edges)
            result.loc[positions, f"{kind}.neighbor_count"] = counts
            for number, column in enumerate(NEIGHBOR_FEATURES):
                result.loc[positions, f"{kind}.mean.{column}"] = means[:, number]
    return result


def prepare_partner_features(all_features, *, graphs=None, shuffle_seed=None):
    """Construct partner-only relative statistics before joining any labels.

    ``graphs`` maps complex IDs to dictionaries loaded from the existing NPZ
    organization graphs. It is required only for a neighbor-shuffle control.
    The returned table includes native alanine and unmeasured WT interface rows.
    """
    if all_features.duplicated(KEYS).any():
        raise ValueError("Each WT residue must occur once per complex.")
    if shuffle_seed is not None:
        if graphs is None:
            raise ValueError("Neighbor source shuffling requires native and spatial graphs.")
        all_features = _shuffle_neighbor_sources(all_features, graphs, shuffle_seed)
    frame = all_features.loc[all_features["contact.partner_distance"] <= 6.].copy()
    frame = frame.sort_values(KEYS).reset_index(drop=True)
    frame["query_id"] = _query_ids(frame)
    descriptors = list(SELF_CONTINUOUS)
    for kind in ("spatial", "topology"):
        descriptors += [f"{kind}.neighbor_count"] + [f"{kind}.mean.{source}" for source in NEIGHBOR_FEATURES]
    if not np.isfinite(frame[list(SELF) + descriptors].to_numpy(float)).all():
        raise ValueError("Partner ranking requires finite WT structural descriptors.")
    grouped = frame.groupby(QUERY_COLUMNS, sort=False)
    extra = {}
    for source in descriptors:
        values = grouped[source]
        mean, std = values.transform("mean"), values.transform("std", ddof=0)
        extra[f"context.partner.mean.{source}"] = mean
        extra[f"context.partner.deviation.{source}"] = frame[source] - mean
        extra[f"relative.partner.zscore.{source}"] = (frame[source] - mean) / std.where(std > 0, 1.)
        extra[f"relative.partner.percentile.{source}"] = (values.rank(method="average") - .5) / values.transform("size")
    return pd.concat([frame, pd.DataFrame(extra, index=frame.index)], axis=1)
