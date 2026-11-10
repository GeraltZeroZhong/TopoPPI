"""Chemical-contact, solvent-accessibility and local-context residue descriptors.

The feature scheme uses complete 6 Å WT interfaces, 10 Å same-partner spatial
neighbors, and raw/partner-relative values. It preserves the published ranking
model's amino-acid classes and normalization population.
"""
from __future__ import annotations

from itertools import combinations_with_replacement

import numpy as np
import pandas as pd

from .features import BASE_FEATURES, CHEMICAL, neighbor_means
from .solvent_accessibility import SASA_COLUMNS

# Tien et al. (2013), Table 1 theoretical allowed maxima, Å².
# https://doi.org/10.1371/journal.pone.0080635
MAX_ASA = dict(zip("ARNDCEQGHILKMFPSTWYV", [129., 274., 195., 193., 167., 223., 225., 104., 224., 197.,
                  201., 236., 224., 240., 159., 155., 172., 285., 263., 174.], strict=True))
CHEMISTRY = {"charged": "DEKR", "polar": "CHNQSTY", "apolar": "AFGILMPVW"}
SASA_EXTRA = ["sasa.buried_fraction", "sasa.rasa_isolated", "sasa.rasa_bound", "sasa.rasa_buried",
              "sasa.region_core", "sasa.region_rim", "sasa.region_support", "sasa.region_noninterface"]
CONTACT_TYPES = ["typed_contact." + "_".join(pair)
                 for pair in combinations_with_replacement(sorted(CHEMISTRY), 2)]
NEIGHBOR_SOURCES = CHEMICAL + ["contact.degree", "contact.atom_pairs"] + SASA_COLUMNS + SASA_EXTRA
SCHEMA = "hotspot_physical_spatial_165_v1"


def _expanded(columns):
    return list(columns) + [f"relative.partner.{role}.{column}" for column in columns
                           if not column.startswith("aa.") for role in ("zscore", "percentile")]


def feature_columns():
    """Ordered predictor names for the physical-context ranking model."""
    return (_expanded(BASE_FEATURES + SASA_COLUMNS + SASA_EXTRA + CONTACT_TYPES)
            + _expanded(["spatial_baseline.neighbor_count"]
                        + [f"spatial_baseline.mean.{column}" for column in NEIGHBOR_SOURCES]))


def _relative(frame, columns):
    grouped, extra = frame.groupby(["complex_id", "side"], sort=False), {}
    for column in columns:
        values = grouped[column]
        mean, std = values.transform("mean"), values.transform("std", ddof=0)
        extra[f"relative.partner.zscore.{column}"] = (frame[column] - mean) / std.where(std > 0, 1.)
        extra[f"relative.partner.percentile.{column}"] = (values.rank(method="average") - .5) / values.transform("size")
    return pd.concat([frame.drop(columns=list(extra), errors="ignore"), pd.DataFrame(extra, index=frame.index)], axis=1)


def _mapped_edges(part, graph, name):
    lookup = pd.Series(np.arange(len(part)), index=part.residue_key)
    mapping = lookup.reindex(np.asarray(graph["residue_keys"]).astype(str)).to_numpy()
    edges = mapping[np.asarray(graph[name], int).reshape(-1, 2)]
    if not np.isfinite(edges).all():
        raise ValueError("Physical-context edges reference a residue outside the full WT interface.")
    return edges.astype(int)


def add_physical_context(frame, graphs):
    """Add physical features to all WT interface rows, including native alanine.

    ``frame`` includes per-residue isolated, bound and buried SASA. ``graphs``
    contains the corresponding cross-partner contact and spatial edges. Relative
    statistics use all interface residues of each partner before ranking filters.
    """
    frame = frame.copy().reset_index(drop=True)
    maxima = frame.wt.map(MAX_ASA)
    if maxima.isna().any():
        raise ValueError("The physical-context model requires standard amino-acid identities for interface residues.")
    isolated, bound, buried = (frame[column].to_numpy(float) for column in SASA_COLUMNS)
    frame["sasa.buried_fraction"] = np.divide(buried, isolated, out=np.zeros_like(buried), where=isolated > 0)
    for name, values in (("isolated", isolated), ("bound", bound), ("buried", buried)):
        frame[f"sasa.rasa_{name}"] = values / maxima
    interface = buried > 0.
    support = interface & (frame["sasa.rasa_isolated"] < .25)
    rim = interface & (frame["sasa.rasa_bound"] >= .25)
    for name, mask in (("support", support), ("rim", rim), ("core", interface & ~support & ~rim),
                       ("noninterface", ~interface)):
        frame[f"sasa.region_{name}"] = np.asarray(mask, float)
    columns = CONTACT_TYPES + ["spatial_baseline.neighbor_count"] + [f"spatial_baseline.mean.{c}" for c in NEIGHBOR_SOURCES]
    extra = pd.DataFrame(0., index=frame.index, columns=columns)
    classes = {aa: group for group, aas in CHEMISTRY.items() for aa in aas}
    for complex_id, positions in frame.groupby("complex_id", sort=True).indices.items():
        part, graph = frame.iloc[positions], graphs[complex_id]
        chemistry, side = part.wt.map(classes).to_numpy(), part.side.to_numpy()
        cross = _mapped_edges(part, graph, "cross_edges")
        if len(cross) and np.any(side[cross[:, 0]] == side[cross[:, 1]]):
            raise ValueError("Chemical-contact edges must cross partners.")
        for i, j in cross:
            column = "typed_contact." + "_".join(sorted((chemistry[i], chemistry[j])))
            extra.loc[positions[[i, j]], column] += 1.
        edges = _mapped_edges(part, graph, "spatial_edges")
        if len(edges) and np.any(side[edges[:, 0]] != side[edges[:, 1]]):
            raise ValueError("Spatial context must stay within each partner.")
        means, counts = neighbor_means(part[NEIGHBOR_SOURCES].to_numpy(float), edges)
        extra.loc[positions, "spatial_baseline.neighbor_count"] = counts
        extra.loc[positions, [f"spatial_baseline.mean.{c}" for c in NEIGHBOR_SOURCES]] = means
    if not np.array_equal(extra[CONTACT_TYPES].sum(axis=1).to_numpy(), frame["contact.degree"].to_numpy()):
        raise ValueError("Chemical contact counts disagree with the interface contact degree.")
    return _relative(pd.concat([frame, extra], axis=1), SASA_COLUMNS + SASA_EXTRA + columns)
