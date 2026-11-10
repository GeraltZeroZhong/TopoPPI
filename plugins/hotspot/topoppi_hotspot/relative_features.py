"""Label-free within-interface feature decomposition for residue ranking.

The normalization population is every WT residue in the contact interface,
including unmeasured sites. It is never the subset with energetic labels. All
statistics are computed separately for each structure's known interface, so a
held-out interface needs no labels or training-interface population statistics.

Centering alone does not change a linear model's within-interface order at fixed
slopes: the subtracted interface mean is a constant score offset. Percentiles and
interface/partner scales introduce additional relative information. Explicit
interface means are retained for explanation, not included as model predictors.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .features import FEATURE_SETS

SELF_DESCRIPTORS = [
    "surface.area", "surface.footprint_area", "surface.footprint_fraction",
    "contact.degree", "contact.atom_pairs", "contact.neighbor_count",
    "contact.betweenness", "contact.pagerank", "chemical.aromatic", "chemical.volume",
]


def _neighbor_descriptors(kind):
    return [f"{kind}.neighbor_count"] + [
        f"{kind}.mean.{source}" for source in (
            "chemical.aromatic", "contact.degree", "contact.atom_pairs",
            "surface.footprint_area",
        )
    ]


def add_relative_features(all_features, *, contact_cutoff=6.):
    """Return interface residues, model feature sets, and JSON-ready provenance.

    ``all_features`` is the all-residue WT structural table, with ``complex_id``
    identifying one interface and ``side`` identifying its binding partner. Keep
    all WT contact residues when calling this function, then join energetic
    labels afterwards. Rows outside the contact cutoff are omitted before any
    statistics are computed. Input columns and row order are otherwise retained.

    For every selected descriptor, x = interface mean + residue deviation. The
    population standard deviation uses ddof=0; constant strata have z-score 0.
    Percentiles use average tie rank, (rank - 0.5) / n, making singleton and
    entirely constant strata 0.5. Partner strata are (complex_id, side).

    The returned feature sets retain their original absolute predictors and add
    residue deviations, interface z-scores/percentiles and partner z-scores/
    percentiles. Interface mean columns are explanatory context only. No label,
    mutation, family, or interface identifier is a predictor.
    """
    frame = all_features.loc[
        all_features["contact.partner_distance"] <= contact_cutoff
    ].copy()
    descriptors = SELF_DESCRIPTORS + [
        column for kind in ("spatial", "topology", "shuffled")
        for column in _neighbor_descriptors(kind)
    ]
    if not np.isfinite(frame[descriptors].to_numpy(float)).all():
        raise ValueError("Relative features require finite WT structural descriptors.")

    interface = frame.groupby("complex_id", sort=False)
    partner = frame.groupby(["complex_id", "side"], sort=False)
    extra, feature_map = {}, {}
    for source in descriptors:
        mapping = {"absolute": source,
                   "context": f"context.interface.mean.{source}",
                   "deviation": f"relative.interface.deviation.{source}",
                   "interface_zscore": f"relative.interface.zscore.{source}",
                   "interface_percentile": f"relative.interface.percentile.{source}",
                   "partner_zscore": f"relative.partner.zscore.{source}",
                   "partner_percentile": f"relative.partner.percentile.{source}"}
        mean = interface[source].transform("mean")
        deviation = frame[source] - mean
        extra[mapping["context"]] = mean
        extra[mapping["deviation"]] = deviation
        for name, group in (("interface", interface), ("partner", partner)):
            grouped = group[source]
            group_mean = grouped.transform("mean")
            std = grouped.transform("std", ddof=0)
            extra[mapping[f"{name}_zscore"]] = (
                (frame[source] - group_mean) / std.where(std > 0, 1.)
            )
            extra[mapping[f"{name}_percentile"]] = (
                (grouped.rank(method="average") - .5) / grouped.transform("size")
            )
        feature_map[source] = mapping
    frame = pd.concat([frame, pd.DataFrame(extra, index=frame.index)], axis=1)

    feature_sets = {name: list(columns) for name, columns in FEATURE_SETS.items()}
    for kind in ("self", "spatial", "topology", "shuffled"):
        sources = SELF_DESCRIPTORS + ([] if kind == "self" else _neighbor_descriptors(kind))
        feature_sets[f"relative_{kind}"] = list(FEATURE_SETS[kind]) + [
            column for source in sources for role, column in feature_map[source].items()
            if role not in ("absolute", "context")
        ]
    populations = [
        {"complex_id": str(complex_id), "n_interface_residues": int(len(part)),
         "n_by_partner": {str(side): int(count) for side, count in part.groupby("side").size().items()}}
        for complex_id, part in frame.groupby("complex_id", sort=True)
    ]
    metadata = {
        "protocol": "hotspot_relative_v2",
        "normalization_population": "all WT contact-interface residues, including unmeasured sites",
        "contact_cutoff_angstrom": float(contact_cutoff),
        "n_input_residues": int(len(all_features)),
        "n_interface_residues": int(len(frame)),
        "n_interfaces": int(frame.complex_id.nunique()),
        "interface_group": ["complex_id"],
        "partner_group": ["complex_id", "side"],
        "zscore": "population standard deviation (ddof=0); zero variance gives 0",
        "percentile": "(average tie rank - 0.5) / population size",
        "decomposition": "absolute = interface mean + residue deviation; no rounding",
        "context_in_model": False,
        "labels_used": False,
        "centering_limit": "At fixed linear slopes, interface centering alone cannot change within-interface ranking.",
        "feature_columns": feature_map,
        "interface_populations": populations,
    }
    return frame, feature_sets, metadata
