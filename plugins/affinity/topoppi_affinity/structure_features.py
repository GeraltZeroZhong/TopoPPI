"""Matched frozen ESM-C/contact, pure-RDE WT context, and TopoPPI controls.

RDE embeddings come from the unsupervised rotamer-density checkpoint, not an
affinity-trained DDG network. Mutation readouts are odd conditional on the same
WT structural context; they do not assert a recomputed mutant structure.
"""
from __future__ import annotations

import numpy as np

from topoppi_affinity.embedding_baselines import zero_surface_slots
from topoppi_affinity.graph_features import FEATURE_NAMES, mutant_node_features

REPRESENTATIONS = {
    "rde_control_contact": (False, False),
    "rde_context": (True, False),
    "rde_control_surface": (False, True),
    "rde_context_surface": (True, True),
}


def align_structure_embedding(graph, bundle):
    """Map author-encoded rows to exact graph residue keys, preserving omissions."""
    lookup = {str(key): index for index, key in enumerate(bundle["residue_key"])}
    indices = [lookup[str(key)] for key in graph["residue_key"]]
    return bundle["embedding"][indices], bundle["encoded_mask"][indices].astype(bool)


def structure_pool(graph, embedding, encoded_mask, *, target, mutations=""):
    """Pool valid WT context; require every changed mutation site to be encoded."""
    side = graph["side"]
    contact = graph["x"][:, FEATURE_NAMES.index("contact.partner_distance")] <= 10.
    valid = np.asarray(encoded_mask, dtype=bool)
    coverage = {"n_nodes": len(valid), "n_encoded": int(valid.sum()),
                "encoded_fraction": float(valid.mean())}
    for partner in (0, 1):
        for scope, scope_mask in (("global", side == partner), ("contact", (side == partner) & contact)):
            coverage[f"side{partner}_{scope}_nodes"] = int(scope_mask.sum())
            coverage[f"side{partner}_{scope}_encoded"] = int((scope_mask & valid).sum())
    if target == "pkd":
        parts, names = [], []
        for scope, scope_mask in (("global", np.ones(len(valid), dtype=bool)), ("contact", contact)):
            means = []
            for partner in (0, 1):
                selected = valid & scope_mask & (side == partner)
                if not selected.any():
                    raise ValueError(f"RDE has no encoded side {partner} {scope} residues")
                means.append(np.asarray(embedding[selected], dtype=np.float64).mean(axis=0))
            values = np.concatenate([means[0] + means[1], np.abs(means[0] - means[1])])
            parts.append(values)
            names.extend(f"rde.{scope}.sum_absdiff.{index}" for index in range(len(values)))
        vector = np.concatenate(parts)
    else:
        mutant, _, mutation_mask = mutant_node_features(graph, mutations)
        missing = mutation_mask & ~valid
        if missing.any():
            keys = ", ".join(map(str, graph["residue_key"][missing]))
            raise ValueError(f"RDE missing mutation residues: {keys}")
        chemical_delta = mutant[mutation_mask, :26] - graph["x"][mutation_mask, :26]
        vector = np.einsum("ni,nj->ij", chemical_delta.astype(np.float64),
                           np.asarray(embedding[mutation_mask], dtype=np.float64)).reshape(-1)
        names = [f"rde.chemical_delta_context.{chemical}.{channel}"
                 for chemical in FEATURE_NAMES[:26] for channel in range(embedding.shape[1])]
        coverage["n_mutation_nodes"] = int(mutation_mask.sum())
        coverage["n_mutation_encoded"] = int((mutation_mask & valid).sum())
    return vector.astype(np.float32), np.asarray(names), coverage


def matched_structure_features(base, base_names, structure, structure_names):
    """Four equal-dimensional views with exactly shared non-ablated entries."""
    without_surface = zero_surface_slots(base, base_names)
    names = np.concatenate([base_names, structure_names])
    values = {}
    for representation, (use_rde, use_surface) in REPRESENTATIONS.items():
        values[representation] = np.concatenate([
            base if use_surface else without_surface,
            structure if use_rde else np.zeros_like(structure),
        ], axis=-1)
    return values, names
