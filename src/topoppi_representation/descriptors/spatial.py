"""Global surface organization and mutation-conditioned spatial descriptors."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..chemistry import HYDROPATHY
from .surface import local_queries, physical_targets

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWYX"


def global_spatial_descriptors(surface):
    """Return 60 partner-symmetric descriptors of a decoded surface.

    Areas use transported source mass, the radius is an area-weighted second
    moment in Å², and chemical proximity uses the common three-dimensional
    frame. The returned names define the column order.
    """
    targets = physical_targets(surface)
    weights = surface.weights
    one_hot = np.array([[float(letter == aa) for aa in AMINO_ACIDS] for letter in surface.amino_acids])
    node_aa = np.asarray(surface.mass @ one_hot)
    parts = []
    local_names = ["mass", "radius_squared", "L_H"] + ["AA_fraction_" + aa for aa in AMINO_ACIDS]
    for side in (0, 1):
        take = np.asarray(surface.side) == side
        mass = float(weights[take].sum())
        if mass:
            center = np.average(surface.xyz[take], axis=0, weights=weights[take])
            radius = float(np.average(np.square(surface.xyz[take] - center).sum(axis=1), weights=weights[take]))
            composition = node_aa[take].sum(axis=0) / mass
        else:
            radius, composition = 0.0, np.zeros(len(AMINO_ACIDS))
        parts.append(np.r_[mass, radius, targets["O_by_side"][str(side)]["L_H"], composition])
    values = np.r_[parts[0] + parts[1], np.abs(parts[0] - parts[1]), targets["X"], np.log1p(targets["W"])]
    names = (["partner_sum." + n for n in local_names] + ["partner_absdiff." + n for n in local_names]
             + ["proximity." + n for n in targets["X_states"]] + ["proximity.log1p_weight"])
    return values, names


@dataclass
class SpatialDescriptorCache:
    """Share a decoded WT surface across substitution queries.

    ``describe`` accepts a sequence of ``(wt, chain, number, insertion, mutant)``
    tuples, including multicharacter author chain IDs. It returns values, names
    and footprint coverage. ``None`` requests the 60 global descriptors; a
    substitution query appends 61 fields. Local kernels use graph distances
    at 4 and 12 Å. Missing footprints produce NaN local summaries.
    """
    surface: object
    local: dict = field(default_factory=dict)

    def __post_init__(self):
        self.global_values, self.global_names = global_spatial_descriptors(self.surface)
        self.residue_lookup = {}
        for i, key in enumerate(self.surface.residue_keys):
            fields = str(key).split(":")
            self.residue_lookup[fields[-3], fields[-1]] = i

    def describe(self, mutations=None):
        values, names = self.global_values.copy(), list(self.global_names)
        if mutations is None:
            return values, names, {"task_context": "absolute_binding"}
        changes = list(mutations)
        if isinstance(mutations, str):
            raise TypeError("Pass parsed (wt, chain, number, insertion, mutant) tuples to describe().")
        wt_counts, mt_counts = np.zeros(len(AMINO_ACIDS)), np.zeros(len(AMINO_ACIDS))
        context, covered, missing = [], [], []
        for wt, chain, number, insertion, mt in changes:
            wt_counts[AMINO_ACIDS.index(wt)] += 1
            mt_counts[AMINO_ACIDS.index(mt)] += 1
            idx = self.residue_lookup.get((chain, f"{number}{insertion}"))
            if idx is None:
                context.append([np.nan] * 5)
                missing.append(f"{chain}:{number}{insertion}")
                continue
            if self.surface.amino_acids[idx] != wt:
                raise ValueError(f"WT residue mismatch for {chain}:{number}{insertion}")
            key = self.surface.residue_keys[idx]
            if key not in self.local:
                self.local[key] = local_queries(self.surface, [key], (4.0, 12.0))["queries"][key]
            row = self.local[key]
            if not all(row[str(scale)]["covered"] for scale in (4.0, 12.0)):
                context.append([np.nan] * 5)
                missing.append(key)
                continue
            covered.append(key)
            context.append([row["4.0"]["support_mass"], row["4.0"]["value"], row["12.0"]["value"],
                            row["4.0"]["kernel_mass"], row["12.0"]["kernel_mass"]])
        context_names = ["footprint_mass", "hydropathy_4A", "hydropathy_12A", "kernel_mass_4A", "kernel_mass_12A"]
        observed = np.asarray(context, dtype=float).reshape(-1, 5)
        summaries = []
        for reduction in (np.mean, np.min, np.max):
            summaries.extend([float(reduction(col[np.isfinite(col)])) if np.isfinite(col).any() else np.nan
                              for col in observed.T])
        extra = np.r_[len(changes), len(covered), len(missing), wt_counts, mt_counts,
                      sum(HYDROPATHY[mt] - HYDROPATHY[wt] for wt, _, _, _, mt in changes), summaries]
        extra_names = (["query.n_mutations", "query.covered_sites", "query.missing_sites"]
                       + ["query.wt_" + aa for aa in AMINO_ACIDS] + ["query.mt_" + aa for aa in AMINO_ACIDS]
                       + ["query.delta_hydropathy"]
                       + [f"query.{stat}.{name}" for stat in ("mean", "min", "max") for name in context_names])
        return np.r_[values, extra], names + extra_names, {
            "task_context": "fixed_WT_geometry_plus_mutation_identity", "missing_sites": missing,
            "covered_sites": covered, "relaxed_mutant_structure": False}
