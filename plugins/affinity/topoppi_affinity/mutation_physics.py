"""Atomic solvent-accessibility context for substitutions on a fixed WT structure."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from topoppi_representation.chemistry import classify_aa

from topoppi_affinity.features import parse_mutations as _parse_mutations


def parse_mutations(value):
    return [(m['wt'], m['chain'], m['position'], m['insertion'], m['mt']) for m in _parse_mutations(value)]


SASA_COLUMNS = ("sasa.isolated_partner", "sasa.bound_complex", "sasa.buried")

RESIDUE_VALUES = SASA_COLUMNS + ("rasa.isolated", "rasa.bound", "rasa.buried", "buried_fraction")

LOCAL_RADII = (4.0, 8.0)

REGIONS = ("core", "rim", "support", "noninterface", "unknown_interface")

def compute_residue_context(prepared):
    """Return one canonical selected-partner atomic table, independent of query."""
    from Bio.Data.PDBData import protein_letters_3to1_extended, residue_sasa_scales
    from topoppi.residues import source_atom_residue_labels
    from topoppi_representation.descriptors.sasa import residue_sasa

    groups = sorted(p.chains for p in prepared.partners)
    frame = residue_sasa(prepared.source.get("path", ""), *groups, prepared=prepared,
                         n_points=960, probe_radius=1.4)
    context = {}
    atom_count = 0
    for side, group in enumerate(groups):
        partner = prepared.partner(group)
        coords, atoms = partner.coords, partner.atoms
        labels = source_atom_residue_labels(atoms)
        for coord, atom, key in zip(coords, atoms, labels, strict=True):
            residue = atom.get_parent()
            identity = (side, str(key))
            if identity not in context:
                context[identity] = {
                    "chain": str(residue.get_parent().id),
                    "residue_number": int(residue.id[1]),
                    "insertion_code": str(residue.id[2]).strip(),
                    "aa": protein_letters_3to1_extended.get(residue.resname.upper(), "X"),
                    "coords": [],
                }
            context[identity]["coords"].append(coord)
            atom_count += 1
    extra = []
    for record in frame.to_dict("records"):
        item = context[(int(record["side"]), str(record["residue_key"]))]
        center = np.mean(item["coords"], axis=0)
        extra.append({k: item[k] for k in ("chain", "residue_number", "insertion_code", "aa")} |
                     {"center_x": center[0], "center_y": center[1], "center_z": center[2],
                      "n_heavy_atoms": len(item["coords"])})
    frame = pd.concat([frame.reset_index(drop=True), pd.DataFrame(extra)], axis=1)
    maximum = {protein_letters_3to1_extended[name]: value
               for name, value in residue_sasa_scales["Wilke"].items()}
    frame["max_asa_A2"] = frame.aa.map(maximum)
    frame["rasa.isolated"] = frame["sasa.isolated_partner"] / frame.max_asa_A2
    frame["rasa.bound"] = frame["sasa.bound_complex"] / frame.max_asa_A2
    frame["rasa.buried"] = frame["sasa.buried"] / frame.max_asa_A2
    isolated = frame["sasa.isolated_partner"].to_numpy(float)
    frame["buried_fraction"] = np.divide(frame["sasa.buried"], isolated,
                                          out=np.zeros(len(frame)), where=isolated > 0)
    interface = frame["sasa.buried"] > 1e-8
    frame["region"] = "noninterface"
    frame.loc[interface & frame.max_asa_A2.isna(), "region"] = "unknown_interface"
    frame.loc[interface & (frame["rasa.isolated"] < .25), "region"] = "support"
    frame.loc[interface & (frame["rasa.isolated"] >= .25) & (frame["rasa.bound"] < .25), "region"] = "core"
    frame.loc[interface & (frame["rasa.bound"] >= .25), "region"] = "rim"
    frame["nis"] = frame["rasa.bound"] >= .05
    frame["chemical_class"] = frame.aa.map(classify_aa)
    frame["n_points"] = 960
    frame["probe_radius_A"] = 1.4
    if frame.duplicated(["chain", "residue_number", "insertion_code"]).any():
        raise ValueError("WT selected partners contain duplicate chain/residue/insertion identities")
    return frame


def _reduce(values, kind):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan
    return float({"sum": np.sum, "mean": np.mean, "max": np.max}[kind](values))

def _site_context(table):
    """Per-WT-residue local means; the atomic table is the entire input."""
    values = table[list(RESIDUE_VALUES)].to_numpy(float)
    centers = table[["center_x", "center_y", "center_z"]].to_numpy(float)
    sides = table.side.to_numpy(int)
    tree = cKDTree(centers)
    site = table[list(RESIDUE_VALUES)].copy()
    for radius in LOCAL_RADII:
        neighbors = tree.query_ball_point(centers, radius)
        for domain in ("all", "opposite"):
            groups = [np.asarray([j for j in js if domain == "all" or sides[j] != sides[i]], dtype=int)
                      for i, js in enumerate(neighbors)]
            site[f"local{radius:g}A.{domain}.n_residues"] = [len(g) for g in groups]
            for column, name in enumerate(RESIDUE_VALUES):
                site[f"local{radius:g}A.{domain}.mean.{name}"] = [_reduce(values[g, column], "mean") for g in groups]
    return site

def _global_context(table):
    result = {"atomic.global.n_residues": float(len(table)),
              "atomic.global.partner_residue_count_absdiff": float(abs((table.side == 0).sum() - (table.side == 1).sum())),
              "atomic.global.unknown_max_asa_count": float(table.max_asa_A2.isna().sum())}
    for name in RESIDUE_VALUES:
        result[f"atomic.global.sum.{name}"] = _reduce(table[name], "sum")
        result[f"atomic.global.mean.{name}"] = _reduce(table[name], "mean")
        for kind in ("sum", "mean"):
            a, b = [_reduce(table.loc[table.side == side, name], kind) for side in (0, 1)]
            result[f"atomic.global.partner_{kind}_absdiff.{name}"] = abs(a - b)
    for region in REGIONS:
        result[f"atomic.global.region_fraction.{region}"] = float((table.region == region).mean())
    nis = table["nis"].astype(bool)
    for chemical in ("P", "A", "C", "X"):
        result[f"atomic.global.nis_fraction.{chemical}"] = float((table.loc[nis, "chemical_class"] == chemical).mean()) if nis.any() else np.nan
    return result

class AtomicQueryContext:
    def __init__(self, table):
        self.table = table.copy()
        self.table["insertion_code"] = self.table["insertion_code"].fillna("").astype(str)
        self.sites = _site_context(self.table)
        self.global_values = _global_context(self.table)
        self.lookup = {(str(r.chain), int(r.residue_number), str(r.insertion_code)): i
                       for i, r in enumerate(self.table.itertuples())}
        if len(self.lookup) != len(self.table):
            raise ValueError("Atomic table has duplicate mutation position identities")

    def describe(self, mutations):
        changes = parse_mutations(mutations)
        matched, missing, mismatches = [], [], []
        for wt, chain, number, insertion, _ in changes:
            index = self.lookup.get((chain, number, insertion))
            if index is None:
                missing.append(f"{chain}:{number}{insertion}")
            elif self.table.iloc[index].aa != wt:
                mismatches.append(f"{chain}:{number}{insertion}:expected_{wt}:observed_{self.table.iloc[index].aa}")
            else:
                matched.append(index)
        result = dict(self.global_values)
        result.update({"atomic.query.n_sites": float(len(changes)), "atomic.query.matched_sites": float(len(matched)),
                       "atomic.query.missing_sites": float(len(missing)), "atomic.query.wt_mismatch_sites": float(len(mismatches))})
        for kind in ("sum", "mean", "max"):
            for name in self.sites:
                result[f"atomic.query.{kind}.{name}"] = _reduce(self.sites.iloc[matched][name], kind)
        for region in REGIONS:
            result[f"atomic.query.region_count.{region}"] = float((self.table.iloc[matched].region == region).sum()) if matched else np.nan
        result["atomic.query.nis_count"] = float(self.table.iloc[matched].nis.astype(bool).sum()) if matched else np.nan
        quality = {"status": "wt_mismatch" if mismatches else "missing_sites" if missing else "ok",
                   "n_sites": len(changes), "matched_sites": len(matched), "missing": missing, "wt_mismatches": mismatches,
                   "geometry": "fixed_WT", "relaxed_mutant": False}
        return result, quality
