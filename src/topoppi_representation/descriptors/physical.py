"""The 477-field physicochemical dictionary computed from complete partners.

SASA: Shrake--Rupley, 1.4 Å probe; relative SASA: Tien et al. (2013)
theoretical maxima (Biopython's Wilke table). Chemical contact classes follow
Vangone and Bonvin (2015). Both partner buried-area losses are summed.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from Bio.Data.PDBData import protein_letters_3to1_extended, residue_sasa_scales
from scipy.spatial import cKDTree

from topoppi.residues import source_atom_residue_labels

from .sasa import residue_sasa

SCHEMA = "physicochemical_477_v1"
AA = "ACDEFGHIKLMNPQRSTVWYX"
CLASSES = {"P": "CHNQSTYW", "A": "AFGIVLMP", "C": "DEKR", "X": "X"}
MAX_ASA = {protein_letters_3to1_extended[name]: value for name, value in residue_sasa_scales["Wilke"].items()}
SASA_FIELDS = ("sasa.isolated_partner", "sasa.bound_complex", "sasa.buried")


def chemistry(aa):
    return next((key for key, members in CLASSES.items() if aa in members), "X")


def summarize(frame, pairs):
    side_values = []
    for side in (0, 1):
        part = frame[frame.side.eq(side)]
        masks = {
            "all": np.ones(len(part), bool),
            "interface": part["sasa.buried"].gt(1e-8),
            **{reg: part.region.eq(reg) for reg in ("core", "rim", "support")},
            "nis": part.nis,
        }
        values = {}
        for scope, mask in masks.items():
            selected = part[mask]
            values[f"{scope}.n_residues"] = len(selected)
            for col in (
                "sasa.isolated_partner",
                "sasa.bound_complex",
                "sasa.buried",
                "rasa.isolated",
                "rasa.bound",
                "buried_fraction",
                "contact.degree_5.5A",
            ):
                x = selected[col].to_numpy(float)
                x = x[np.isfinite(x)]
                for stat in ("mean", "sum", "std"):
                    values[f"{scope}.{col}.{stat}"] = float(getattr(np, stat)(x)) if len(x) else 0.0
            for chem in CLASSES:
                values[f"{scope}.fraction_{chem}"] = (
                    float(selected.chemical_class.eq(chem).mean()) if len(selected) else 0.0
                )
            if scope in ("all", "interface", "nis"):
                for aa in AA:
                    values[f"{scope}.aa_{aa}"] = float(selected.aa.eq(aa).mean()) if len(selected) else 0.0
        side_values.append(values)
    result = {}
    for key in side_values[0]:
        a, b = side_values[0][key], side_values[1][key]
        result["physical." + key + ".sum"] = a + b
        result["physical." + key + ".absdiff"] = abs(a - b)
    result["physical.contact_pairs.total"] = float(len(pairs))
    for pair in ("AA", "AC", "AP", "CC", "CP", "PP", "AX", "CX", "PX", "XX"):
        result["physical.contact_pairs." + pair] = float(pairs.chemical_pair.eq(pair).sum())
    for chem in CLASSES:
        nis = frame[frame.nis]
        result["physical.nis_combined.fraction_" + chem] = (
            float(nis.chemical_class.eq(chem).mean()) if len(nis) else 0.0
        )
    return result


def annotate_residues(frame):
    """Add chemical classes, relative accessibility and interface regions."""
    frame = frame.copy()
    frame["chemical_class"] = frame.aa.map(chemistry)
    maximum = frame.aa.map(MAX_ASA).to_numpy(float)
    frame["rasa.isolated"] = frame[SASA_FIELDS[0]].to_numpy(float) / maximum
    frame["rasa.bound"] = frame[SASA_FIELDS[1]].to_numpy(float) / maximum
    frame["buried_fraction"] = frame[SASA_FIELDS[2]] / frame[SASA_FIELDS[0]].clip(lower=1e-8)
    frame["region"] = "noninterface"
    interface = frame[SASA_FIELDS[2]].gt(1e-8)
    frame.loc[interface & ~np.isfinite(frame["rasa.isolated"]), "region"] = "unknown_interface"
    frame.loc[interface & frame["rasa.isolated"].lt(0.25), "region"] = "support"
    frame.loc[interface & frame["rasa.isolated"].ge(0.25) & frame["rasa.bound"].lt(0.25), "region"] = "core"
    frame.loc[interface & frame["rasa.bound"].ge(0.25), "region"] = "rim"
    frame["nis"] = frame["rasa.bound"].ge(0.05)
    return frame


def residue_physics(interface, *, n_points=960, probe_radius=1.4):
    """Return complete-partner residue and 5.5 Å contact-pair tables.

    ``interface`` is a Core Interface or GeometryRepresentation. Partner sides
    follow sorted chain groups, preserving the dictionary's canonical ordering.
    Selected heavy atoms are reused in memory. SASA has units of Å².
    """
    groups = sorted(tuple(partner.chains) for partner in interface.partners)
    frame = residue_sasa(None, *groups, prepared=interface, n_points=n_points, probe_radius=probe_radius)
    selected = [interface.partner(group) for group in groups]
    labels = [source_atom_residue_labels(partner.atoms) for partner in selected]
    amino_acids = {}
    for partner, keys in zip(selected, labels, strict=False):
        for atom, key in zip(partner.atoms, keys, strict=False):
            amino_acids[str(key)] = protein_letters_3to1_extended.get(atom.parent.resname, "X")
    frame["aa"] = frame.residue_key.map(amino_acids)
    frame = annotate_residues(frame)
    pairs, counts = set(), {}
    other = cKDTree(selected[1].coords)
    for i, neighbors in enumerate(other.query_ball_point(selected[0].coords, 5.5)):
        for j in neighbors:
            a, b = str(labels[0][i]), str(labels[1][j])
            pairs.add((a, b))
            counts[a] = counts.get(a, 0) + 1
            counts[b] = counts.get(b, 0) + 1
    degree = {}
    for a, b in pairs:
        degree[a] = degree.get(a, 0) + 1
        degree[b] = degree.get(b, 0) + 1
    frame["contact.degree_5.5A"] = frame.residue_key.map(degree).fillna(0.0)
    frame["contact.atom_pairs_5.5A"] = frame.residue_key.map(counts).fillna(0.0)
    frame["n_points"] = int(n_points)
    contacts = pd.DataFrame(sorted(pairs), columns=["residue_a", "residue_b"])
    contacts["chemical_pair"] = [
        "".join(sorted((chemistry(amino_acids[a]), chemistry(amino_acids[b]))))
        for a, b in contacts.itertuples(index=False, name=None)
    ]
    return frame, contacts


def chemical_area_descriptors(frame):
    """Partner-symmetric chemically typed SASA and buried-area summaries."""
    values = {}
    for chemical in CLASSES:
        for field in SASA_FIELDS:
            sides = [
                float(frame.loc[(frame.side == side) & (frame.chemical_class == chemical), field].sum())
                for side in (0, 1)
            ]
            values[f"physical.chemical_area.{chemical}.{field}.sum"] = sum(sides)
            values[f"physical.chemical_area.{chemical}.{field}.absdiff"] = abs(sides[0] - sides[1])
    return values


def physical_descriptors(interface, *, n_points=960, probe_radius=1.4):
    """Compute the named 477-field dictionary from selected Core atoms."""
    residues, contacts = residue_physics(interface, n_points=n_points, probe_radius=probe_radius)
    return {**summarize(residues, contacts), **chemical_area_descriptors(residues)}


def recover_from_sasa(frame):
    """Recover the 430 fields determined by side, amino acid and three SASAs.

    Required columns are side, aa, sasa.isolated_partner, sasa.bound_complex,
    and sasa.buried. Atomic-contact descriptors require an additional input.
    """
    frame = annotate_residues(frame[["side", "aa", *SASA_FIELDS]])
    frame["contact.degree_5.5A"] = np.nan
    values = summarize(frame, pd.DataFrame(columns=["chemical_pair"]))
    values = {
        name: value for name, value in values.items() if ".contact." not in name and ".contact_pairs." not in name
    }
    return {**values, **chemical_area_descriptors(frame)}


def recover_from_carrier(surface):
    """Recover SASA-derived fields from the decoded, explicitly carried table."""
    missing = set(SASA_FIELDS) - surface.residue_fields.keys()
    if missing:
        raise ValueError(f"Required residue fields are absent: {', '.join(sorted(missing))}")
    frame = pd.DataFrame(
        {
            "side": [int(key.split(":", 1)[0][4:]) for key in surface.residue_keys],
            "aa": surface.amino_acids,
            **{field: surface.residue_fields[field] for field in SASA_FIELDS},
        }
    )
    return recover_from_sasa(frame)
