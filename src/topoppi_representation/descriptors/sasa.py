"""Isolated, bound, and buried solvent-accessible residue areas."""

from __future__ import annotations

import numpy as np
import pandas as pd
from Bio.PDB import Chain, Model, Residue
from Bio.PDB.SASA import ShrakeRupley

from topoppi.interface import chain_group, load_partner_data
from topoppi.io.io_loader import PDBLoader
from topoppi.residues import source_atom_residue_labels

SASA_COLUMNS = ["sasa.isolated_partner", "sasa.bound_complex", "sasa.buried"]


def residue_sasa(structure_path, chains_a, chains_b, *, prepared=None, n_points=960, probe_radius=1.4):
    """Return SASA in Å², preserving core residue keys and canonical side order.

    The isolated value removes only the other selected partner; intra-partner
    occlusion remains, including between chains in a multichain partner. Delta
    SASA is isolated minus complex SASA per residue, with no division by two.
    No relaxation, mutant modeling, solvent, or unselected chain enters the result.
    """
    groups = sorted((chain_group(chains_a), chain_group(chains_b)))
    if set(groups[0]) & set(groups[1]):
        raise ValueError("Partner chain groups overlap.")
    loader = PDBLoader(str(structure_path)) if prepared is None else None
    model = Model.Model(0)
    labels, sides = [], []
    for side, group in enumerate(groups):
        # Synthetic chains keep all atoms of one partner in a single SR entity.
        # Original chain and insertion identities remain in the lookup labels.
        chain = Chain.Chain(str(side))
        model.add(chain)
        atoms = prepared.partner(group).atoms if prepared is not None else load_partner_data(loader, group)[1]
        residue_lookup = {}
        for atom, key in zip(atoms, source_atom_residue_labels(atoms), strict=True):
            if key not in residue_lookup:
                residue = Residue.Residue((" ", len(labels), " "), atom.parent.resname, " ")
                chain.add(residue)
                residue_lookup[key] = residue
                labels.append(str(key))
                sides.append(side)
            residue_lookup[key].add(atom.copy())
    sr = ShrakeRupley(probe_radius=probe_radius, n_points=n_points)
    sr.compute(model, level="R")
    bound = np.zeros(len(labels))
    isolated = np.zeros(len(labels))
    for residue in model.get_residues():
        bound[residue.id[1]] = residue.sasa
    for chain in model:
        sr.compute(chain, level="R")
        for residue in chain:
            isolated[residue.id[1]] = residue.sasa
    buried = isolated - bound
    if np.any(buried < -1e-8):
        raise ValueError("Identical-coordinate complex SASA exceeds isolated-partner SASA.")
    return pd.DataFrame(
        {
            "residue_key": labels,
            "side": sides,
            "sasa.isolated_partner": isolated,
            "sasa.bound_complex": bound,
            "sasa.buried": np.maximum(buried, 0.0),
        }
    )
