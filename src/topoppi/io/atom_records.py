"""Portable author residue and atom identities for saved TopoPPI views."""
import numpy as np
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue


def atom_records(atoms):
    records = []
    for atom in atoms or ():
        residue = atom.get_parent()
        records.append({
            "name": atom.get_name(), "fullname": atom.get_fullname(), "element": atom.element,
            "serial": atom.get_serial_number(), "bfactor": atom.get_bfactor(),
            "occupancy": atom.get_occupancy(), "altloc": atom.get_altloc(),
            "chain": str(residue.get_parent().id), "residue_id": list(residue.get_id()),
            "resname": residue.get_resname(),
        })
    return records


def restore_atoms(records, coords):
    if len(records) != len(coords):
        raise ValueError("Atom records and coordinates have different lengths.")
    chains, residues, atoms = {}, {}, []
    for record, xyz in zip(records, coords, strict=True):
        chain_id = record["chain"]
        if chain_id not in chains:
            chains[chain_id] = Chain(chain_id)
        residue_id = tuple(record["residue_id"])
        key = (chain_id, residue_id)
        if key not in residues:
            residues[key] = Residue(residue_id, record["resname"], "")
            chains[chain_id].add(residues[key])
        atom = Atom(record["name"], np.array(xyz, dtype=float), record["bfactor"], record["occupancy"],
                    record["altloc"], record["fullname"], record["serial"], element=record["element"])
        residues[key].add(atom)
        atoms.append(atom)
    return atoms


