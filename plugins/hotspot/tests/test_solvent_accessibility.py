"""Independent synthetic checks for reusable numerical primitives."""
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBIO, PDBParser
from topoppi_hotspot.solvent_accessibility import SASA_COLUMNS, residue_sasa

FIXTURE = Path(__file__).resolve().parent / "fixtures/tiny_complex.pdb"

def test_sasa_burial_conservation_and_partner_symmetry():
    frame = residue_sasa(FIXTURE, "A", "B")
    assert frame.residue_key.tolist() == ["A:GLY:1", "A:ALA:2", "B:GLY:1", "B:ALA:2"]
    assert frame["sasa.buried"].sum() > 0
    np.testing.assert_allclose(frame["sasa.isolated_partner"] - frame["sasa.bound_complex"], frame["sasa.buried"])
    assert (frame[SASA_COLUMNS] >= 0).all().all()
    pd.testing.assert_frame_equal(frame, residue_sasa(FIXTURE, "B", "A"))

def test_separating_partner_eliminates_burial_without_changing_isolated_sasa(tmp_path):
    structure = PDBParser(QUIET=True).get_structure("tiny", str(FIXTURE))
    before = residue_sasa(FIXTURE, "A", "B")
    for atom in structure[0]["B"].get_atoms():
        atom.coord += np.array([100., 0., 0.])
    writer = PDBIO()
    writer.set_structure(structure)
    path = tmp_path / "separated.pdb"
    writer.save(str(path))
    after = residue_sasa(path, "A", "B")
    np.testing.assert_allclose(after["sasa.buried"], 0., atol=1e-8)
    np.testing.assert_allclose(after["sasa.isolated_partner"], before["sasa.isolated_partner"])

def test_multichain_partner_keeps_internal_occlusion(tmp_path):
    structure = PDBParser(QUIET=True).get_structure("tiny", str(FIXTURE))
    original_pair = residue_sasa(FIXTURE, "A", "B")
    chain_c = structure[0]["A"].copy()
    chain_c.id = "C"
    for atom in chain_c.get_atoms():
        atom.coord += np.array([100., 0., 0.])
    structure[0].add(chain_c)
    writer = PDBIO()
    writer.set_structure(structure)
    path = tmp_path / "three_chains.pdb"
    writer.save(str(path))
    result = residue_sasa(path, "A,B", "C").set_index("residue_key")
    first = result.loc[original_pair.residue_key]
    # A and B are now one partner: their isolated SASA is the old bound SASA.
    np.testing.assert_allclose(first["sasa.isolated_partner"], original_pair["sasa.bound_complex"])
    np.testing.assert_allclose(result["sasa.buried"], 0., atol=1e-8)
