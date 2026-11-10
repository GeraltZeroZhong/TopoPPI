import json

import numpy as np
import pytest
from Bio.PDB import MMCIFIO, PDBIO, PDBParser
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from topoppi_affinity.external_energy import parse_binding_energy
from topoppi_affinity.external_inputs import external_partner_groups, prepare_numbered_complex, remap_mutations


def test_numbering_handles_insertions_and_skips_calcium_ions(tmp_path):
    structure, model = Structure('example'), Model(0)
    structure.add(model)
    for chain_id, residues in [('A', [('ALA', (' ', 10, ' ')), ('GLY', (' ', 10, 'A')), ('CA', ('H_CA', 900, ' '))]),
                                ('B', [('SER', (' ', 20, ' '))])]:
        chain = Chain(chain_id)
        model.add(chain)
        for index, (name, identity) in enumerate(residues):
            residue = Residue(identity, name, '')
            atom = Atom('CA', np.array([index, 0., 0.]), 0., 1., ' ', ' CA ', index + 1, element='CA' if name == 'CA' else 'C')
            residue.add(atom)
            chain.add(residue)
    original = tmp_path / 'original.pdb'
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(original))
    row = {'structure_path': str(original), 'chains_a': 'A', 'chains_b': 'B', 'sequence_a': 'AG', 'sequence_b': 'S'}
    mapping = prepare_numbered_complex(row, tmp_path / 'mapped.pdb')
    assert len(mapping) == 3
    assert remap_mutations('GA10AV,SB20A', mapping) == 'GA2V,SB1A'
    with pytest.raises(ValueError, match='WT identity'):
        remap_mutations('AA10AV', mapping)


def test_binding_energy_parser_requires_grouped_total():
    assert parse_binding_energy('Some term = 2\nTotal = -12.32\n') == -12.32
    with pytest.raises(ValueError):
        parse_binding_energy('Total = 1\nTotal = 2\n')
    with pytest.raises(ValueError):
        parse_binding_energy('Total = nan\n')


def test_mmcif_author_chains_use_consistent_pdb_aliases_and_insertion_mapping(tmp_path):
    structure, model = Structure('author_ids'), Model(0)
    structure.add(model)
    for number, (chain_id, name, resid, insertion) in enumerate(
            [('A', 'ALA', 9, ' '), ('alpha', 'GLY', 10, 'A'), ('1', 'SER', 20, ' ')]):
        chain = Chain(chain_id)
        residue = Residue((' ', resid, insertion), name, '')
        residue.add(Atom('CA', np.array([float(number), 0., 0.]), 0., 1., ' ', ' CA ', number + 1, element='C'))
        chain.add(residue)
        model.add(chain)
    original = tmp_path / 'author.cif'
    writer = MMCIFIO()
    writer.set_structure(structure)
    writer.save(str(original))
    mapped = tmp_path / 'mapped.pdb'
    mapping = prepare_numbered_complex({'structure_path': original, 'chains_a': 'alpha,1', 'chains_b': 'A'}, mapped)
    assert external_partner_groups('alpha,1', 'A') == ('B1', 'A')
    assert remap_mutations('alpha:G10AV,1:S20A', mapping) == 'S11A,GB1V'
    assert set(PDBParser(QUIET=True).get_structure('mapped', mapped)[0].child_dict) == {'A', 'B', '1'}
    saved = json.loads(mapped.with_suffix('.numbering.json').read_text())
    alpha = next(row for row in saved if row['chain'] == 'alpha')
    assert alpha == {'chain': 'alpha', 'external_chain': 'B', 'source_number': 10,
                     'source_insertion': 'A', 'number': 1, 'aa': 'G'}
