"""Tests for failures that would corrupt labels or leak experimental variants."""
from pathlib import Path

import pandas as pd
import pytest

from topoppi_affinity.data import (
    attach_structure,
    canonical_mutations,
    complex_key,
    exact_positive,
    parse_kd,
    parse_mutations,
)
from topoppi_affinity.splits import audit_splits, group_split, within_complex_split


def test_exact_kd_units_reject_bounds_and_other_endpoints():
    assert parse_kd('Kd=22.5nM') == pytest.approx(22.5e-9)
    assert parse_kd('Kd=1.2e-3uM') == pytest.approx(1.2e-9)
    for value in ['Ki=1nM', 'IC50=1nM', 'Kd<1nM', 'Kd~1nM', 'Kd=0M']:
        with pytest.raises(ValueError):
            parse_kd(value)
    for value in ['<1e-9', '298(assumed)', 'nan', '-1']:
        with pytest.raises(ValueError):
            exact_positive(value)


def test_cleaned_mutations_preserve_chain_position_and_insertion_code():
    assert parse_mutations('LI38G,AH12AV') == [('A', 'H', 12, 'A', 'V'), ('L', 'I', 38, '', 'G')]
    assert canonical_mutations('WI2F,AH12V') == 'AH12V,WI2F'
    assert complex_key('1JRH', ['L', 'H'], ['I']) == complex_key('1jrh', ['I'], ['H', 'L'])
    with pytest.raises(ValueError, match='same position'):
        parse_mutations('LI38G,LI38A')


def _write_pdb(path: Path):
    lines = []
    for serial, resname, chain, number, x in [(1, 'LEU', 'I', 38, 0), (2, 'ALA', 'H', 12, 3), (3, 'TYR', 'L', 7, 4)]:
        lines.append(f'ATOM  {serial:5d}  CA  {resname:3s} {chain}{number:4d}    {x:8.3f}{0:8.3f}{0:8.3f}{1:6.2f}{20:6.2f}           C  \n')
    path.write_text(''.join(lines) + 'END\n')


def test_multichain_template_identity_and_cleaned_numbering(tmp_path):
    path = tmp_path / '1ABC.pdb'
    _write_pdb(path)
    row = {'structure_path': str(path), 'chains_a': 'H,L', 'chains_b': 'I'}
    result = attach_structure(row.copy(), 'LI38G')
    assert result['sequence_a'] == 'A|Y'
    assert result['sequence_b'] == 'L'
    with pytest.raises(ValueError, match='wildtype_mismatch'):
        attach_structure(row.copy(), 'AI38G')
    with pytest.raises(ValueError, match='residue_missing'):
        attach_structure(row.copy(), 'LI45G')


def test_within_complex_split_keeps_repeated_variant_measurements_together():
    frame = pd.DataFrame({'complex_id': ['a'] * 12 + ['b'] * 3,
                          'mutations': [f'AA{i}V' for i in range(6)] * 2 + ['AB1V', 'AB2V', 'AB3V']})
    frame['split'] = within_complex_split(frame)
    assert frame.groupby(['complex_id', 'mutations'])['split'].nunique().max() == 1
    assert set(frame[frame.complex_id == 'a']['split']) == {'train', 'validation', 'test'}
    assert set(frame[frame.complex_id == 'b']['split']) == {'train'}


def test_group_assignment_stable_and_audit_detects_same_pdb_leakage():
    assert group_split(['a', 'a', 'b', 'c', 'd']) == group_split(['d', 'a', 'c', 'b', 'a'])
    frame = pd.DataFrame({'pdb_id': ['1abc', '1abc'], 'sequence_group': ['g', 'g'],
                          'split_complex': ['train', 'test'], 'split_sequence': ['train', 'train'],
                          'n_mutations': [0, 0]})
    with pytest.raises(ValueError, match='groups occur'):
        audit_splits(frame)


def test_exact_interaction_groups_join_swapped_alternate_pdb_pairs():
    from topoppi_affinity.splits import exact_interaction_groups
    frame = pd.DataFrame({'pdb_id': ['1aaa', '2aaa', '2aaa', '3aaa'],
                          'sequence_a': ['AAA|BBB', 'CCC', 'XXX', 'AAA'],
                          'sequence_b': ['CCC', 'BBB|AAA', 'YYY', 'ZZZ']})
    groups = exact_interaction_groups(frame)
    assert groups.iloc[0] == groups.iloc[1] == groups.iloc[2]
    assert groups.iloc[3] != groups.iloc[0]
