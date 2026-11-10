"""Prevent wrong-interface/cache reuse and mixing absolute/mutation references."""
import math
import pickle

import numpy as np
import pandas as pd
import pytest

from topoppi_affinity.author_references import read_author_partners
from topoppi_affinity.cached_author_references import cache_coverage, import_gearbind_table
from topoppi_affinity.proaffinity_inputs import (
    align_observed_residues,
    atom_histograms,
    interaction_pairs,
    read_fasta_chains,
    write_proaffinity_pdb,
)


def test_cache_identity_uses_partners_and_unordered_mutations_not_row_numbers(tmp_path):
    lmdb = pytest.importorskip('lmdb')
    entries = [dict(id=91, complex='1ABC', group_ligand=['A'], group_receptor=['B', 'C'],
                    mutstr='AA2V,GC5S,AA2V', pdbcode='91_1ABC', num_muts=2, ddG='unreadable')]
    with (tmp_path/'entries.pkl').open('wb') as stream:
        pickle.dump(entries, stream)
    lmdb.open(str(tmp_path/'structures.lmdb'), subdir=False, map_size=1024*1024).close()
    pd.DataFrame([dict(sample_id='matching', pdb_id='1abc', chains_a='C,B', chains_b='A',
                       mutations='GC5S,AA2V', ddg='do not read'),
                  dict(sample_id='wrong-partners', pdb_id='1abc', chains_a='A', chains_b='B',
                       mutations='AA2V,GC5S', ddg='do not read')]).to_csv(tmp_path/'manifest.csv', index=False)
    result = cache_coverage(tmp_path/'manifest.csv', tmp_path, tmp_path/'coverage.csv')
    assert result.status.tolist() == ['ok', 'error']
    assert result.iloc[0].cache_entry == '91_1ABC'
    assert len(result) == 2


def test_published_predictions_keep_author_scope_and_never_match_pdb_only(tmp_path):
    source = pd.DataFrame({'#Pdb': ['1ABC_A_B', '1ABC_A_B', '1ABC_A_C'],
                            'Mutation(s)_cleaned': ['AA2V']*3,
                            'GearBind': [1., 3., 20.], 'fold': [0, 0, 3]})
    source.to_csv(tmp_path/'published.csv', index=False)
    pd.DataFrame([dict(sample_id='a', pdb_id='1abc', chains_a='B', chains_b='A', mutations='AA2V'),
                  dict(sample_id='b', pdb_id='1abc', chains_a='A', chains_b='D', mutations='AA2V')]).to_csv(
                      tmp_path/'manifest.csv', index=False)
    result = import_gearbind_table(tmp_path/'manifest.csv', tmp_path/'published.csv', tmp_path/'out.csv')
    assert result.predicted.iloc[0] == 2.0
    assert result.status.tolist() == ['ok', 'error']
    assert not result.local_model_executed.any()
    assert result.comparison_scope.eq('published_author_split_reference').all()


def test_complete_fasta_mapping_auth_chains_and_case_sensitive_graph_keys(tmp_path):
    sequences = read_fasta_chains('>1abc_1|Chains X [auth A], Y [auth b]|protein\nACDE\n')
    assert sequences == {'A': 'ACDE', 'b': 'ACDE'}
    residues = [{'type': c} for c in 'ACDE']
    align_observed_residues(residues, sequences['A'])
    assert [r['number'] for r in residues] == [0, 1, 2, 3]
    index = tmp_path/'partners.txt'
    index.write_text('1abc_A_b\tA; b;\n')
    assert read_author_partners(index)['1abc_A_b'] == (frozenset('A'), frozenset('b'))


def test_proaffinity_preserves_pdb_element_columns_and_selects_first_model_partners(tmp_path):
    # Two-letter elements in real benchmark structures fail in ADFR when the
    # locally installed BioPython formatter drops a blank before column 73.
    calcium = 'HETATM    1 CA    CA A   1       1.000   2.000   3.000  1.00 20.00          CA  '
    carbon = 'ATOM      2  CA  ALA B   2       4.000   5.000   6.000  1.00 20.00           C  '
    excluded = carbon[:21] + 'C' + carbon[22:]
    source, target = tmp_path/'original.pdb', tmp_path/'partners.pdb'
    source.write_text('MODEL        1\n' + '\n'.join([calcium, carbon, excluded])
                      + '\nENDMDL\nMODEL        2\n' + carbon + '\nENDMDL\nEND\n')
    write_proaffinity_pdb(dict(structure_path=source, chains_a='A', chains_b='B'), target)
    assert target.read_text().splitlines() == [calcium, carbon, 'END']
    assert target.read_text().splitlines()[0][76:78] == 'CA'


def test_contact_and_histogram_match_original_distance_bins_at_cutoffs():
    def residue(number, xs, types):
        return {'number': number, 'atoms': [dict(type='CA' if i == 0 else 'N',
                x=str(x), y='0', z='0', pdbqt_type=t) for i, (x, t) in enumerate(zip(xs, types, strict=True))]}
    a = residue(0, [0., 1.5], ['C', 'N'])
    b = residue(1, [3.5, 4.], ['C', 'N'])
    far = residue(2, [8.], ['C'])
    pairs = interaction_pairs([a, b, far], [a, b, far], 3.5, intra=True)
    assert [(x['number'], y['number']) for x, y, _ in pairs] == [(0, 1), (1, 0)]
    classes = ['C_C', 'C_N', 'N_N']
    expected = np.zeros(30)
    for first in a['atoms']:
        for second in b['atoms']:
            distance = abs(float(first['x'])-float(second['x']))
            category = '_'.join(sorted([first['pdbqt_type'], second['pdbqt_type']]))
            index = (min(math.ceil(distance/(3.5/10)), 10)-1)*3+classes.index(category)
            expected[index] += 1
    np.testing.assert_array_equal(atom_histograms([(a, b, 3.5)], classes, 3.5).numpy()[0], expected)




