from __future__ import annotations

import pandas as pd
import pytest

from topoppi_affinity.benchmark_data import (
    TASKS,
    add_exposure,
    attach_global_groups,
    balanced_group_folds,
    canonical_state_pair,
    make_splits,
    merge_variant_interactions,
    physical_variant_groups,
)


def write_pdb(path, chains):
    names = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'V': 'VAL'}
    lines = []
    for chain, (start, sequence) in chains.items():
        for offset, aa in enumerate(sequence):
            serial = len(lines) + 1
            lines.append(f'ATOM  {serial:5d}  CA  {names[aa]:3s} {chain}{start + offset:4d}    '
                         f'{serial * 3.8:8.3f}{1.0:8.3f}{2.0:8.3f}  1.00 20.00           C  \n')
    path.write_text(''.join(lines) + 'END\n')


def test_physical_variant_groups_use_actual_mapping_and_reverse_states(tmp_path):
    first, second, reverse = [tmp_path / f'{name}.pdb' for name in ['first', 'second', 'reverse']]
    write_pdb(first, {'A': (1, 'AC'), 'B': (1, 'DE')})
    write_pdb(second, {'X': (50, 'DE'), 'Z': (101, 'AC')})
    write_pdb(reverse, {'R': (6, 'VC'), 'S': (8, 'DE')})
    frame = pd.DataFrame([
        dict(sample_id='one', structure_path=str(first), chains_a='A', chains_b='B',
             sequence_a='AC', sequence_b='DE', mutations='AA1V', ddg=2),
        dict(sample_id='two', structure_path=str(second), chains_a='X', chains_b='Z',
             sequence_a='DE', sequence_b='AC', mutations='AZ101V', ddg=9),
        dict(sample_id='reverse', structure_path=str(reverse), chains_a='R', chains_b='S',
             sequence_a='VC', sequence_b='DE', mutations='VR6A', ddg=-2),
    ])
    result, identities = physical_variant_groups(frame)
    assert result.variant_group.nunique() == 1
    assert len(identities) == 1
    assert result.mutant_sequence_a.tolist() == ['VC', 'DE', 'AC']
    assert result.ddg.tolist() == [2, 9, -2]
    bad = frame.iloc[[0]].copy()
    bad['mutations'] = 'CA1V'
    with pytest.raises(ValueError, match='WT identity'):
        physical_variant_groups(bad)


def test_canonical_variants_preserve_partner_membership_and_chain_multiplicity():
    wt = {'A': 'AC', 'B': 'DE', 'C': 'VA'}
    mt = dict(wt, A='VC')
    original = canonical_state_pair(wt, mt, [['A', 'B'], ['C']])
    assert original == canonical_state_pair(wt, mt, [['C'], ['B', 'A']])
    assert original != canonical_state_pair(wt, mt, [['A'], ['B', 'C']])
    assert original != canonical_state_pair(wt, mt, [['A', 'A', 'B'], ['C']])


def artificial_frame():
    rows = []
    for number in range(6):
        common = dict(pdb_id=f'p{number}', complex_id=f'p{number}:A:B',
                      family_group=f'family_{number}', interaction_group=f'interaction_{number}',
                      initial_interaction_group=f'interaction_{number}')
        rows.append(dict(common, sample_id=f'a{number}', source='absolute_source',
                         benchmark_target='absolute', n_mutations=0, variant_group=''))
        for variant in range(1 if number == 0 else 2):
            rows.append(dict(common, sample_id=f'm{number}_{variant}', source='mutation_source',
                             benchmark_target='mutation', n_mutations=variant + 1,
                             variant_group=f'v{number}_{variant}'))
    return pd.DataFrame(rows)


def test_joint_folds_exposure_and_known_complex_flag_are_label_blind():
    frame = artificial_frame()
    used = frame.loc[frame.sample_id.isin(['a1', 'm2_0']), ['sample_id', 'source', 'complex_id']]
    marked = add_exposure(frame, {'used.csv': ('confirmation', used)})
    assert marked.ever_used_row.sum() == 2
    assert not marked.loc[marked.sample_id.eq('m1_0'), 'fresh_family'].item()
    assert marked.loc[marked.sample_id.eq('m1_0'), 'fresh_row'].item()
    assert marked.loc[marked.benchmark_target.eq('absolute'), 'fresh_variant'].isna().all()
    full, splits = make_splits(marked, n_folds=5, seed=9)
    assert set(splits) == set(TASKS)
    for column, group in [('fold_family', 'family_group'), ('fold_interaction', 'interaction_group')]:
        assert full.groupby(group)[column].nunique().eq(1).all()
    variants = splits['mutation_variant']
    assert not variants.loc[variants.sample_id.eq('m0_0'), 'complex_seen_in_training'].item()
    assert not variants.loc[variants.sample_id.eq('m0_0'), 'evaluation_eligible'].item()
    assert not variants.loc[variants.sample_id.eq('m2_0'), 'fresh_only'].item()
    shuffled = marked.sample(frac=1, random_state=3).copy()
    shuffled['ddg'] = range(len(shuffled))
    again = balanced_group_folds(shuffled, 'family_group', n_folds=5, seed=9)
    expected = full.groupby('family_group').fold_family.first().to_dict()
    assert again == expected


def test_global_group_lookup_retains_construct_specific_sequences():
    groups = pd.DataFrame([
        dict(complex_id='x:A:B', pdb_id='x', sequence_a='AC', sequence_b='DE',
             family_group='f1', initial_interaction_group='i1'),
        dict(complex_id='x:A:B', pdb_id='x', sequence_a='ACC', sequence_b='DEE',
             family_group='f1', initial_interaction_group='i1'),
    ])
    frame = groups.drop(columns=['family_group', 'initial_interaction_group']).copy()
    frame['sample_id'] = ['a', 'm']
    assert len(attach_global_groups(frame, groups)) == 2
    frame.loc[0, 'sequence_a'] = 'AV'
    with pytest.raises(ValueError, match='missing from global groups'):
        attach_global_groups(frame, groups)


def test_variant_bridges_merge_absolute_and_mutation_interactions():
    frame = artificial_frame()
    # The second mutation system contains the inverse of a physical variant in
    # the first; its absolute rows must follow the same merged component.
    frame.loc[frame.sample_id.eq('m1_0'), 'variant_group'] = 'v0_0'
    frame.loc[frame.family_group.eq('family_1'), 'family_group'] = 'family_0'
    merged, evidence = merge_variant_interactions(frame, n_folds=5, seed=9)
    linked = merged[merged.complex_id.isin(['p0:A:B', 'p1:A:B'])]
    assert linked.interaction_group.nunique() == 1
    assert linked.initial_interaction_group.nunique() == 2
    assert len(evidence) == 2
    assert evidence.variant_group.nunique() == 1
    assert merged.family_group.equals(frame.family_group)
