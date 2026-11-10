"""Real protocol hazards: row identities, targets, state mapping and heldout use."""
import json

import numpy as np
import pandas as pd
import pytest

from topoppi_affinity.benchmark_models import LoadedFeatures
from topoppi_affinity.cath_protocol import prepare_cath_protocol
from topoppi_affinity.ppb_protocol import match_cached_structures, prepare_ppb_protocol, prepare_ppb_sequences
from topoppi_affinity.representation_controls import nested_features


def test_cath_raw_row_mapping_targets_and_official_val_is_test(tmp_path):
    raw = pd.DataFrame({'#Pdb': ['1AAA_A_B', '2BBB_A_B', '3CCC_A_B'],
                        'Mutation(s)_cleaned': ['AA1C', 'CA2D', 'DA3E']})
    raw.to_csv(tmp_path/'raw.csv', sep=';', index=False)
    released = pd.DataFrame({'#Pdb': ['2_3CCC', '0_1AAA', '1_2BBB'],
                            'complex_PPI': ['3CCC_A_B', '1AAA_A_B', '2BBB_A_B'],
                            'Mutation(s)_cleaned': ['DA3E', 'AA1C', 'CA2D'],
                            'cath_fold': ['', 'train', 'val'], 'ddG': [3.1, 1.1, 2.1],
                            'ddG_FoldX': [3., 1., 'missing'], 'Flex-ddG': [3., 1., 2.]})
    released.to_csv(tmp_path/'author.csv', index=False)
    manifest = pd.DataFrame({'sample_id': ['skempi2:00001', 'skempi2:00002', 'skempi2:00003'],
                             'complex_id': ['a', 'b', 'c'], 'ddg': [1., 2., 3.],
                             'variant_group': ['v1', 'v2', 'v3'], 'interaction_group': ['i1', 'i2', 'i3'],
                             'family_group': ['f1', 'f2', 'f3'],
                             'sequence_a': ['AA', 'CC', 'DD'], 'sequence_b': ['EE', 'FF', 'GG']})
    manifest.to_csv(tmp_path/'manifest.csv', index=False)
    report = prepare_cath_protocol(tmp_path/'manifest.csv', tmp_path/'raw.csv', tmp_path/'author.csv', tmp_path/'out')
    mapped = pd.read_csv(tmp_path/'out/manifest.csv')
    assert mapped.partition.tolist() == ['train', 'test', 'unassigned']
    np.testing.assert_allclose(mapped.ddg, [1.1, 2.1, 3.1])
    np.testing.assert_allclose(mapped.ddg_at_reported_temperature, [1., 2., 3.])
    assert report['target_mean_absolute_difference_from_reported_temperature'] == pytest.approx(0.1)
    assert report['models_trained'] is False
    ref = pd.read_csv(tmp_path/'out/reference_predictions.csv')
    assert ref.status.eq('error').sum() == 1
    # The released prefix is positional, so a different raw release must fail.
    raw.loc[0, 'Mutation(s)_cleaned'] = 'AA1D'
    raw.to_csv(tmp_path/'raw.csv', sep=';', index=False)
    with pytest.raises(ValueError, match='row numbering'):
        prepare_cath_protocol(tmp_path/'manifest.csv', tmp_path/'raw.csv', tmp_path/'author.csv', tmp_path/'wrong')


def test_ppb_preserves_already_mutated_sequences_labels_and_reports_overlap(tmp_path):
    raw = pd.DataFrame([{'Source Data Set': 'SKEMPI', 'Complex ID': f'c{i}', 'PDB': '1AAA',
        'Ligand Chains': 'A', 'Receptor Chains': 'B', 'Ligand Sequences': 'AYD',
        'Receptor Sequences': seq, 'Mutations': 'CA2Y', 'KD(M)': 1e-8,
        'split': split, 'Affinity PubMed ID': 'test', 'Temperature(K)': '298'}
        for i, (seq, split) in enumerate([('FGH', 'train'), ('KKL', 'val'), ('KKL', 'test')])])
    raw.to_csv(tmp_path/'ppb.csv', index=False)
    report = prepare_ppb_protocol(tmp_path/'ppb.csv', tmp_path/'out')
    manifest = pd.read_csv(tmp_path/'out/manifest.csv')
    assert manifest.sequence_a.tolist() == ['AYD']*3
    assert manifest.complex_id.nunique() == 1
    assert manifest.source_complex_id.tolist() == ['c0', 'c1', 'c2']
    np.testing.assert_allclose(manifest.pkd, 8.)
    assert report['exact_overlap']['train/test']['identical_individual_chain_sequences'] == 1
    assert report['exact_overlap']['validation/test']['exact_sequence_pairs'] == 1
    assert manifest.exact_chain_seen_in_official_train.all()
    register = prepare_ppb_sequences(tmp_path/'out/manifest.csv', tmp_path/'seq')
    assert register['unique_sequences'] == 3
    assert 'AYD' in set(pd.read_csv(tmp_path/'seq/sequences.csv').sequence)


def test_ppb_cached_geometry_requires_correct_chain_specific_state(tmp_path):
    (tmp_path/'manifests').mkdir()
    original = pd.DataFrame([dict(sample_id='old', pdb_id='1aaa', chains_a='A,B', chains_b='C',
                                  sequence_a='AAA|CCC', sequence_b='DDD')])
    original.to_csv(tmp_path/'manifests/absolute.csv', index=False)
    frame = pd.concat([original.assign(sample_id='match'), original.assign(sample_id='wrong', sequence_a='CCC|AAA'),
                       original.assign(sample_id='swapped', chains_a='C', chains_b='B,A', sequence_a='DDD', sequence_b='CCC|AAA')])
    matched = match_cached_structures(frame, tmp_path).set_index('sample_id')
    assert matched.loc['match', 'status'] == matched.loc['swapped', 'status'] == 'ok'
    assert matched.loc['wrong', 'status'] == 'excluded'


def test_nested_surface_ablation_preserves_every_other_slot():
    ids = np.array(['a', 'b'])
    seq = LoadedFeatures(ids, np.arange(6).reshape(2, 3), np.array(['sequence.a', 'sequence.b', 'sequence.c']))
    native = LoadedFeatures(ids, np.array([[2., 3.], [4., 5.]]), np.array(['plm.contact', 'surface.x']))
    a, b = nested_features(seq, native, view='surface'), nested_features(seq, native, view='contact')
    np.testing.assert_array_equal(a.X[:, :-1], b.X[:, :-1])
    np.testing.assert_array_equal(b.X[:, -1], 0.)
    pure = nested_features(seq, native, view='sequence')
    np.testing.assert_array_equal(pure.X[:, :3], seq.X)
    np.testing.assert_array_equal(pure.X[:, 3:], 0.)




@pytest.mark.parametrize("limit_complexes,expected", [(None, ['a', 'b', 'c']), (1, ['a', 'b'])])
def test_evoef_uses_supplied_manifest_rows(tmp_path, monkeypatch, limit_complexes, expected):
    from topoppi_affinity.external_energy import run_evoef
    frame = pd.DataFrame({'sample_id': ['a', 'b', 'c'], 'complex_id': ['c1', 'c1', 'c2']})
    frame.to_csv(tmp_path/'manifest.csv', index=False)
    def pretend_complex(binary, output, group):
        return [{'sample_id': sid, 'status': 'ok', 'predicted': .5} for sid in group.sample_id]
    monkeypatch.setattr('topoppi_affinity.external_energy._run_complex', pretend_complex)
    result = run_evoef(tmp_path/'manifest.csv', tmp_path/'binary', tmp_path/'out',
                       workers=1, limit_complexes=limit_complexes)
    assert result.sample_id.tolist() == expected
    metadata = json.loads((tmp_path/'out/evoef.json').read_text())
    assert metadata['cohort'] == 'all_manifest_rows' and metadata['n_attempted'] == len(expected)
    assert metadata['limit_complexes'] == limit_complexes


@pytest.mark.parametrize("states,summary", [
    (["possible", "possible"], "possible"),
    (["known", "possible"], "known"),
    (["excluded", "excluded"], "excluded"),
    (["excluded", "unknown"], "unknown"),
    (["known", "unknown"], "known"),
])
def test_external_mutation_saved_overlap_reflects_input_rows(tmp_path, states, summary):
    from topoppi_affinity.external_mutation import _save_predictions

    records = [dict(sample_id=f"s{i}", training_overlap=state,
                    status="error" if state == "unknown" else "ok",
                    predicted=np.nan if state == "unknown" else .5)
               for i, state in enumerate(states)]
    result = _save_predictions(records, tmp_path, "predictions", {"method": "external"})
    metadata = json.loads((tmp_path / "predictions.json").read_text())
    saved = pd.read_csv(tmp_path / "predictions.csv")
    assert result.training_overlap.tolist() == saved.training_overlap.tolist() == states
    assert metadata["training_overlap"] == summary
    assert metadata["training_overlap_counts"] == {state: states.count(state) for state in set(states)}
    assert metadata["n_attempted"] == len(states)
    assert metadata["n_success"] == len(states) - states.count("unknown")

