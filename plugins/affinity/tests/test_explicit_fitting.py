"""Explicit fitting inputs preserve alignment and prediction contracts."""
import importlib.util
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from topoppi_affinity.external_mutation import read_mutation_inputs
from topoppi_affinity.residual_artifacts import fit_residual_model, predict_residual


def test_mutation_adapter_reads_all_declared_rows_without_split_columns(tmp_path):
    frame = pd.DataFrame({'sample_id': ['a', 'b', 'c'], 'complex_id': ['two', 'one', 'two']})
    path = tmp_path / 'inputs.csv'
    frame.to_csv(path, index=False)
    assert read_mutation_inputs(path).sample_id.tolist() == ['a', 'b', 'c']
    assert read_mutation_inputs(path, limit_complexes=1).sample_id.tolist() == ['b']


def test_absolute_fit_tool_aligns_two_independently_ordered_feature_files(tmp_path, monkeypatch):
    tool = Path(__file__).parents[1] / 'tools/models/fit_absolute_model.py'
    spec = importlib.util.spec_from_file_location('fit_absolute_tool', tool)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    frame = pd.DataFrame({'sample_id': ['a', 'b', 'c', 'd'],
                          'complex_id': ['a', 'b', 'c', 'd'],
                          'family_group': ['x', 'x', 'y', 'y'], 'pkd': [6., 8., 7., 9.]})
    table = tmp_path / 'training.csv'
    frame.to_csv(table, index=False)
    sequence = np.array([[0., 1.], [1., 0.], [2., 0.], [0., 2.]])
    descriptors = np.array([[1.], [2.], [3.], [4.]])
    plm, phys, artifact = (tmp_path / p for p in ['sequence.npz', 'descriptors.npz', 'model.joblib'])
    np.savez(plm, sample_ids=frame.sample_id.to_numpy(str)[::-1], X=sequence[::-1],
             feature_names=['sequence.a', 'sequence.b'])
    np.savez(phys, sample_ids=frame.sample_id.to_numpy(str), X=descriptors,
             feature_names=['contact.a'], feature_groups=['contact'])
    monkeypatch.setattr('sys.argv', ['fit', '--training-table', str(table),
        '--sequence-features', str(plm), '--descriptor-features', str(phys), '--output', str(artifact)])
    module.main()
    bundle = joblib.load(artifact)
    from topoppi_affinity.legacy_kernel import fit_legacy_bundle
    expected = fit_legacy_bundle(sequence, ['sequence.a', 'sequence.b'], descriptors,
        ['contact.a'], ['contact'], frame)
    matrix = np.column_stack([sequence, descriptors])
    np.testing.assert_allclose(bundle['model'].predict(matrix), expected['model'].predict(matrix), atol=1e-12)


def test_residual_fit_accepts_explicit_table_and_reordered_features(tmp_path):
    rng = np.random.default_rng(5)
    count = 8
    x = rng.normal(size=(count, 2))
    frame = pd.DataFrame({'sample_id': [f's{i}' for i in range(count)],
        'complex_id': [f'c{i // 2}' for i in range(count)],
        'family_group': [f'f{i // 2}' for i in range(count)],
        'mutations': [f'AA{i + 1}V' for i in range(count)], 'ddg': x[:, 0] + .2 * x[:, 1]})
    prior = pd.DataFrame({'sample_id': frame.sample_id, 'predicted': x[:, 0],
        'target': 'ddg', 'method': 'synthetic_prior', 'status': 'ok', 'error': ''})
    features, prior_path, artifact = (tmp_path / p for p in ['features.npz', 'prior.csv', 'model.joblib'])
    np.savez(features, X=x[::-1], sample_ids=frame.sample_id.to_numpy(str)[::-1], feature_names=['a', 'b'])
    prior.to_csv(prior_path, index=False)
    report = fit_residual_model(frame, features, prior_path, artifact,
                                task='mutation_family', representation='surface')
    result = predict_residual(artifact, features, prior)
    assert report['n_train'] == count
    assert report['train_sample_ids'] == frame.sample_id.tolist()
    assert result.sample_id.tolist() == frame.sample_id.tolist()
    assert np.isfinite(result.predicted).all()
    np.testing.assert_allclose(result.predicted, result.scaled_prior + result.residual)
