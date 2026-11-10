"""Serialization, ID alignment and failed-input behavior for deployable residuals."""
import joblib
import numpy as np
import pandas as pd
import pytest

from topoppi_affinity.embedding_baselines import fit_frozen_regressor, zero_surface_slots
from topoppi_affinity.residual_artifacts import predict_residual, predict_residual_affinity
from topoppi_affinity.residual_baselines import ResidualRegressor


def fitted_bundle():
    generator = np.random.default_rng(9)
    matrix = generator.normal(size=(12, 3))
    names = np.asarray(['plm.a', 'contact.b', 'surface.c'])
    matrix = zero_surface_slots(matrix, names)
    labels = matrix[:, 0] + matrix[:, 1] * .3
    model = fit_frozen_regressor(matrix, labels, np.repeat(['a', 'b', 'c', 'd'], 3),
                                family='rbf', params={'alpha': 1., 'gamma': .1}, antisymmetric=True)
    return {'model': ResidualRegressor(1.2, model), 'feature_names': names,
            'report': {'representation': 'contact', 'prior_method': 'frozen-prior',
                       'method': 'frozen-prior-Topo-contact-residual', 'units': 'kcal/mol'}}, matrix


@pytest.mark.parametrize('ids', [[f's{i}' for i in range(8)], [f'{i:03d}' for i in range(8)]])
def test_unlabeled_roundtrip_alignment_components_and_contact_mask(tmp_path, ids):
    bundle, matrix = fitted_bundle()
    ids = np.asarray(ids)
    prior = pd.DataFrame({'sample_id': ids, 'predicted': np.arange(8) / 10., 'target': 'ddg',
                          'method': 'frozen-prior', 'status': 'ok', 'error': ''})
    expected = bundle['model'].predict(prior.predicted, matrix[:8])
    # A contact artifact must mask surface slots even when supplied a surface NPZ.
    raw = matrix[:8].copy()
    raw[:, 2] = 100.
    np.savez(tmp_path / 'features.npz', sample_ids=ids[::-1], X=raw[::-1], feature_names=bundle['feature_names'])
    shuffled = prior.iloc[[4, 1, 7, 0, 6, 2, 5, 3]]
    shuffled.to_csv(tmp_path / 'prior.csv', index=False)
    joblib.dump(bundle, tmp_path / 'artifact.joblib')
    result = predict_residual_affinity(tmp_path / 'artifact.joblib', tmp_path / 'prior.csv',
                                      tmp_path / 'features.npz', tmp_path / 'prediction.csv')
    assert result.sample_id.tolist() == shuffled.sample_id.tolist()
    ordered = result.set_index('sample_id').loc[ids]
    np.testing.assert_allclose(ordered.predicted, expected, rtol=0, atol=1e-12)
    np.testing.assert_allclose(ordered.predicted, ordered.scaled_prior+ordered.residual, rtol=0, atol=1e-12)
    np.testing.assert_allclose(ordered.scaled_prior, 1.2*prior.predicted)
    assert ordered.status.eq('ok').all()


def test_failed_prior_not_filled_and_wrong_prior_or_schema_rejected(tmp_path):
    bundle, matrix = fitted_bundle()
    path = tmp_path / 'features.npz'
    np.savez(path, sample_ids=['a', 'b'], X=matrix[:2], feature_names=bundle['feature_names'])
    prior = pd.DataFrame({'sample_id': ['a', 'b'], 'predicted': [np.nan, .2], 'target': 'ddg',
                          'method': 'frozen-prior', 'status': ['error', 'ok'], 'error': ['missing residue', '']})
    result = predict_residual(bundle, path, prior)
    assert result.loc[0, 'error'] == 'missing residue'
    assert result.loc[0, ['predicted', 'scaled_prior', 'residual']].isna().all()
    with pytest.raises(ValueError, match='requires prior method'):
        predict_residual(bundle, path, prior.assign(method='binding-supervised-other-model'))
    np.savez(path, sample_ids=['a', 'b'], X=matrix[:2], feature_names=bundle['feature_names'][::-1])
    with pytest.raises(ValueError, match='Feature columns differ'):
        predict_residual(bundle, path, prior)
