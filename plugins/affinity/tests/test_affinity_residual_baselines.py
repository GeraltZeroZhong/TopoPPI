"""Coverage and training-boundary checks for the stability residual experiment."""
import numpy as np
import pandas as pd
import pytest

from topoppi_affinity import residual_baselines as residual


def test_prior_alignment_excludes_failed_and_nonfinite_without_zero_filling():
    frame = pd.DataFrame({'sample_id': ['c', 'a', 'b', 'd'], 'complex_id': ['p'] * 4})
    predictions = pd.DataFrame({
        'sample_id': ['a', 'b', 'c'], 'predicted': [1., 0., np.nan], 'target': ['ddg'] * 3,
        'status': ['ok', 'error', 'ok'], 'error': ['', 'unsupported residue', ''],
    })
    matched, prior, usable, coverage = residual.align_prior(frame, predictions)
    assert matched.sample_id.tolist() == ['a']
    np.testing.assert_array_equal(prior, [1.])
    np.testing.assert_array_equal(usable, [False, True, False, False])
    assert coverage.loc[0, 'error'] == 'Nonfinite prior prediction'
    assert coverage.loc[3, 'prior_status'] == 'missing'


def test_origin_slope_uses_equal_complex_weights_and_preserves_zero():
    slope = residual.fit_origin_slope([1., 1., 1.], [1., 3., 10.], ['a', 'a', 'b'])
    assert slope == pytest.approx(6.)
    model = residual.ResidualRegressor(slope)
    np.testing.assert_allclose(model.predict([0., 2., -2.]), [0., 12., -12.])


def test_inner_split_keeps_family_and_repeated_mutants_together():
    frame = pd.DataFrame({
        'complex_id': ['a', 'a', 'b', 'c', 'd', 'e'],
        'mutations': ['AA1V', 'AA1V', 'BB2C', 'CC3D', 'DD4F', 'EE5G'],
        'family_group': ['f1', 'f1', 'f1', 'f2', 'f3', 'f4'],
        'initial_interaction_group': ['i1', 'i1', 'i2', 'i3', 'i4', 'i5'],
    })
    for task, group in [('mutation_family', frame.family_group),
                        ('mutation_within', frame.complex_id + '|' + frame.mutations)]:
        for train, validation in residual.inner_partitions(frame, task):
            assert set(group[train]).isdisjoint(group[validation])
            assert np.all(train ^ validation)
