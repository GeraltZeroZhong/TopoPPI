import numpy as np

from topoppi_affinity.group_dro_baselines import (
    FixedTrainingCoordinates,
    adversarial_update,
    family_risk_matrix,
    training_path,
)


def test_family_risk_balances_complexes_not_variant_counts():
    names, matrix = family_risk_matrix(['a', 'a', 'a', 'b', 'c'], ['f', 'f', 'f', 'f', 'g'])
    assert names.tolist() == ['f', 'g']
    np.testing.assert_allclose(matrix @ np.array([2., 2., 2., 8., 10.]), [5., 10.])
    np.testing.assert_allclose(matrix.sum(axis=1), 1.)


def test_adversarial_update_matches_author_equation():
    q, risks = np.array([.2, .8]), np.array([4., 1.])
    expected = q * np.exp(.01 * risks)
    np.testing.assert_allclose(adversarial_update(q, risks, .01), expected / expected.sum())
    assert adversarial_update(q, risks, .01)[0] > q[0]


def test_training_path_is_odd_and_reweights_training_risk():
    x = np.array([[1., 2.], [2., 1.], [3., -1.], [4., 2.]])
    y = np.array([1., 2., 8., -3.])
    groups, families = ['a', 'a', 'b', 'c'], ['f', 'f', 'g', 'g']
    coords = FixedTrainingCoordinates(x, y, groups, True)
    path = list(training_path(coords, groups, families, alpha=1000., steps=2))
    np.testing.assert_allclose(path[1][3], adversarial_update(path[0][3], path[0][4], .01))
    for _, model, _, _, _ in path:
        np.testing.assert_allclose(model.predict(-x), -model.predict(x), atol=1e-12)
        np.testing.assert_allclose(model.predict(np.zeros((1, 2))), 0., atol=1e-12)
