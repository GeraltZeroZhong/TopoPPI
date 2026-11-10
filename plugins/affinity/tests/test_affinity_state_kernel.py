"""Check linear-observation regression and physical prediction algebra."""
import numpy as np
from scipy.sparse import csr_matrix
from sklearn.kernel_ridge import KernelRidge

from topoppi_affinity.state_kernel import fit_state_kernel


def test_absolute_observations_equal_standard_kernel_ridge():
    x = np.random.default_rng(4).normal(size=(8, 3))
    labels = np.linspace(5., 9., 8)
    weights = np.linspace(.5, 1.5, 8)
    model = fit_state_kernel(x, np.eye(8), labels, weights, alpha=.1, gamma_factor=1.)
    standard = KernelRidge(alpha=.1, kernel='rbf', gamma=1./3.).fit(
        model.scaler.transform(x), labels - np.average(labels, weights=weights), sample_weight=weights)
    np.testing.assert_allclose(model.predict(x), standard.predict(model.scaler.transform(x)) + model.mean, atol=1e-12)


def test_mixed_state_constraints_reverse_identity_and_cycle():
    x = np.array([[0.], [1.], [2.]])
    operator = csr_matrix([[1., 0., 0.], [1., -1., 0.], [0., 1., -1.]])
    labels = np.array([8., 1., -2.])
    model = fit_state_kernel(x, operator, labels, np.ones(3), alpha=1e-8, gamma_factor=1.)
    np.testing.assert_allclose(operator @ model.predict(x), labels, atol=1e-6)
    ab, bc, ac = [model.difference(x[[i]], x[[j]], 298.15) for i, j in [(0, 1), (1, 2), (0, 2)]]
    np.testing.assert_allclose(ab + bc, ac, atol=1e-12)
    np.testing.assert_allclose(ab, -model.difference(x[[1]], x[[0]], 298.15))
    np.testing.assert_array_equal(model.difference(x, x, 298.15), np.zeros(3))
