import numpy as np

from topoppi_affinity.geometric_scattering import SCALES, lazy_walk, scattering_signals, wavelet_signals


def test_sparse_wavelets_match_author_dense_definition():
    adjacency = np.array([[0., 2., 1.], [2., 0., 0.], [1., 0., 0.]])
    index = np.array(np.nonzero(adjacency))
    walk = lazy_walk(index, adjacency[tuple(index)], 3)
    expected = .5 * (np.eye(3) + adjacency / adjacency.sum(0))
    np.testing.assert_allclose(walk.toarray(), expected)
    signal = np.arange(6).reshape(3, 2)
    actual = wavelet_signals(walk, signal)
    for scale, value in zip(SCALES, actual, strict=True):
        reference = (np.linalg.matrix_power(expected, scale) - np.linalg.matrix_power(expected, 2 * scale)) @ signal
        np.testing.assert_allclose(value, reference, atol=1e-12)
    permutation = np.array([2, 0, 1])
    before = scattering_signals(walk, signal)
    after = scattering_signals(walk[permutation][:, permutation], signal[permutation])
    np.testing.assert_allclose(after, before[permutation], atol=1e-12)
    assert before.shape == (3, 22)
