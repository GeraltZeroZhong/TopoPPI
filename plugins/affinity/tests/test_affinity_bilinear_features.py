import numpy as np
import torch
from torch_geometric.data import Data

from topoppi_affinity.bilinear_features import bilinear_features, symmetric_bilinear
from topoppi_affinity.embedding_baselines import zero_surface_slots


def test_bilinear_original_definition_and_partner_exchange():
    rng = np.random.default_rng(2)
    x, edges, weights = rng.normal(size=(8, 7)), np.array([[0, 1, 2], [4, 5, 6]]), np.array([1., 2., 3.])
    raw = sum(w * np.outer(x[i], x[j]) for (i, j), w in zip(edges.T, weights, strict=True))
    raw = (raw + raw.T) / 2
    raw = np.sign(raw) * np.sqrt(abs(raw))
    raw /= np.linalg.norm(raw)
    i, j = np.triu_indices(7)
    expected = raw[i, j] * np.where(i == j, 1., np.sqrt(2.))
    np.testing.assert_allclose(symmetric_bilinear(x, edges, weights), expected, atol=1e-12)
    np.testing.assert_allclose(symmetric_bilinear(x, edges[::-1], weights), expected, atol=1e-12)


def test_mutation_pair_symmetry_and_no_surface_control():
    torch.manual_seed(5)
    x, esm = torch.rand(6, 45), torch.rand(6, 8)
    data = Data(x=x, x_mut=x.clone(), esm=esm, esm_mut=esm.clone(),
                mutation_mask=torch.tensor([True, False, False, False, False, False]),
                contact_mask=torch.ones(6, dtype=torch.bool), side=torch.tensor([0, 0, 0, 1, 1, 1]),
                edge_index=torch.tensor([[0, 1, 2, 3, 4, 5], [3, 4, 5, 0, 1, 2]]),
                edge_attr=torch.ones(6, 8))
    projection = np.random.default_rng(11).normal(size=(8, 3))
    noop, _, _ = bilinear_features(data, projection, target='ddg')
    np.testing.assert_array_equal(noop, np.zeros_like(noop))
    data.x_mut[0, :26] += .4
    data.esm_mut[0] += .1
    forward, names, _ = bilinear_features(data, projection, target='ddg')
    reverse = data.clone()
    reverse.x, reverse.x_mut = data.x_mut, data.x
    reverse.esm, reverse.esm_mut = data.esm_mut, data.esm
    backward, _, _ = bilinear_features(reverse, projection, target='ddg')
    np.testing.assert_allclose(forward, -backward, atol=1e-7)
    altered = data.clone()
    altered.x[:, 30:] *= 19
    altered.x_mut[:, 30:] *= 19
    changed, _, _ = bilinear_features(altered, projection, target='ddg')
    np.testing.assert_array_equal(zero_surface_slots(forward, names), zero_surface_slots(changed, names))
