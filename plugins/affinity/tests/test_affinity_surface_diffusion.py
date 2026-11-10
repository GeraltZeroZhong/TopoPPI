"""Mapping, native author-layer gradients and exchange invariance on small surfaces."""

import numpy as np
import pytest
import torch

from topoppi_affinity.surface_diffusion import (
    pool_quadrature,
    quadrature_mapping,
)


def surface_fixture():
    """Two spherical point surfaces, with repeated whole/interface quadratures."""
    count = 80
    index = np.arange(count, dtype=np.float64)
    z = 1 - 2 * (index + .5) / count
    angle = index * np.pi * (3 - np.sqrt(5))
    normal = np.column_stack((np.sqrt(1 - z * z) * np.cos(angle), np.sqrt(1 - z * z) * np.sin(angle), z))
    points = np.concatenate([normal * 4, normal * 4 + [12, 0, 0]]).astype(np.float32)
    nodes = np.concatenate([np.arange(count) // 20, 4 + np.arange(count) // 20])
    repeat = np.arange(0, count * 2, 3)
    return {
        'x': np.zeros((8, 3), dtype=np.float32), 'side': np.repeat([0, 1], 4),
        'residue_key': np.asarray([f'A:{i}' for i in range(8)]),
        'point_pos': np.concatenate([points, points[repeat]]),
        'point_normal': np.concatenate([normal, normal, np.tile(normal, (2, 1))[repeat]]).astype(np.float32),
        'point_node_index': np.concatenate([nodes, nodes[repeat]]),
        'point_weight': np.concatenate([np.ones(160), np.zeros(len(repeat))]).astype(np.float32),
        'point_interface_weight': np.concatenate([np.zeros(160), np.full(len(repeat), .5)]).astype(np.float32),
    }


def test_duplicate_point_geometry_preserves_both_quadratures():
    graph = surface_fixture()
    for side in (0, 1):
        mapping = quadrature_mapping(graph, side)
        assert len(mapping['positions']) == 80
        original = mapping['quadrature_original_index'].numpy()
        assert np.array_equal(mapping['positions'][mapping['quadrature_inverse']].numpy(), graph['point_pos'][original])
        assert np.array_equal(mapping['quadrature_node_index'].numpy(), graph['point_node_index'][original])
        values = torch.ones((80, 2), requires_grad=True)
        for weight, original_weight in [('surface_weight', 'point_weight'), ('interface_weight', 'point_interface_weight')]:
            pooled, masses = pool_quadrature(values, mapping, 8, weight)
            expected = np.bincount(graph['point_node_index'][original], weights=graph[original_weight][original], minlength=8)
            np.testing.assert_allclose(masses.numpy(), expected)
            assert torch.all(pooled[masses > 0] == 1)
            assert torch.all(pooled[masses == 0] == 0)
            pooled.sum().backward(retain_graph=True)
        assert torch.isfinite(values.grad).all()
        assert values.grad.abs().sum() > 0


def test_inconsistent_duplicate_residue_ownership_is_rejected():
    graph = surface_fixture()
    graph['point_node_index'][-1] = 4
    with pytest.raises(ValueError, match='residue ownership'):
        quadrature_mapping(graph, 1)


