"""Existing bilinear pooling applied to interacting TopoPPI residue pairs.

Lin et al., ICCV 2015: outer-product aggregation, signed square root and L2
normalization. The two streams here are residues on opposite binding partners.
This input adaptation retains the original frozen ESMC pooled control.
"""
from __future__ import annotations

import numpy as np

from topoppi_affinity.embedding_baselines import pooled_features


def symmetric_bilinear(nodes, edges, weights):
    """Normalized symmetric cross-stream outer product, with no ordered partner."""
    left, right = edges
    matrix = nodes[left].T @ (nodes[right] * weights[:, None])
    matrix = (matrix + matrix.T) * .5
    # sqrt(2) accounts for the two identical off-diagonal matrix entries.
    i, j = np.triu_indices(matrix.shape[0])
    vector = np.sign(matrix[i, j]) * np.sqrt(np.abs(matrix[i, j]))
    vector[i != j] *= np.sqrt(2.)
    return vector / max(np.linalg.norm(vector), 1e-12)


def pair_state(data, projection, *, mutant=False):
    x = (data.x_mut if mutant else data.x).numpy().astype(float)
    embedding = (data.esm_mut if mutant else data.esm).numpy().astype(float) @ projection
    embedding /= np.maximum(np.linalg.norm(embedding, axis=1, keepdims=True), 1e-12)
    chemistry = np.sign(x[:, 21:26]) * np.log1p(np.abs(x[:, 21:26]))
    chemistry /= np.maximum(np.linalg.norm(chemistry, axis=1, keepdims=True), 1.)
    surface = np.sign(x[:, 30:]) * np.log1p(np.abs(x[:, 30:]))
    surface /= np.maximum(np.linalg.norm(surface, axis=1, keepdims=True), 1.)
    common = np.concatenate([embedding, x[:, :21], chemistry], axis=1)
    contact_nodes = np.concatenate([common, np.zeros_like(surface)], axis=1)
    surface_nodes = np.concatenate([common, surface], axis=1)
    edges, attr = data.edge_index.numpy(), data.edge_attr.numpy()
    keep = (attr[:, 3] > 0) & (edges[0] < edges[1])
    edges, attr = edges[:, keep], attr[keep]
    weights = np.log1p(attr[:, 5]) * np.exp(-attr[:, 0] ** 2 / 128.)
    area = np.sqrt(np.maximum(x[edges[0], 31] * x[edges[1], 31], 0.))
    return np.concatenate([symmetric_bilinear(contact_nodes, edges, weights),
                           symmetric_bilinear(surface_nodes, edges, weights * area)])


def bilinear_features(data, projection, *, target):
    base, names = pooled_features(data, target=target, use_surface=True)
    state = pair_state(data, projection)
    if target == 'ddg':
        state = pair_state(data, projection, mutant=True) - state
    half = len(state) // 2
    names = np.concatenate([names, np.asarray([f'pair.contact.{i}' for i in range(half)] +
                                               [f'surface.pair.{i}' for i in range(half)])])
    return np.concatenate([base, state]).astype(np.float32), names, len(base)
