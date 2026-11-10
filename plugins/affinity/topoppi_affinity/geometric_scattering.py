"""Geometric graph scattering adapted to TopoPPI residue/footprint graphs.

Gao, Wolf & Hirn, ICML 2019, https://proceedings.mlr.press/v97/gao19e.html.
The author's column-normalized lazy walk and diffusion wavelets are evaluated
with sparse matrix-vector products. Four scales and two orders are fixed before
evaluation. No scattering weights or PLM projection are fitted to labels.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

SCALES = (1, 2, 4, 8)


def lazy_walk(edge_index, weights, n_nodes):
    adjacency = sp.csr_matrix((weights, edge_index), shape=(n_nodes, n_nodes), dtype=np.float64)
    adjacency = adjacency.maximum(adjacency.T)
    degree = np.asarray(adjacency.sum(axis=0)).reshape(-1)
    inverse = np.zeros_like(degree)
    np.divide(1., degree, out=inverse, where=degree > 0)
    return .5 * (sp.eye(n_nodes, format='csr') + adjacency @ sp.diags(inverse))


def wavelet_signals(walk, signals, scales=SCALES):
    states = {}
    wanted = set(scales) | {2 * t for t in scales}
    current = np.asarray(signals, dtype=np.float64)
    for step in range(1, max(wanted) + 1):
        current = walk @ current
        if step in wanted:
            states[step] = current
    return [states[t] - states[2 * t] for t in scales]


def scattering_signals(walk, signals):
    first = [np.abs(value) for value in wavelet_signals(walk, signals)]
    second = []
    for i, value in enumerate(first[:-1]):
        second.extend(np.abs(w) for w in wavelet_signals(walk, value, SCALES[i+1:]))
    return np.concatenate([signals, *first, *second], axis=1)


def _physical(x):
    return np.concatenate([x[:, :21], np.sign(x[:, 21:]) * np.log1p(np.abs(x[:, 21:]))], axis=1)


def representation(data, *, use_surface, target, projection):
    x, xm = data.x.numpy().copy(), data.x_mut.numpy().copy()
    edges, attributes = data.edge_index.numpy(), data.edge_attr.numpy()
    if not use_surface:
        x[:, 30:], xm[:, 30:] = 0., 0.
        keep = ((attributes[:, 2] > 0) & (attributes[:, 0] <= 10.)) | (attributes[:, 3] > 0) | (attributes[:, 4] > 0)
        edges, attributes = edges[:, keep], attributes[keep]
    weights = np.exp(-attributes[:, 0] / 10.)
    if use_surface:
        weights += np.log1p(attributes[:, 6])
    walk = lazy_walk(edges, weights, len(x))
    wt = np.concatenate([_physical(x), data.esm.numpy() @ projection], axis=1)
    wt = scattering_signals(walk, wt)
    contact, side, mutation = data.contact_mask.numpy(), data.side.numpy(), data.mutation_mask.numpy()
    parts = []

    def mean(values, mask):
        return values[mask].mean(0) if mask.any() else np.zeros(values.shape[1])

    if target == 'pkd':
        for mask in [np.ones(len(wt), dtype=bool), contact]:
            for moment in [1, 2]:
                values = wt if moment == 1 else wt ** 2
                a, b = mean(values, mask & (side == 0)), mean(values, mask & (side == 1))
                parts.extend([a + b, np.abs(a - b)])
    else:
        mt = np.concatenate([_physical(xm), data.esm_mut.numpy() @ projection], axis=1)
        mt = scattering_signals(walk, mt)
        delta, context = mt - wt, (mt + wt) / 2.
        parts.extend([delta.mean(0), mean(delta, contact), delta[mutation].sum(0),
                      (delta[mutation] * context[mutation]).sum(0)])
    return np.concatenate(parts).astype(np.float32)
