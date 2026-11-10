"""Author DiffusionNet operators and layers on TopoPPI surface quadrature points.

Coordinates retain Angstrom units. Operators use the original point-cloud
Laplacian, mass matrix and tangent-gradient implementation. Two quadratures
(whole surface and interface) share geometry but retain independent weights.
No affinity labels, sequence embeddings or split assignments enter this cache.
"""
from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import numpy as np
import torch
from torch import nn

DEFAULT_AUTHOR_REPO = Path('external/repos/diffusion-net')
OPERATOR_PROTOCOL = 'author_diffusionnet_pointcloud_angstrom_v1'
AUTHOR_REVISION = 'b1019b049597711d10259ce28250f7dfc4335a2b'


def author_module(repo=DEFAULT_AUTHOR_REPO):
    """Import the separately downloaded MIT-licensed author implementation."""
    source = Path(repo) / 'src'
    if not (source / 'diffusion_net' / 'geometry.py').is_file():
        raise FileNotFoundError(f'DiffusionNet author source missing: {source}')
    if str(source) not in sys.path:
        sys.path.insert(0, str(source))
    return importlib.import_module('diffusion_net')


def quadrature_mapping(graph, side):
    """Deduplicate geometry without discarding either area measure.

    Whole-surface and interface samples may select the same mesh vertex. The
    unique geometry is used once in the eigensolve; inverse indices restore
    every original quadrature entry and its original residue ownership.
    """
    point_nodes = np.asarray(graph['point_node_index'], dtype=np.int64)
    original_indices = np.flatnonzero(np.asarray(graph['side'])[point_nodes] == side)
    coords = np.asarray(graph['point_pos'])[original_indices]
    unique, first, inverse = np.unique(coords, axis=0, return_index=True, return_inverse=True)
    nodes = point_nodes[original_indices]
    if not np.array_equal(nodes, nodes[first][inverse]):
        raise ValueError('Identical surface coordinates have inconsistent residue ownership')
    normals = np.asarray(graph['point_normal'])[original_indices[first]].copy()
    lengths = np.linalg.norm(normals, axis=1, keepdims=True)
    if np.any(lengths == 0):
        raise ValueError('Surface point has a zero normal')
    normals /= lengths
    return {
        'positions': torch.as_tensor(unique, dtype=torch.float32),
        'normals': torch.as_tensor(normals, dtype=torch.float32),
        'node_index': torch.as_tensor(nodes[first], dtype=torch.long),
        'quadrature_node_index': torch.as_tensor(nodes, dtype=torch.long),
        'quadrature_inverse': torch.as_tensor(inverse, dtype=torch.long),
        'quadrature_original_index': torch.as_tensor(original_indices, dtype=torch.long),
        'surface_weight': torch.as_tensor(np.asarray(graph['point_weight'])[original_indices], dtype=torch.float32),
        'interface_weight': torch.as_tensor(np.asarray(graph['point_interface_weight'])[original_indices], dtype=torch.float32),
        'side': int(side),
    }


def compute_surface_operators(graph, *, source, complex_id, k_eig=64, author_repo=DEFAULT_AUTHOR_REPO):
    """Build both partners independently with the author's exact point-cloud path."""
    author = author_module(author_repo)
    sides = []
    for side in np.unique(graph['side']):
        values = quadrature_mapping(graph, side)
        count = len(values['positions'])
        if count <= max(30, k_eig):
            raise ValueError(f'{source}/{complex_id} side {side}: {count} unique points insufficient for {k_eig} modes')
        tensors = author.geometry.compute_operators(
            values['positions'], torch.empty((0, 3), dtype=torch.long),
            k_eig=k_eig, normals=values['normals'],
        )
        for key, tensor in zip(('frames', 'mass', 'L', 'evals', 'evecs', 'gradX', 'gradY'), tensors, strict=True):
            values[key] = tensor.coalesce() if tensor.is_sparse else tensor
        sides.append(values)
    return {
        'protocol': OPERATOR_PROTOCOL, 'source': str(source), 'complex_id': str(complex_id),
        'author_revision': AUTHOR_REVISION, 'k_eig': int(k_eig),
        'n_nodes': len(graph['x']), 'coordinate_units': 'angstrom',
        'residue_key': [str(value) for value in graph['residue_key']],
        'sides': sides,
    }


def load_surface_operators(path, device='cpu'):
    """Load tensor-only operators; callers may lazily move them to their device."""
    return torch.load(path, map_location=device, weights_only=True)


def operators_to(bundle, device):
    """Move one cached complex to a device, preserving ordinary identity metadata."""
    output = {key: value for key, value in bundle.items() if key != 'sides'}
    output['sides'] = [
        {key: value.to(device) if isinstance(value, torch.Tensor) else value for key, value in side.items()}
        for side in bundle['sides']
    ]
    return output


def pool_quadrature(point_embeddings, side, n_nodes, weight_key):
    """Area-weighted residue mean with exactly zero output for uncovered nodes."""
    weights = side[weight_key].to(point_embeddings.dtype)
    values = point_embeddings[side['quadrature_inverse']]
    nodes = side['quadrature_node_index']
    pooled = point_embeddings.new_zeros((n_nodes, point_embeddings.shape[-1]))
    mass = point_embeddings.new_zeros(n_nodes)
    pooled.index_add_(0, nodes, values * weights[:, None])
    mass.index_add_(0, nodes, weights)
    pooled = pooled / mass.clamp_min(torch.finfo(pooled.dtype).tiny)[:, None]
    return pooled, mass


class SurfaceEncoder(nn.Module):
    """Shared original DiffusionNet on each partner, followed by footprint pooling.

    Input: full residue tensor [N,C_in], matching bundle['residue_key'] order.
    Output: dict surface/interface [N,C_out] and corresponding area masses [N].
    Geometry is fixed for WT/MT; mutation chemistry may change the input tensor.
    Use the complete graph before any residue cropping. No coordinates enter
    the learned input, avoiding an artificial global rotation dependence.
    """

    def __init__(self, c_in, c_out, width=32, n_blocks=2, dropout=False, author_repo=DEFAULT_AUTHOR_REPO):
        super().__init__()
        author = author_module(author_repo)
        self.network = author.layers.DiffusionNet(
            C_in=c_in, C_out=c_out, C_width=width, N_block=n_blocks,
            outputs_at='vertices', dropout=dropout, diffusion_method='spectral',
            with_gradient_features=True, with_gradient_rotations=True,
        )
        self.c_out = c_out

    def forward(self, node_features_full, operator_bundle):
        if node_features_full.ndim != 2 or len(node_features_full) != operator_bundle['n_nodes']:
            raise ValueError('SurfaceEncoder requires all residues in the cached graph order')
        n_nodes = len(node_features_full)
        result = {
            'surface': node_features_full.new_zeros((n_nodes, self.c_out)),
            'interface': node_features_full.new_zeros((n_nodes, self.c_out)),
            'surface_mass': node_features_full.new_zeros(n_nodes),
            'interface_mass': node_features_full.new_zeros(n_nodes),
        }
        for side in operator_bundle['sides']:
            point_features = node_features_full[side['node_index']]
            embedding = self.network(
                point_features, side['mass'], L=side['L'], evals=side['evals'],
                evecs=side['evecs'], gradX=side['gradX'], gradY=side['gradY'],
            )
            for name, weight in (('surface', 'surface_weight'), ('interface', 'interface_weight')):
                pooled, mass = pool_quadrature(embedding, side, n_nodes, weight)
                result[name] = result[name] + pooled
                result[f'{name}_mass'] = result[f'{name}_mass'] + mass
        return result


def write_provenance(path, author_repo=DEFAULT_AUTHOR_REPO):
    """Document the transplanted algorithm separately from a paper reproduction."""
    import importlib.metadata

    payload = {
        'protocol': OPERATOR_PROTOCOL,
        'source_url': 'https://github.com/nmwsharp/diffusion-net',
        'author_repo': str(author_repo), 'author_revision': AUTHOR_REVISION, 'license': 'MIT',
        'implementation': 'Unmodified author compute_operators and DiffusionNet; TopoPPI quadrature mapping and residue pooling',
        'coordinate_units': 'angstrom', 'faces': 'empty; original point-cloud mode',
        'geometry': 'WT native TopoPPI surface; same fixed operators for WT and MT',
        'weights': 'Author robust-Laplacian mass for spectral diffusion; original TopoPPI whole/interface quadrature weights for residue pooling',
        'pretraining': 'None; these operators and author modules contain no pretrained model weights',
        'training_labels_in_operator_cache': False,
        'robust_laplacian_version': importlib.metadata.version('robust-laplacian'),
        'potpourri3d_version': importlib.metadata.version('potpourri3d'),
        'notes': 'Algorithm transfer, not a reproduction of DiffusionNet or AtomSurf protein benchmark results. No coordinate rescaling.',
    }
    Path(path).write_text(json.dumps(payload, indent=2) + '\n')
