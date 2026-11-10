"""Fixed cross-partner context and paired native-surface descriptors.

The normal/distance field is an SC-inspired descriptor on TopoPPI Gaussian
meshes, not the CCP4 Sc statistic. Contact products are descriptive channels,
not interaction energies, mutation predictions, or physical decomposition.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from topoppi.atlas.footprints import source_atom_residue_labels
from topoppi.interface import chain_group, prepare_interface

from .features import AA, CHEMICAL
from .mesh_diffusion import native_footprint_pooling, pool_channels


def cross_partner_context(features, graphs):
    """Atom-pair-count weighted opposite-partner residue context, no labels."""
    frame = features.copy().reset_index(drop=True)
    source_columns = [f'aa.{a}' for a in AA + 'X'] + list(CHEMICAL)
    extras = {f'contact_context.mean.{c}': np.zeros(len(frame)) for c in source_columns}
    for column in CHEMICAL:
        extras[f'contact_context.product.{column}'] = np.zeros(len(frame))
    for complex_id, part in frame.groupby('complex_id', sort=True):
        graph = graphs[str(complex_id)]
        keys = np.asarray(graph['residue_keys']).astype(str)
        ordered = part.set_index('residue_key').loc[keys]
        source = ordered[source_columns].to_numpy(float)
        side = ordered.side.to_numpy()
        edges = np.asarray(graph['cross_edges'], int).reshape(-1, 2)
        weights = np.asarray(graph['cross_weights'], float)
        if len(edges) and np.any(side[edges[:, 0]] == side[edges[:, 1]]):
            raise ValueError('Cross-partner context received within-partner contacts')
        totals, count = np.zeros_like(source), np.zeros(len(source))
        if len(edges):
            a, b = edges.T
            np.add.at(totals, a, source[b] * weights[:, None])
            np.add.at(totals, b, source[a] * weights[:, None])
            np.add.at(count, a, weights)
            np.add.at(count, b, weights)
        means = totals / np.maximum(count[:, None], 1.)
        positions = pd.Series(part.index.to_numpy(), index=part.residue_key).loc[keys].to_numpy()
        for i, column in enumerate(source_columns):
            extras[f'contact_context.mean.{column}'][positions] = means[:, i]
        for column in CHEMICAL:
            i = source_columns.index(column)
            extras[f'contact_context.product.{column}'][positions] = source[:, i] * means[:, i]
    return pd.concat([frame, pd.DataFrame(extras)], axis=1), list(extras)


def paired_vertex_fields(vertices, normals, other_vertices, other_normals, *, width=1.):
    """Return invariant nearest-opposite-vertex geometric fields.

    The stored mesh normals are used consistently. Current native Gaussian
    meshes have negative signed volume and inward winding. Reversing both
    normal fields leaves their dot product unchanged, while oriented_gap
    reverses sign. It must be interpreted using the stored winding, not as an
    outward gap or physical solvent depth.
    This uses nearest vertices and an area mean, unlike CCP4's full surface
    sampling, rim exclusion and two-sided median statistic.
    """
    vertices, normals = np.asarray(vertices, float), np.asarray(normals, float)
    other_vertices, other_normals = np.asarray(other_vertices, float), np.asarray(other_normals, float)
    distance, nearest = cKDTree(other_vertices).query(vertices)
    dot = np.clip(np.einsum('ij,ij->i', normals, other_normals[nearest]), -1., 1.)
    separation = other_vertices[nearest] - vertices
    oriented_gap = np.einsum('ij,ij->i', separation, normals)
    attenuation = np.exp(-.5 * (distance / width) ** 2)
    return np.column_stack([distance, -dot, -dot * attenuation, oriented_gap,
                            (dot < 0).astype(float), attenuation])


def extract_paired_surface(structure_path, chains_a, chains_b, *, prepared=None):
    groups = sorted((chain_group(chains_a), chain_group(chains_b)))
    prepared = prepared or prepare_interface(structure_path, *groups)
    partners = [prepared.partner(group) for group in groups]
    channel_names = ['gap', 'normal_opposition', 'weighted_normal_opposition',
                     'oriented_gap', 'opposing_normal_fraction', 'proximity']
    results, diagnostics = [], []
    for side, partner in enumerate(partners):
        other = partners[1 - side]
        labels = source_atom_residue_labels(partner.atoms)
        keys = list(dict.fromkeys(map(str, labels)))
        values = paired_vertex_fields(partner.mesh.vertices, partner.mesh.vertex_normals,
                                      other.mesh.vertices, other.mesh.vertex_normals)
        pooling = native_footprint_pooling(partner.mesh, partner.patches, labels, keys)
        means, std, area = pool_channels(pooling, values)
        extra = {'residue_key': keys, 'side': side}
        for i, channel in enumerate(channel_names):
            extra[f'surface_pair.mean.{channel}'] = means[:, i]
            extra[f'surface_pair.std.{channel}'] = std[:, i]
        extra['surface_pair.missing_footprint'] = (area == 0).astype(float)
        results.append(pd.DataFrame(extra))
        signed_volume = float(partner.mesh.volume)
        diagnostics.append({'side': side, 'n_vertices': len(partner.mesh.vertices),
                            'mesh_signed_volume': signed_volume,
                            'stored_normal_orientation': ('inward_negative_signed_volume' if signed_volume < 0
                                                          else 'outward_positive_signed_volume' if signed_volume > 0
                                                          else 'undetermined_zero_signed_volume'),
                            'n_missing_footprints': int((area == 0).sum())})
    return pd.concat(results, ignore_index=True), {
        'protocol': 'paired_native_gaussian_surface_context_v1', 'labels_used': False,
        'normal_field': '-dot(stored_normal, nearest_opposite_stored_normal) * exp(-distance^2/2)',
        'normal_convention': 'stored parent-mesh winding; current native meshes have inward normals (negative signed volume); no orientation correction applied',
        'oriented_gap_convention': '(nearest_opposite_vertex - vertex) dot stored_normal; sign reverses if this normal is reversed; not an outward gap or solvent depth',
        'normal_opposition_invariance': 'dot-product channels are unchanged when both partner normal fields are reversed',
        'width_angstrom': 1., 'pooling': 'native footprint face-corner physical area mean and SD',
        'interpretation': 'SC-inspired local descriptor, not CCP4 Sc or an energy',
        'reference': 'https://www.ccp4.ac.uk/html/sc.html', 'partners': diagnostics,
    }
