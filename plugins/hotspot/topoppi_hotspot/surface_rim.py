"""Native surface distance-to-rim descriptors for interface ranking.

This is a graph-edge approximation of a geodesic distance transform, not an
atomic solvent depth or energetic calculation. Distances follow retained native
surface faces; footprint pooling preserves parent vertex provenance.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import trimesh
from scipy import sparse
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree

from topoppi.atlas.footprints import source_atom_residue_labels
from topoppi.interface import chain_group, prepare_interface
from topoppi.mesh.provenance import connected_face_components

from .mesh_diffusion import native_footprint_pooling, pool_channels

CHANNELS = ('geodesic', 'euclidean', 'normalized_geodesic', 'normalized_euclidean',
            'geodesic_excess', 'deep_2A', 'deep_4A', 'deep_6A', 'no_boundary')


def rim_distance_fields(vertices, faces):
    """Multi-source shortest paths from physical boundary vertices.

    Each shared-edge face component uses its own rim and area. A shared vertex
    alone does not merge physically distinct patches. Unused vertices and
    closed components receive zeros; a separate channel marks closed-component
    vertices. Unique original edges avoid counting internal patch seams as rim.
    """
    vertices, faces = np.asarray(vertices, float), np.asarray(faces, int)
    values = np.zeros((len(vertices), len(CHANNELS)))
    if not len(faces):
        return values, {'n_boundary_vertices': 0, 'n_closed_components': 0}
    face_components = connected_face_components(trimesh.Trimesh(vertices=vertices, faces=faces, process=False))
    active = np.unique(faces)
    triangle = vertices[faces]
    face_area = .5 * np.linalg.norm(np.cross(triangle[:, 1] - triangle[:, 0],
                                            triangle[:, 2] - triangle[:, 0]), axis=1)
    closed, boundary = 0, set()
    vertex_component_count = np.zeros(len(vertices), int)
    for face_rows in face_components:
        current_faces = faces[face_rows]
        rows, inverse = np.unique(current_faces, return_inverse=True)
        vertex_component_count[rows] += 1
        local_faces = inverse.reshape(-1, 3)
        edges = np.sort(np.concatenate([local_faces[:, [0, 1]], local_faces[:, [1, 2]], local_faces[:, [2, 0]]]), axis=1)
        edges, counts = np.unique(edges, axis=0, return_counts=True)
        rim = np.unique(edges[counts == 1])
        boundary.update(map(int, rows[rim]))
        if not len(rim):
            values[rows, -1] = 1.
            closed += 1
            continue
        local_vertices = vertices[rows]
        lengths = np.linalg.norm(local_vertices[edges[:, 0]] - local_vertices[edges[:, 1]], axis=1)
        adjacency = sparse.coo_matrix((np.tile(lengths, 2),
                                       (np.r_[edges[:, 0], edges[:, 1]], np.r_[edges[:, 1], edges[:, 0]])),
                                      shape=(len(rows), len(rows))).tocsr()
        depth = dijkstra(adjacency, directed=False, indices=rim, min_only=True)
        euclidean = cKDTree(local_vertices[rim]).query(local_vertices)[0]
        area = face_area[face_rows].sum()
        scale = np.sqrt(area)
        values[rows] = np.column_stack([depth, euclidean, depth / scale, euclidean / scale,
                                       depth - euclidean, depth >= 2., depth >= 4., depth >= 6.,
                                       np.zeros(len(rows))])
    return values, {'n_boundary_vertices': len(boundary), 'n_closed_components': closed,
                    'n_face_components': len(face_components),
                    'n_shared_component_vertices': int((vertex_component_count > 1).sum()),
                    'n_interface_vertices': len(active), 'n_interface_faces': len(faces),
                    'interface_area_angstrom2': float(face_area.sum())}


def extract_surface_rim(structure_path, chains_a, chains_b, *, prepared=None):
    groups = sorted((chain_group(chains_a), chain_group(chains_b)))
    prepared = prepared or prepare_interface(structure_path, *groups)
    parts, diagnostics = [], []
    for side, group in enumerate(groups):
        partner = prepared.partner(group)
        mesh = partner.mesh
        source_faces = np.asarray(mesh.metadata['source_face_ids'], int)
        retained = {int(face) for patch in partner.patches for face in patch.metadata['source_face_ids']}
        faces = np.asarray(mesh.faces)[np.isin(source_faces, list(retained))]
        fields, report = rim_distance_fields(mesh.vertices, faces)
        labels = source_atom_residue_labels(partner.atoms)
        keys = list(dict.fromkeys(map(str, labels)))
        pooling = native_footprint_pooling(mesh, partner.patches, labels, keys)
        means, std, area = pool_channels(pooling, fields)
        extra = {'residue_key': keys, 'side': side}
        for index, channel in enumerate(CHANNELS):
            extra[f'surface_rim.mean.{channel}'] = means[:, index]
            extra[f'surface_rim.std.{channel}'] = std[:, index]
        extra['surface_rim.missing_footprint'] = (area == 0).astype(float)
        parts.append(pd.DataFrame(extra))
        diagnostics.append({'side': side, **report})
    return pd.concat(parts, ignore_index=True), {
        'protocol': 'native_surface_rim_distance_v2', 'labels_used': False,
        'operator': 'multi-source Dijkstra over unique retained parent-mesh edges, physical lengths',
        'rim': 'edges belonging to exactly one retained original face; no artificial cut/patch seam',
        'normalization': 'divide distances by square root of shared-edge face-component physical area; vertex-only contacts do not merge components',
        'thresholds_angstrom': [2., 4., 6.], 'pooling': 'native footprint area-weighted mean and SD',
        'interpretation': 'edge-graph geodesic approximation; not solvent accessibility or atomic depth',
        'partners': diagnostics,
    }
