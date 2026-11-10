"""Physical area and native boundary integrals shared by geometry consumers."""

from __future__ import annotations

from collections import defaultdict

import numpy as np


def vertex_area_masses(mesh):
    """Assign each triangle's area equally to its three incident vertices."""
    return np.bincount(np.asarray(mesh.faces).reshape(-1),
                       weights=np.repeat(np.asarray(mesh.area_faces) / 3., 3),
                       minlength=len(mesh.vertices))


def residue_dual_boundaries(mesh, node_indices):
    """Integrate native barycentric-dual boundaries between residue indices."""
    faces = np.asarray(mesh.faces)
    centers = np.asarray(mesh.triangles_center)
    result = defaultdict(float)
    for corner in range(3):
        left, right = faces[:, corner], faces[:, (corner + 1) % 3]
        changed = node_indices[left] != node_indices[right]
        left, right = left[changed], right[changed]
        if not len(left):
            continue
        lengths = np.linalg.norm((mesh.vertices[left] + mesh.vertices[right]) / 2. - centers[changed], axis=1)
        pairs = np.sort(np.stack([node_indices[left], node_indices[right]], axis=1), axis=1)
        unique, inverse = np.unique(pairs, return_inverse=True, axis=0)
        for pair, length in zip(unique, np.bincount(inverse, weights=lengths), strict=True):
            result[tuple(map(int, pair))] += float(length)
    return result
