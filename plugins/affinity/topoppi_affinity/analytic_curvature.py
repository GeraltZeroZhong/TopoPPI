"""Existing implicit-surface curvature formulas for TopoPPI Gaussian fields.

Gaussian derivatives are analytic; atom neighborhoods use the same axis-aligned
four-sigma truncation as SurfaceGenerator. Goldman (2005), equations 4.1/4.2,
provide Gaussian/mean curvature at the mesh sample positions. No curvature
clipping, triangle-area denominator, or change to the cached mesh is applied.
"""
from __future__ import annotations

import json

import numpy as np
from scipy.spatial import cKDTree

CURVATURE_PROTOCOL = "analytic_unit_peak_gaussian_box4_v1"
GRAPH_PROTOCOL = "native_residue_surface_graph_analytic_gaussian_v1"


def gaussian_derivatives(points, atom_coords, *, sigma=1., truncate_sigma=4., chunk_size=2048):
    """Density, gradient, Hessian from sparse nearby atom-point pairs.

The native truncation is a box, not a Euclidean ball. Derivatives describe the
smooth Gaussian terms active locally; the numerical hard cutoff itself has
discontinuous boundaries. This distinction and mesh interpolation error are
recorded in the protocol rather than hidden by a fitted curvature transform.
"""
    points, atom_coords = np.asarray(points, dtype=np.float64), np.asarray(atom_coords, dtype=np.float64)
    tree = cKDTree(atom_coords)
    density = np.empty(len(points), dtype=np.float64)
    gradient = np.empty((len(points), 3), dtype=np.float64)
    hessian = np.empty((len(points), 3, 3), dtype=np.float64)
    for start in range(0, len(points), chunk_size):
        query = points[start:start + chunk_size]
        neighbors = tree.query_ball_point(query, truncate_sigma * sigma, p=np.inf)
        counts = np.fromiter((len(items) for items in neighbors), dtype=np.int64)
        rows = np.repeat(np.arange(len(query)), counts)
        columns = np.concatenate(neighbors).astype(np.int64)
        delta = query[rows] - atom_coords[columns]
        phi = np.exp(-np.einsum("ij,ij->i", delta, delta) / (2. * sigma**2))
        total = np.bincount(rows, weights=phi, minlength=len(query))
        density[start:start + len(query)] = total
        for axis in range(3):
            gradient[start:start + len(query), axis] = np.bincount(
                rows, weights=-phi * delta[:, axis] / sigma**2, minlength=len(query))
            for other in range(axis, 3):
                values = np.bincount(rows, weights=phi * delta[:, axis] * delta[:, other] / sigma**4,
                                     minlength=len(query))
                if axis == other:
                    values -= total / sigma**2
                hessian[start:start + len(query), axis, other] = values
                hessian[start:start + len(query), other, axis] = values
    return density, gradient, hessian


def implicit_curvatures(gradient, hessian):
    """Goldman Gaussian and signed mean curvature, without Hessian inversion."""
    cofactors = np.stack([np.cross(hessian[:, 1], hessian[:, 2]),
                          np.cross(hessian[:, 2], hessian[:, 0]),
                          np.cross(hessian[:, 0], hessian[:, 1])], axis=1)
    norm2 = np.einsum("ij,ij->i", gradient, gradient)
    if np.any(norm2 == 0):
        raise ValueError("Implicit curvature is undefined at zero-gradient points")
    gaussian = np.einsum("ni,nij,nj->n", gradient, cofactors, gradient) / norm2**2
    mean = (np.einsum("ni,nij,nj->n", gradient, hessian, gradient)
            - norm2 * np.trace(hessian, axis1=1, axis2=2)) / (2. * norm2**1.5)
    return gaussian, mean


def gaussian_curvature(points, atom_coords, *, sigma=1., truncate_sigma=4.):
    density, gradient, hessian = gaussian_derivatives(points, atom_coords, sigma=sigma, truncate_sigma=truncate_sigma)
    gaussian, _ = implicit_curvatures(gradient, hessian)
    return gaussian, {"protocol": CURVATURE_PROTOCOL, "sigma_angstrom": sigma,
                      "truncate_sigma_box": truncate_sigma, "gradient_norm_min": float(np.linalg.norm(gradient, axis=1).min()),
                      "density_min": float(density.min()), "density_max": float(density.max()),
                      "gaussian_min": float(gaussian.min()), "gaussian_max": float(gaussian.max()),
                      "sampling": "native mesh vertices; no level-set projection or clipping"}


def replace_graph_curvature(graph, structure_path):
    """Replace only the two Gaussian-curvature columns on a supplied WT/MT graph."""
    from topoppi.atlas.footprints import mesh_vertex_residue_labels, source_atom_residue_labels
    from topoppi.config import SurfaceConfig
    from topoppi.io.io_loader import PDBLoader
    from topoppi.mesh.surface import SurfaceGenerator

    from topoppi_affinity.features import _group_data
    from topoppi_affinity.graph_features import FEATURE_NAMES, _vertex_mass

    result = {key: value.copy() for key, value in graph.items()}
    mean_column = FEATURE_NAMES.index("surface.gaussian_curvature_mean")
    std_column = FEATURE_NAMES.index("surface.gaussian_curvature_std")
    result["x"][:, [mean_column, std_column]] = 0.
    lookup = {str(key): index for index, key in enumerate(graph["residue_key"])}
    loader = PDBLoader(str(structure_path))
    reports = []
    for partner in (0, 1):
        chains = tuple(np.unique(graph["chain"][graph["side"] == partner]).tolist())
        coordinates, atoms = _group_data(loader, chains)
        config = SurfaceConfig()
        mesh = SurfaceGenerator(coordinates, config=config).generate_mesh()
        if mesh is None:
            raise ValueError(f"Cannot reconstruct native side {partner} mesh for analytic curvature")
        mass = _vertex_mass(mesh)
        labels = mesh_vertex_residue_labels(mesh, source_atom_residue_labels(atoms))
        nodes = np.asarray([lookup[str(label)] for label in labels], dtype=np.int64)
        curvature, report = gaussian_curvature(mesh.vertices, coordinates, sigma=config.sigma)
        for node in np.unique(nodes):
            selected = nodes == node
            weights, values = mass[selected], curvature[selected]
            if not weights.sum():
                continue
            mean = np.average(values, weights=weights)
            result["x"][node, mean_column] = mean
            result["x"][node, std_column] = np.sqrt(np.average((values - mean)**2, weights=weights))
        reports.append(report)
    metadata = json.loads(str(graph["metadata_json"]))
    metadata.update(protocol=GRAPH_PROTOCOL, curvature_method="analytic_gaussian", analytic_curvature=reports,
                    replaced_columns=[FEATURE_NAMES[mean_column], FEATURE_NAMES[std_column]])
    result["metadata_json"] = np.asarray(json.dumps(metadata))
    return result


def update_pooled_curvature(vector, names, graph, *, target, mutations=""):
    """Change exactly the pooled slots affected by the two replaced graph columns."""
    from topoppi_affinity.graph_features import mutant_node_features

    result = np.asarray(vector).copy()
    lookup = {str(name): index for index, name in enumerate(names)}
    curvature = graph["x"][:, 34:36].astype(np.float64)
    if target == "pkd":
        for statistic in ("mean", "sum", "std"):
            left, right = [getattr(curvature[graph["side"] == partner], statistic)(axis=0) for partner in (0, 1)]
            values = np.concatenate([left + right, np.abs(left - right)])
            for index, value in zip((4, 5, 19, 20), values, strict=True):
                result[lookup[f"surface.{statistic}.sum_absdiff.{index}"]] = value
    else:
        mutant, _, mask = mutant_node_features(graph, mutations)
        delta = mutant[mask, :26].astype(np.float64) - graph["x"][mask, :26].astype(np.float64)
        products = np.einsum("ni,nj->ij", delta, curvature[mask])
        for chemical in range(26):
            for column in range(2):
                result[lookup[f"surface.chemical_delta_context.{chemical * 15 + 4 + column}"]] = products[chemical, column]
    return result
