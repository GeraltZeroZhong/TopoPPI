"""Label-free residue graphs and native surface point clouds for affinity models.

All nodes use the supplied WT structure. Mutant inputs replace amino-acid
identities; coordinates and Gaussian surfaces remain the WT geometry. No
experimental affinity, annotation weight, split, or mutation target enters
this module.
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from topoppi import prepare_interface, residue_dual_boundaries, vertex_area_masses
from topoppi.atlas.footprints import mesh_vertex_residue_labels, source_atom_residue_labels
from topoppi.config import SurfaceConfig

from topoppi_affinity.features import (
    AA,
    PROPERTIES,
    _chains,
    extract_complex_features,
    parse_mutations,
    residue_properties,
)

# Previous private helper names retain their numerical contract.
_vertex_mass = vertex_area_masses
_dual_boundaries = residue_dual_boundaries

GRAPH_PROTOCOL = "native_residue_surface_graph_v1"
CHEMICAL_FEATURES = [f"aa.{aa}" for aa in AA + "X"] + [f"chemical.{p}" for p in PROPERTIES]
CONTACT_FEATURES = ["contact.partner_distance", "contact.degree", "contact.atom_pairs", "contact.neighbor_count"]
SURFACE_FEATURES = [
    "surface.area", "surface.footprint_area", "surface.footprint_fraction", "surface.patch_count",
    "surface.gaussian_curvature_mean", "surface.gaussian_curvature_std", "surface.bending_density",
    "surface.normal_dispersion", "surface.extent_eig0", "surface.extent_eig1", "surface.extent_eig2",
    "surface.footprint_normal_dispersion", "surface.footprint_extent_eig0",
    "surface.footprint_extent_eig1", "surface.footprint_extent_eig2",
]
FEATURE_NAMES = CHEMICAL_FEATURES + CONTACT_FEATURES + SURFACE_FEATURES
EDGE_FEATURE_NAMES = [
    "distance.ca", "distance.centroid", "edge.same_side", "edge.cross_contact", "edge.peptide",
    "contact.atom_pairs", "surface.boundary_length", "surface.footprint_boundary_length",
]
BACKBONE_ATOMS = ("N", "CA", "C", "O")




def _vertex_bending(mesh):
    """Unsigned integrated dihedral bending assigned equally to edge endpoints."""
    values = np.zeros(len(mesh.vertices), dtype=np.float64)
    adjacency = np.asarray(mesh.face_adjacency)
    if len(adjacency):
        edges = np.asarray(mesh.face_adjacency_edges)
        dot = (mesh.face_normals[adjacency[:, 0]] * mesh.face_normals[adjacency[:, 1]]).sum(axis=1)
        integral = np.arccos(np.clip(dot, -1., 1.)) * np.linalg.norm(
            mesh.vertices[edges[:, 0]] - mesh.vertices[edges[:, 1]], axis=1)
        np.add.at(values, edges.reshape(-1), np.repeat(integral / 2., 2))
    return values




def _weighted_shape(points, normals, weights):
    total = float(weights.sum())
    if total <= 0:
        return np.zeros(3), np.zeros(3), 0., np.zeros(3)
    weights = weights / total
    center = (points * weights[:, None]).sum(axis=0)
    normal = (normals * weights[:, None]).sum(axis=0)
    centered = points - center
    eigenvalues = np.linalg.eigvalsh((centered * weights[:, None]).T @ centered)
    return center, normal, float(1. - np.linalg.norm(normal)), np.maximum(eigenvalues, 0.)


def _mass_sample(weights, maximum):
    """Deterministic equal-mass quadrature; repeated vertices retain their mass."""
    if maximum <= 0:
        return np.empty(0, dtype=np.int64), np.empty(0)
    active = np.flatnonzero(weights > 0)
    if len(active) <= maximum:
        return active, weights[active]
    cumulative = np.cumsum(weights)
    mass_per_point = cumulative[-1] / maximum
    selected = np.searchsorted(cumulative, (np.arange(maximum) + .5) * mass_per_point)
    selected, counts = np.unique(selected, return_counts=True)
    return selected, counts * mass_per_point


def _surface_data(prepared, lookup, node_count, maximum, curvature_method="angle_defect"):
    mesh, coords, atoms = prepared.mesh, prepared.coords, prepared.atoms
    source_labels = source_atom_residue_labels(atoms)
    labels = mesh_vertex_residue_labels(mesh, source_labels)
    nodes = np.asarray([lookup[str(label)] for label in labels], dtype=np.int64)
    mass = _vertex_mass(mesh)
    if curvature_method == "analytic_gaussian":
        from topoppi_affinity.analytic_curvature import gaussian_curvature

        curvature, _ = gaussian_curvature(mesh.vertices, coords, sigma=SurfaceConfig().sigma)
    elif curvature_method == "angle_defect":
        curvature = np.asarray(mesh.vertex_defects) / np.maximum(mass, 1e-12)
    else:
        raise ValueError(f"Unknown curvature method: {curvature_method}")
    bending = _vertex_bending(mesh)
    patches = prepared.patches
    interface_mass = np.zeros(len(mesh.vertices), dtype=np.float64)
    footprint_boundaries = defaultdict(float)
    for patch in patches:
        source_vertices = np.asarray(patch.metadata["source_vertex_ids"], dtype=np.int64)
        np.add.at(interface_mass, source_vertices, _vertex_mass(patch))
        patch_nodes = np.asarray([lookup[str(label)] for label in mesh_vertex_residue_labels(patch, source_labels)])
        for pair, length in _dual_boundaries(patch, patch_nodes).items():
            footprint_boundaries[pair] += length
    geometry = np.zeros((node_count, 11), dtype=np.float64)
    normals = np.zeros((node_count, 3), dtype=np.float64)
    footprint_centers = np.zeros((node_count, 3), dtype=np.float64)
    footprint_normals = np.zeros((node_count, 3), dtype=np.float64)
    for node in np.unique(nodes):
        selected = nodes == node
        weights = mass[selected]
        _, normal, dispersion, eigenvalues = _weighted_shape(mesh.vertices[selected], mesh.vertex_normals[selected], weights)
        curv = curvature[selected]
        curv_mean = np.average(curv, weights=weights)
        curv_std = np.sqrt(np.average((curv - curv_mean) ** 2, weights=weights))
        foot_center, foot_normal, foot_dispersion, foot_eigenvalues = _weighted_shape(
            mesh.vertices[selected], mesh.vertex_normals[selected], interface_mass[selected])
        geometry[node] = [curv_mean, curv_std, bending[selected].sum() / weights.sum(),
                          dispersion, *eigenvalues, foot_dispersion, *foot_eigenvalues]
        normals[node] = normal
        footprint_centers[node] = foot_center
        footprint_normals[node] = foot_normal
    # Half the quota samples complete surface, half samples the accepted interface.
    # These are independent quadratures; each has its own physical area weight.
    whole_indices, whole_weights = _mass_sample(mass, maximum // 2)
    foot_indices, foot_weights = _mass_sample(interface_mass, maximum - maximum // 2)
    indices = np.concatenate([whole_indices, foot_indices])
    points = {
        "point_pos": mesh.vertices[indices], "point_normal": mesh.vertex_normals[indices],
        "point_node_index": nodes[indices],
        "point_weight": np.concatenate([whole_weights, np.zeros(len(foot_indices))]),
        "point_interface_weight": np.concatenate([np.zeros(len(whole_indices)), foot_weights]),
    }
    report = {"generation": prepared.surface_report, "topology_status": prepared.topology_report["status"],
              "n_vertices": len(mesh.vertices), "n_faces": len(mesh.faces), "n_points": len(indices),
              "surface_area": float(mesh.area), "footprint_area": float(interface_mass.sum())}
    return geometry, normals, footprint_centers, footprint_normals, _dual_boundaries(mesh, nodes), footprint_boundaries, points, report


def _build_edges(pos, centroids, side, chain, backbone, backbone_mask, records, keys,
                 boundaries, footprint_boundaries, ca_cutoff):
    lookup = {key: i for i, key in enumerate(keys)}
    edges = set(boundaries)
    for i, j in cKDTree(pos).query_pairs(ca_cutoff):
        if side[i] == side[j]:
            edges.add((i, j))
    peptide = set()
    for i in range(len(pos) - 1):
        j = i + 1
        if (chain[i] == chain[j] and backbone_mask[i, 2] and backbone_mask[j, 0]
                and np.linalg.norm(backbone[i, 2] - backbone[j, 0]) < 2.):
            peptide.add((i, j))
            edges.add((i, j))
    contacts = {}
    for i, record in enumerate(records):
        for key, count in record["partners"].items():
            pair = tuple(sorted((i, lookup[key])))
            contacts[pair] = count
            edges.add(pair)
    pairs = sorted(edges)
    attrs = [[np.linalg.norm(pos[i] - pos[j]), np.linalg.norm(centroids[i] - centroids[j]),
              side[i] == side[j], (i, j) in contacts, (i, j) in peptide,
              contacts.get((i, j), 0), boundaries.get((i, j), 0), footprint_boundaries.get((i, j), 0)]
             for i, j in pairs]
    undirected = np.asarray(pairs, dtype=np.int64).reshape(-1, 2).T
    attributes = np.asarray(attrs, dtype=np.float32).reshape(-1, len(EDGE_FEATURE_NAMES))
    return np.concatenate([undirected, undirected[::-1]], axis=1), np.concatenate([attributes, attributes], axis=0)


def extract_complex_graph(structure_path, chains_a, chains_b, feature_cache=None,
                          max_surface_points=2048, ca_cutoff=10., curvature_method="angle_defect", *, prepared=None, descriptor_bundle=None) -> dict[str, np.ndarray]:
    """Extract a complete-partner graph; partner swap produces the same graph."""
    group_a, group_b = sorted((_chains(chains_a), _chains(chains_b)))
    prepared = prepared or prepare_interface(structure_path, group_a, group_b)
    bundle = (descriptor_bundle if descriptor_bundle is not None else
              extract_complex_features(structure_path, group_a, group_b, feature_cache, prepared=prepared))
    left, right = prepared.partner(group_a), prepared.partner(group_b)
    keys = list(bundle["residues"])
    records = [bundle["residues"][key] for key in keys]
    lookup = {key: i for i, key in enumerate(keys)}
    count = len(keys)
    centroids = np.asarray([r["centroid"] for r in records], dtype=np.float32)
    backbone = np.repeat(centroids[:, None, :], 4, axis=1)
    backbone_mask = np.zeros((count, 4), dtype=bool)
    side = np.asarray([r["side"] for r in records], dtype=np.int64)
    chain = np.asarray([r["chain"] for r in records])
    atoms_a = left.atoms
    atoms_b = right.atoms
    for atom, label in zip(atoms_a + atoms_b, source_atom_residue_labels(atoms_a + atoms_b), strict=True):
        if atom.get_name() in BACKBONE_ATOMS:
            row, col = lookup[str(label)], BACKBONE_ATOMS.index(atom.get_name())
            backbone[row, col] = atom.get_coord()
            backbone_mask[row, col] = True
    pos = backbone[:, 1].copy()  # Explicit centroid fallback when CA is unresolved.
    chemical = []
    for record in records:
        aa_index = AA.find(record["aa"])
        identity = np.zeros(21)
        identity[aa_index if aa_index >= 0 else 20] = 1.
        chemical.append([*identity, *residue_properties(record["aa"]).values()])
    contact = [[r["partner_distance"], r["contact_degree"], r["contact_atom_pairs"], r["local_neighbor_count"]]
               for r in records]
    surface_base = [[r["surface_area"], r["footprint_area"], r["footprint_area"] / max(r["surface_area"], 1e-12),
                     r["footprint_patch_count"]] for r in records]
    geometry = np.zeros((count, 11))
    normals, footprint_centers, footprint_normals = (np.zeros((count, 3)) for _ in range(3))
    boundaries, footprint_boundaries, point_parts, reports = {}, {}, [], []
    for partner in (left, right):
        geom, normal, centers, foot_normal, boundary, foot_boundary, points, report = _surface_data(
            partner, lookup, count, max_surface_points, curvature_method=curvature_method)
        geometry += geom
        normals += normal
        footprint_centers += centers
        footprint_normals += foot_normal
        boundaries.update(boundary)
        footprint_boundaries.update(foot_boundary)
        point_parts.append(points)
        reports.append(report)
    surface = np.concatenate([surface_base, geometry], axis=1)
    x = np.concatenate([chemical, contact, surface], axis=1).astype(np.float32)
    edge_index, edge_attr = _build_edges(pos, centroids, side, chain, backbone, backbone_mask, records,
                                       keys, boundaries, footprint_boundaries, ca_cutoff)
    footprint_mask = surface[:, 1] > 0
    footprint_offset = np.zeros((count, 3))
    footprint_offset[footprint_mask] = footprint_centers[footprint_mask] - pos[footprint_mask]
    protocol = GRAPH_PROTOCOL
    if curvature_method == "analytic_gaussian":
        from topoppi_affinity.analytic_curvature import GRAPH_PROTOCOL as analytic_protocol

        protocol = analytic_protocol
    result = {
        "x": x, "feature_names": np.asarray(FEATURE_NAMES),
        "surface_feature_mask": np.asarray([name.startswith("surface.") for name in FEATURE_NAMES]),
        "aa_index": np.asarray([AA.find(r["aa"]) if r["aa"] in AA else 20 for r in records], dtype=np.int64),
        "pos": pos, "centroid": centroids, "backbone": backbone,
        "backbone_mask": backbone_mask, "side": side, "chain": chain,
        "residue_key": np.asarray(keys), "residue_number": np.asarray([r["residue_number"] for r in records]),
        "insertion": np.asarray([r["insertion"] for r in records]),
        "edge_index": edge_index, "edge_attr": edge_attr,
        "edge_feature_names": np.asarray(EDGE_FEATURE_NAMES),
        "edge_surface_feature_mask": np.asarray([name.startswith("surface.") for name in EDGE_FEATURE_NAMES]),
        "surface_normal": normals.astype(np.float32), "footprint_normal": footprint_normals.astype(np.float32),
        "footprint_offset": footprint_offset.astype(np.float32), "footprint_mask": footprint_mask,
        "metadata_json": np.asarray(json.dumps({"protocol": protocol, "curvature_method": curvature_method, "source": bundle["metadata"]["source"],
            "surface_reports": reports, "ca_cutoff": ca_cutoff, "max_surface_points_per_side": max_surface_points,
            "geometry_protocol": "fixed_WT_identity_substitution", "coordinates": "angstrom_original_structure",
            "normal_orientation": "native_mesh_winding", "missing_CA": "centroid_with_backbone_mask_false"})),
    }
    for name in point_parts[0]:
        result[name] = np.concatenate([points[name] for points in point_parts]).astype(
            np.int64 if name == "point_node_index" else np.float32)
    return result


def mutant_node_features(graph, mutations):
    """Return joint mutant node features and mask, validating supplied WT identities."""
    x = graph["x"].copy()
    aa_index = graph["aa_index"].copy()
    mask = np.zeros(len(x), dtype=bool)
    lookup = {(str(chain), int(number), str(insertion)): i for i, (chain, number, insertion) in enumerate(
        zip(graph["chain"], graph["residue_number"], graph["insertion"], strict=True))}
    seen = set()
    for mutation in parse_mutations(mutations):
        key = (mutation["chain"], mutation["position"], mutation["insertion"])
        if key in seen:
            raise ValueError(f"Duplicate mutation site: {key}")
        seen.add(key)
        if key not in lookup:
            raise ValueError(f"Mutation residue absent from graph: {key}")
        node = lookup[key]
        if graph["aa_index"][node] != AA.index(mutation["wt"]):
            raise ValueError(f"WT identity mismatch for {key}")
        if mutation["wt"] == mutation["mt"]:
            continue
        mask[node] = True
        aa_index[node] = AA.index(mutation["mt"])
        x[node, :21] = 0.
        x[node, aa_index[node]] = 1.
        x[node, 21:26] = list(residue_properties(mutation["mt"]).values())
    return x, aa_index, mask


def save_graph(graph, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **graph)


def load_graph(path):
    with np.load(path, allow_pickle=False) as graph:
        return {name: graph[name] for name in graph.files}
