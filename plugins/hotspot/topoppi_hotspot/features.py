"""Label-free residue features from TopoPPI's native, uncut interface surface.

No affinity plugin or learned model is used. Gaussian-surface area is not SASA,
and angle-defect curvature is a mesh descriptor, not a free-energy calculation.
All neighbor readouts exclude the focal residue and contain no hotspot labels.
"""
from __future__ import annotations

from collections import defaultdict

import networkx as nx
import numpy as np
import pandas as pd
from Bio.Data.PDBData import protein_letters_3to1_extended
from scipy.spatial import cKDTree

from topoppi import residue_dual_boundaries, source_atom_residue_labels, vertex_area_masses
from topoppi.atlas.footprints import mesh_vertex_residue_labels
from topoppi.interface import chain_group, prepare_interface

AA = "ACDEFGHIKLMNPQRSTVWY"
# Kyte-Doolittle hydropathy scale; fixed amino-acid reference descriptors.
_HYDRO = dict(zip(AA, [1.8, 2.5, -3.5, -3.5, 2.8, -.4, -3.2, 4.5, -3.9, 3.8,
                      1.9, -3.5, -1.6, -3.5, -4.5, -.8, -.7, 4.2, -.9, -1.3], strict=True))
# Zamyatnin residue volumes in cubic angstroms; independent of fitted models.
_VOLUME = dict(zip(AA, [88.6, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7,
                       168.6, 166.7, 162.9, 114.1, 112.7, 143.8, 173.4, 89.,
                       116.1, 140., 227.8, 193.6], strict=True))
CHEMICAL = ["chemical.hydropathy", "chemical.volume", "chemical.charge",
            "chemical.polar", "chemical.aromatic"]
CONTACT = ["contact.partner_distance", "contact.degree", "contact.atom_pairs",
           "contact.neighbor_count", "contact.degree_centrality",
           "contact.betweenness", "contact.pagerank"]
BASE_FEATURES = [f"aa.{a}" for a in AA + "X"] + CHEMICAL + CONTACT
SURFACE_FEATURES = ["surface.area", "surface.footprint_area", "surface.footprint_fraction",
                    "surface.patch_count", "surface.curvature_mean", "surface.curvature_std",
                    "surface.normal_dispersion", "surface.footprint_patch_hhi",
                    "surface.footprint_largest_patch_fraction"]
NEIGHBOR_FEATURES = CHEMICAL + ["contact.degree", "contact.atom_pairs"] + SURFACE_FEATURES
FEATURE_SETS = {"base": BASE_FEATURES, "self": BASE_FEATURES + SURFACE_FEATURES}
for _name in ("spatial", "topology", "shuffled"):
    FEATURE_SETS[_name] = FEATURE_SETS["self"] + [f"{_name}.neighbor_count"] + [
        f"{_name}.mean.{name}" for name in NEIGHBOR_FEATURES]


def neighbor_means(values, edges):
    """Unweighted means over an undirected simple graph; isolates are zero."""
    values = np.asarray(values, dtype=float)
    pairs = np.asarray(edges, dtype=int).reshape(-1, 2)
    pairs = np.unique(np.sort(pairs[pairs[:, 0] != pairs[:, 1]], axis=1), axis=0)
    total = np.zeros_like(values)
    count = np.zeros(len(values), dtype=int)
    if len(pairs):
        left, right = pairs.T
        np.add.at(total, left, values[right])
        np.add.at(total, right, values[left])
        np.add.at(count, left, 1)
        np.add.at(count, right, 1)
    return total / np.maximum(count[:, None], 1), count


def shuffled_neighbor_means(values, edges, side, interface_mask, seed=42):
    """Keep topology and feature marginals; shuffle joint source rows per stratum.

    Strata are partner side and distance-defined interface membership. This is
    a feature-localization control, not a label permutation or a rewired graph.
    """
    values = np.asarray(values, dtype=float)
    side, interface_mask = np.asarray(side), np.asarray(interface_mask, dtype=bool)
    source = values.copy()
    rng = np.random.default_rng(seed)
    for partner in sorted(np.unique(side)):
        for interface in (False, True):
            indices = np.flatnonzero((side == partner) & (interface_mask == interface))
            source[indices] = values[rng.permutation(indices)]
    return neighbor_means(source, edges)


def _surface_features(partner, lookup, n):
    """Integrate residue corner masses and native barycentric-dual boundaries."""
    labels = source_atom_residue_labels(partner.atoms)
    mesh = partner.mesh
    nodes = np.asarray([lookup[str(key)] for key in mesh_vertex_residue_labels(mesh, labels)])
    mass = vertex_area_masses(mesh)
    area = np.bincount(nodes, weights=mass, minlength=n)
    curvature = np.divide(mesh.vertex_defects, mass, out=np.zeros_like(mass), where=mass > 0)
    mean = np.bincount(nodes, weights=mass * curvature, minlength=n) / np.maximum(area, 1e-12)
    variance = np.bincount(nodes, weights=mass * (curvature - mean[nodes]) ** 2,
                          minlength=n) / np.maximum(area, 1e-12)
    normals = np.stack([np.bincount(nodes, weights=mass * mesh.vertex_normals[:, axis],
                                    minlength=n) for axis in range(3)], axis=1)
    dispersion = np.where(area > 0, 1 - np.linalg.norm(normals, axis=1) / np.maximum(area, 1e-12), 0.)
    patch_masses, boundaries = [], defaultdict(float)
    for patch in partner.patches:
        patch_nodes = np.asarray([lookup[str(key)] for key in mesh_vertex_residue_labels(patch, labels)])
        patch_masses.append(np.bincount(patch_nodes, weights=vertex_area_masses(patch), minlength=n))
        for pair, length in residue_dual_boundaries(patch, patch_nodes).items():
            boundaries[pair] += length
    patches = np.asarray(patch_masses).reshape(-1, n)
    footprint = patches.sum(axis=0)
    fractions = patches / np.maximum(footprint, 1e-12)
    data = {
        "surface.area": area, "surface.footprint_area": footprint,
        "surface.footprint_fraction": footprint / np.maximum(area, 1e-12),
        "surface.patch_count": (patches > 0).sum(axis=0),
        "surface.curvature_mean": mean, "surface.curvature_std": np.sqrt(variance),
        "surface.normal_dispersion": np.clip(dispersion, 0., 1.),
        "surface.footprint_patch_hhi": (fractions ** 2).sum(axis=0),
        "surface.footprint_largest_patch_fraction": fractions.max(axis=0) if len(patches) else np.zeros(n),
    }
    return data, boundaries


def extract_features(structure_path, chains_a, chains_b, *, prepared=None,
                     contact_cutoff=6., neighbor_cutoff=10., shuffle_seed=42):
    """Return all partner residues and an index-aligned organization graph.

    Canonical partner ordering makes the output invariant to partner swapping.
    Residue numbers and insertion codes are those in the supplied structure.
    ``prepared`` optionally reuses a core PreparedInterface with both surfaces.
    Graph edges are unique undirected integer pairs of shape (E, 2). All graph
    arrays may be saved in NPZ; ``metadata`` is a JSON-compatible dictionary.
    Spatial and native neighbor pools contain only same-side interface nodes.
    """
    groups = sorted((chain_group(chains_a), chain_group(chains_b)))
    prepared = prepared or prepare_interface(structure_path, *groups)
    partners = [prepared.partner(group) for group in groups]
    rows = []
    for side, partner in enumerate(partners):
        atom_labels = source_atom_residue_labels(partner.atoms)
        distances = cKDTree(partners[1 - side].coords).query(partner.coords)[0]
        for key in dict.fromkeys(atom_labels):
            indices = np.flatnonzero(atom_labels == key)
            residue = partner.atoms[int(indices[0])].get_parent()
            wt = protein_letters_3to1_extended.get(residue.resname.upper(), "X")
            position = residue["CA"].coord if "CA" in residue else partner.coords[indices].mean(axis=0)
            rows.append({"chain": str(residue.parent.id), "residue_number": int(residue.id[1]),
                         "insertion_code": str(residue.id[2]).strip(), "residue_key": str(key),
                         "wt": wt, "side": side, "x": float(position[0]), "y": float(position[1]),
                         "z": float(position[2]), "ca_resolved": "CA" in residue,
                         **{f"aa.{a}": float(wt == a) for a in AA + "X"},
                         "chemical.hydropathy": _HYDRO.get(wt, 0.),
                         "chemical.volume": _VOLUME.get(wt, 0.),
                         "chemical.charge": float((wt in "KR") - (wt in "DE")),
                         "chemical.polar": float(wt in "CDEHKNQRSTY"),
                         "chemical.aromatic": float(wt in "FWY"),
                         "contact.partner_distance": float(distances[indices].min())})
    frame = pd.DataFrame(rows)
    n = len(frame)
    lookup = {key: i for i, key in enumerate(frame.residue_key)}
    side = frame.side.to_numpy()
    positions = frame[["x", "y", "z"]].to_numpy()
    canonical_first = next(i for i, partner in enumerate(prepared.partners) if partner is partners[0])
    contact_map = prepared.contacts(side=canonical_first, cutoff=contact_cutoff)
    cross = {tuple(sorted((lookup[a], lookup[b]))): count
             for a, neighbors in contact_map.items() for b, count in neighbors.items()}
    cross_edges = np.asarray(sorted(cross), dtype=int).reshape(-1, 2)
    degree, atom_pairs = np.zeros(n), np.zeros(n)
    for (i, j), weight in cross.items():
        degree[[i, j]] += 1
        atom_pairs[[i, j]] += weight
    frame["contact.degree"], frame["contact.atom_pairs"] = degree, atom_pairs
    contact_graph = nx.Graph()
    contact_graph.add_edges_from(cross_edges.tolist())
    for name, scores in [("degree_centrality", nx.degree_centrality(contact_graph)),
                         ("betweenness", nx.betweenness_centrality(contact_graph)),
                         ("pagerank", nx.pagerank(contact_graph) if cross else {})]:
        frame[f"contact.{name}"] = [scores.get(i, 0.) for i in range(n)]
    spatial_all = np.asarray(sorted((i, j) for i, j in cKDTree(positions).query_pairs(neighbor_cutoff)
                                   if side[i] == side[j]), dtype=int).reshape(-1, 2)
    _, all_counts = neighbor_means(np.zeros((n, 1)), spatial_all)
    frame["contact.neighbor_count"] = all_counts
    surface = {name: np.zeros(n) for name in SURFACE_FEATURES}
    boundaries = {}
    for partner in partners:
        values, boundary = _surface_features(partner, lookup, n)
        for name in surface:
            surface[name] += values[name]
        boundaries.update(boundary)
    frame = pd.concat([frame, pd.DataFrame(surface)], axis=1)
    interface = frame["contact.partner_distance"].to_numpy() <= contact_cutoff
    # Both neighborhood controls use the same nodes; only adjacency changes.
    spatial = spatial_all[interface[spatial_all].all(axis=1)]
    native = np.asarray(sorted(pair for pair in boundaries if interface[list(pair)].all()), dtype=int).reshape(-1, 2)
    values = frame[NEIGHBOR_FEATURES].to_numpy()
    readouts = {"spatial": neighbor_means(values, spatial), "topology": neighbor_means(values, native),
                "shuffled": shuffled_neighbor_means(values, native, side, interface, shuffle_seed)}
    extra = {}
    for name, (means, counts) in readouts.items():
        extra[f"{name}.neighbor_count"] = counts
        extra.update({f"{name}.mean.{key}": means[:, i] for i, key in enumerate(NEIGHBOR_FEATURES)})
    frame = pd.concat([frame, pd.DataFrame(extra)], axis=1)
    footprint_edges = np.asarray(sorted(boundaries), dtype=int).reshape(-1, 2)
    graph = {"residue_keys": frame.residue_key.to_numpy(dtype=str), "positions": positions, "side": side,
             "ca_resolved": frame.ca_resolved.to_numpy(dtype=bool),
             "cross_edges": cross_edges, "cross_weights": np.asarray([cross[tuple(e)] for e in cross_edges]),
             "spatial_edges": spatial, "footprint_edges": footprint_edges,
             "footprint_weights": np.asarray([boundaries[tuple(e)] for e in footprint_edges]),
             "interface_mask": interface, "native_neighbor_edges": native,
             "metadata": {"protocol": "hotspot_native_v1", "contact_cutoff_angstrom": contact_cutoff,
                          "spatial_neighbor_cutoff_angstrom": neighbor_cutoff,
                          "curvature": "angle_defect_divided_by_barycentric_vertex_area; area-weighted residue moments",
                          "surface": "TopoPPI Gaussian isosurface; not solvent-accessible surface area",
                          "footprint_concentration": "distribution of a residue's physical area over retained native patches",
                          "centrality": "unweighted bipartite contact graph; contacting nodes only",
                          "neighbor_pool": "same-side residues with partner atom distance <= contact_cutoff; no focal self",
                          "shuffle": "joint source-row permutation within partner/interface strata", "shuffle_seed": shuffle_seed,
                          "surface_area_by_side": [float(p.mesh.area) for p in partners],
                          "footprint_area_by_side": [float(sum(patch.area for patch in p.patches)) for p in partners],
                          "groups": [list(group) for group in groups]}}
    return frame, graph
