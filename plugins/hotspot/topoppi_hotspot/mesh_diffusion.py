"""Label-free implicit diffusion on native TopoPPI triangle surfaces.

This is an adapted fixed mesh-operator baseline, not a DiffusionNet reproduction
or a trained neural network. Physical coordinates remain in angstroms. Cotangent
FEM stiffness and barycentric mass define (M + t L) u_t = M f on the complete
selected partner surface. Native physical footprint areas pool vertex channels
back to residues, without a UV layout or a solvent-accessible-area interpretation.
"""
from __future__ import annotations

import time

import numpy as np
import pandas as pd
from Bio.Data.PDBData import protein_letters_3to1_extended
from scipy import sparse
from scipy.sparse.linalg import splu
from scipy.spatial import cKDTree

from topoppi.atlas.footprints import mesh_vertex_residue_labels, source_atom_residue_labels
from topoppi.interface import chain_group, prepare_interface

from .features import _HYDRO, _VOLUME, CHEMICAL

HEAT_TIMES = (1., 4., 16.)
CHANNELS = tuple(CHEMICAL) + ("partner_distance", "gaussian_curvature", "mean_curvature_magnitude", "normal_dispersion")


def cotangent_operators(vertices, faces):
    """Return cotangent FEM stiffness L, lumped barycentric mass and diagnostics.

    Cotangents retain their signs on obtuse triangles: no graph-weight clipping,
    artificial diagonal regularization or scale normalization is performed.
    Exactly zero-area triangles contribute no mass/stiffness. Zero-mass vertices
    retain their original signal in the resolvent and cannot enter area pooling.
    """
    vertices, faces = np.asarray(vertices, float), np.asarray(faces, int)
    triangle = vertices[faces]
    twice_area = np.linalg.norm(np.cross(triangle[:, 1] - triangle[:, 0], triangle[:, 2] - triangle[:, 0]), axis=1)
    valid = twice_area > 0
    use, triangle, twice_area = faces[valid], triangle[valid], twice_area[valid]
    mass = np.bincount(use.ravel(), weights=np.repeat(twice_area / 6., 3), minlength=len(vertices))
    rows, columns, weights, negative = [], [], [], 0
    for corner in range(3):
        left, right = (corner + 1) % 3, (corner + 2) % 3
        u, v = triangle[:, left] - triangle[:, corner], triangle[:, right] - triangle[:, corner]
        weight = .5 * np.einsum("ij,ij->i", u, v) / twice_area
        negative += int((weight < 0).sum())
        a, b = use[:, left], use[:, right]
        rows.extend([a, b, a, b])
        columns.extend([a, b, b, a])
        weights.extend([weight, weight, -weight, -weight])
    stiffness = sparse.coo_matrix((np.concatenate(weights), (np.concatenate(rows), np.concatenate(columns))),
                                  shape=(len(vertices), len(vertices))).tocsr()
    stiffness.sum_duplicates()
    diagnostics = {"n_vertices": len(vertices), "n_faces": len(faces), "n_zero_area_faces": int((~valid).sum()),
                   "n_zero_mass_vertices": int((mass == 0).sum()), "n_negative_corner_cotangents": negative,
                   "area_angstrom2": float(mass.sum()), "stiffness_nnz": int(stiffness.nnz),
                   "max_abs_constant_residual": float(np.max(np.abs(stiffness @ np.ones(len(vertices)))))}
    return stiffness, mass, diagnostics


def implicit_resolvent(stiffness, mass, values, times=HEAT_TIMES):
    """Sparse implicit heat step at each fixed physical time, from the same f."""
    mass, values = np.asarray(mass, float), np.asarray(values, float)
    active = np.flatnonzero(mass > 0)
    local_stiffness = stiffness[active][:, active].tocsc()
    local_mass = sparse.diags(mass[active], format="csc")
    rhs = mass[active, None] * values[active]
    result = {0.: values.copy()}
    for heat_time in times:
        heat_time = float(heat_time)
        if heat_time <= 0 or not np.isfinite(heat_time):
            raise ValueError("Heat times must be finite positive physical scales.")
        solved = values.copy()
        solved[active] = splu(local_mass + heat_time * local_stiffness).solve(rhs)
        if not np.isfinite(solved).all():
            raise ValueError("Mesh diffusion produced a nonfinite feature.")
        result[heat_time] = solved
    return result


def geometry_channels(vertices, faces, stiffness, mass):
    """Rotation/translation invariant curvature and incident-normal variation."""
    vertices, faces = np.asarray(vertices, float), np.asarray(faces, int)
    triangles = vertices[faces]
    cross = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    twice_area = np.linalg.norm(cross, axis=1)
    valid = twice_area > 0
    triangles, use, cross, twice_area = triangles[valid], faces[valid], cross[valid], twice_area[valid]
    angles = np.zeros(len(vertices))
    for corner in range(3):
        u = triangles[:, (corner + 1) % 3] - triangles[:, corner]
        v = triangles[:, (corner + 2) % 3] - triangles[:, corner]
        np.add.at(angles, use[:, corner], np.arctan2(twice_area, np.einsum("ij,ij->i", u, v)))
    edges = np.sort(np.concatenate([use[:, [0, 1]], use[:, [1, 2]], use[:, [2, 0]]]), axis=1)
    unique_edges, counts = np.unique(edges, axis=0, return_counts=True)
    reference = np.full(len(vertices), 2 * np.pi)
    reference[np.unique(unique_edges[counts == 1])] = np.pi
    gaussian = np.divide(reference - angles, mass, out=np.zeros(len(vertices)), where=mass > 0)
    laplacian_position = np.divide(stiffness @ vertices, mass[:, None], out=np.zeros_like(vertices), where=mass[:, None] > 0)
    mean_magnitude = .5 * np.linalg.norm(laplacian_position, axis=1)
    normal_sum = np.zeros_like(vertices)
    for corner in range(3):
        np.add.at(normal_sum, use[:, corner], cross / 6.)
    normal_dispersion = 1 - np.divide(np.linalg.norm(normal_sum, axis=1), mass, out=np.ones(len(vertices)), where=mass > 0)
    return np.column_stack([gaussian, mean_magnitude, np.clip(normal_dispersion, 0., 1.)])


def native_footprint_pooling(mesh, patches, source_labels, residue_keys):
    """Sparse residue×parent-vertex area operator from native patch provenance.

    Each retained triangle contributes one third of its physical area to each
    labeled corner (TopoPPI's barycentric footprint convention). Source IDs map
    sanitized/cut duplicate vertices back to the original full-partner surface.
    No nearest-coordinate remapping or UV integration is used.
    """
    vertex_ids = np.asarray(mesh.metadata["source_vertex_ids"], int)
    vertex_lookup = {int(source): i for i, source in enumerate(vertex_ids)}
    residue_lookup = {str(key): i for i, key in enumerate(residue_keys)}
    labels = mesh_vertex_residue_labels(mesh, source_labels)
    rows, columns, data, seen_faces = [], [], [], set()
    for patch in patches:
        source_ids = np.asarray(patch.metadata["source_vertex_ids"], int)
        parent = np.asarray([vertex_lookup[int(source)] for source in source_ids])
        if not np.allclose(np.asarray(patch.vertices), np.asarray(mesh.vertices)[parent], atol=1e-8, rtol=0):
            raise ValueError("Native patch source IDs no longer preserve parent 3D coordinates.")
        face_ids = list(map(int, patch.metadata["source_face_ids"]))
        if seen_faces.intersection(face_ids):
            raise ValueError("Native patches repeat a source face and would double-count footprint area.")
        seen_faces.update(face_ids)
        corners = parent[np.asarray(patch.faces)].ravel()
        rows.append(np.asarray([residue_lookup[str(label)] for label in labels[corners]]))
        columns.append(corners)
        data.append(np.repeat(np.asarray(patch.area_faces) / 3., 3))
    if not data:
        return sparse.csr_matrix((len(residue_keys), len(mesh.vertices)))
    return sparse.coo_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(columns))),
                             shape=(len(residue_keys), len(mesh.vertices))).tocsr()


def pool_channels(pooling, values):
    """Area-weighted mean and SD; no-footprint residues get zero plus a flag."""
    area = np.asarray(pooling.sum(axis=1)).ravel()
    values = np.asarray(values, float)
    mean = np.divide(pooling @ values, area[:, None], out=np.zeros((len(area), values.shape[1])), where=area[:, None] > 0)
    entries = pooling.tocoo()
    centered = values[entries.col] - mean[entries.row]
    # Chemical fields at t=0 are constant inside each residue footprint. Avoid
    # turning floating-point cancellation into a spurious standardized feature.
    tolerance = 64 * np.finfo(float).eps * np.maximum(1., np.max(np.abs(values), axis=0))
    centered[np.abs(centered) <= tolerance] = 0.
    variance = np.zeros_like(mean)
    np.add.at(variance, entries.row, entries.data[:, None] * centered ** 2)
    variance = np.divide(variance, area[:, None], out=variance, where=area[:, None] > 0)
    return mean, np.sqrt(variance), area


def extract_mesh_diffusion(structure_path, chains_a, chains_b, *, prepared=None):
    """Return all partner residue rows and label-free extraction provenance."""
    groups = sorted((chain_group(chains_a), chain_group(chains_b)))
    prepared = prepared or prepare_interface(structure_path, *groups)
    partners = [prepared.partner(group) for group in groups]
    parts, diagnostics = [], []
    for side, partner in enumerate(partners):
        started = time.monotonic()
        mesh, other = partner.mesh, partners[1 - side]
        atom_labels = source_atom_residue_labels(partner.atoms)
        residue_keys = list(dict.fromkeys(map(str, atom_labels)))
        distances = cKDTree(other.coords).query(partner.coords)[0]
        residues, chemical = [], {}
        for key in residue_keys:
            positions = np.flatnonzero(atom_labels == key)
            residue = partner.atoms[int(positions[0])].get_parent()
            wt = protein_letters_3to1_extended.get(residue.resname.upper(), "X")
            chemical[key] = [_HYDRO.get(wt, 0.), _VOLUME.get(wt, 0.), float((wt in "KR") - (wt in "DE")),
                             float(wt in "CDEHKNQRSTY"), float(wt in "FWY")]
            residues.append({"residue_key": key, "side": side, "wt": wt, "resname": residue.resname,
                             "contact.partner_distance": float(distances[positions].min())})
        vertex_labels = mesh_vertex_residue_labels(mesh, atom_labels)
        stiffness, mass, report = cotangent_operators(mesh.vertices, mesh.faces)
        channels = np.column_stack([np.asarray([chemical[str(label)] for label in vertex_labels]),
                                    cKDTree(other.coords).query(mesh.vertices)[0],
                                    geometry_channels(mesh.vertices, mesh.faces, stiffness, mass)])
        diffused = implicit_resolvent(stiffness, mass, channels)
        pooling = native_footprint_pooling(mesh, partner.patches, atom_labels, residue_keys)
        extra = {}
        for heat_time, values in diffused.items():
            means, standard_deviation, area = pool_channels(pooling, values)
            for i, channel in enumerate(CHANNELS):
                extra[f"mesh_diffusion.t{heat_time:g}.mean.{channel}"] = means[:, i]
                extra[f"mesh_diffusion.t{heat_time:g}.std.{channel}"] = standard_deviation[:, i]
        # Modulus before pooling retains opposing responses within a footprint.
        for first, second in ((1., 4.), (4., 16.)):
            wavelet, _, _ = pool_channels(pooling, np.abs(diffused[first] - diffused[second]))
            for i, channel in enumerate(CHANNELS):
                extra[f"mesh_wavelet.t{first:g}_t{second:g}.absmean.{channel}"] = wavelet[:, i]
        residue_lookup = {key: i for i, key in enumerate(residue_keys)}
        total_area = np.bincount(np.asarray([residue_lookup[str(label)] for label in vertex_labels]), weights=mass, minlength=len(residue_keys))
        extra.update({"mesh_coverage.footprint_area": area,
                      "mesh_coverage.full_residue_surface_area": total_area,
                      "mesh_coverage.footprint_fraction": np.divide(area, total_area, out=np.zeros_like(area), where=total_area > 0),
                      "mesh_coverage.missing_footprint": (area == 0).astype(int),
                      "mesh_coverage.missing_full_surface": (total_area == 0).astype(int)})
        frame = pd.concat([pd.DataFrame(residues), pd.DataFrame(extra)], axis=1)
        parts.append(frame)
        diagnostics.append({**report, "side": side, "chains": list(partner.chains), "n_residues": len(frame),
                            "n_residues_without_footprint": int((area == 0).sum()), "native_patch_area": float(area.sum()),
                            "n_patches": len(partner.patches), "feature_seconds": time.monotonic() - started,
                            "mass_conservation_max_abs": max(float(np.max(np.abs(mass @ values - mass @ channels))) for values in diffused.values())})
    result = pd.concat(parts, ignore_index=True)
    return result, {"protocol": "hotspot_native_mesh_resolvent_v1", "prepared_interface": prepared.metadata(),
                    "heat_times_angstrom2": list(HEAT_TIMES), "source_channels": list(CHANNELS),
                    "operator": "cotangent FEM stiffness; barycentric mass; u=(M+tL)^-1 M f; full partner mesh",
                    "pooling": "native footprint physical face-corner area/3, mapped by source vertex provenance",
                    "wavelet": "area-weighted mean of vertexwise |u1-u4| and |u4-u16|; modulus before footprint pooling",
                    "missing_footprint": "zero-filled pooled channels with explicit missing flag; no imputed area",
                    "interpretation": "fixed implicit mesh diffusion/scattering-inspired adaptation; not original DiffusionNet; not SASA",
                    "sources": ["https://arxiv.org/abs/2012.00888", "https://github.com/nmwsharp/diffusion-net"],
                    "labels_used": False, "partners": diagnostics,
                    "feature_blocks": {"diffusion": [c for c in result if c.startswith("mesh_diffusion.")],
                                       "wavelet": [c for c in result if c.startswith("mesh_wavelet.")],
                                       "coverage": [c for c in result if c.startswith("mesh_coverage.")]}}
