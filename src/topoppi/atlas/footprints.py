"""Residue-footprint connectivity metrics for cut UV atlases."""

from __future__ import annotations

from collections import defaultdict
from typing import Iterable, Mapping, Sequence

import numpy as np
import trimesh

from topoppi.atlas.seams import UVSeamTopology, uv_seam_topology
from topoppi.atlas.uv import as_corner_uv
from topoppi.contacts import geometric_contact_partner_map as geometric_contact_partner_map
from topoppi.mesh.provenance import (
    OPTCUTS_GEOMETRY_VERTEX_IDS,
    SOURCE_ATOM_IDS,
    SOURCE_VERTEX_IDS,
)
from topoppi.residues import atom_residue_label as atom_residue_label
from topoppi.residues import source_atom_residue_labels as source_atom_residue_labels

_EPS = 1e-12


class _UnionFind:
    def __init__(self, size: int):
        self.parent = np.arange(size, dtype=np.int64)
        self.rank = np.zeros(size, dtype=np.int8)

    def find(self, item: int) -> int:
        parent = self.parent
        root = int(item)
        while parent[root] != root:
            root = int(parent[root])
        while parent[item] != item:
            next_item = int(parent[item])
            parent[item] = root
            item = next_item
        return root

    def union(self, left: int, right: int) -> None:
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left == root_right:
            return
        rank = self.rank
        if rank[root_left] < rank[root_right]:
            root_left, root_right = root_right, root_left
        self.parent[root_right] = root_left
        if rank[root_left] == rank[root_right]:
            rank[root_left] += 1

    def roots(self) -> np.ndarray:
        return np.fromiter((self.find(index) for index in range(len(self.parent))), dtype=np.int64)






def mesh_vertex_residue_labels(
    mesh: trimesh.Trimesh,
    source_labels: Sequence[str],
) -> np.ndarray:
    """Map source-atom residue labels onto current mesh vertices."""

    if SOURCE_ATOM_IDS not in mesh.metadata:
        raise ValueError("Mesh metadata is missing source_atom_indices.")
    source_atom_ids = np.asarray(mesh.metadata[SOURCE_ATOM_IDS], dtype=np.int64)
    if source_atom_ids.shape != (len(mesh.vertices),):
        raise ValueError("source_atom_indices must contain one index per mesh vertex.")
    labels = np.asarray(source_labels, dtype=object)
    if np.any((source_atom_ids < 0) | (source_atom_ids >= len(labels))):
        raise ValueError("source_atom_indices contains an index outside the source atom array.")
    return labels[source_atom_ids]


def contact_partner_counts(contact_pairs: Iterable[tuple[str, str]]) -> dict[str, float]:
    """Count distinct partner residues for each residue on the mapped chain."""

    partners: dict[str, set[str]] = defaultdict(set)
    for residue_a, residue_b in contact_pairs:
        partners[str(residue_a)].add(str(residue_b))
    return {residue: float(len(values)) for residue, values in partners.items()}




def contact_partner_degrees(
    partner_map: Mapping[str, Mapping[str, int]],
) -> dict[str, float]:
    """Return the distinct partner-residue degree encoded by a contact map."""

    return {str(residue): float(len(partners)) for residue, partners in partner_map.items() if partners}


def residue_aware_residue_weights(
    source_labels: Sequence[str],
    interaction_counts: Mapping[str, float],
) -> dict[str, float]:
    """Give every mapped-chain residue unit weight plus its interaction degree."""

    weights: dict[str, float] = {}
    for value in np.asarray(source_labels, dtype=object):
        label = str(value)
        interaction_weight = float(interaction_counts.get(label, 0.0))
        if not np.isfinite(interaction_weight) or interaction_weight < 0.0:
            raise ValueError("Residue interaction weights must be finite and non-negative.")
        weights[label] = 1.0 + interaction_weight
    return weights


def _face_residue_masses(
    mesh: trimesh.Trimesh,
    vertex_labels: np.ndarray,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Integrate piecewise-linear corner indicators over each triangle."""

    face_ids: dict[str, list[int]] = defaultdict(list)
    masses: dict[str, list[float]] = defaultdict(list)
    areas = np.asarray(mesh.area_faces, dtype=np.float64)
    for face_index, face in enumerate(np.asarray(mesh.faces, dtype=np.int64)):
        counts: dict[str, int] = defaultdict(int)
        for vertex_index in face:
            counts[str(vertex_labels[int(vertex_index)])] += 1
        for label, count in counts.items():
            face_ids[label].append(face_index)
            masses[label].append(float(areas[face_index]) * count / 3.0)
    return {
        label: (
            np.asarray(face_ids[label], dtype=np.int64),
            np.asarray(masses[label], dtype=np.float64),
        )
        for label in sorted(face_ids)
    }


def _residue_connectivity(
    mesh: trimesh.Trimesh,
    seam_topology: UVSeamTopology,
    vertex_labels: np.ndarray,
    label: str,
    face_ids: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return baseline and post-seam face-dual connectivity for one residue."""

    internal_edges = np.flatnonzero(seam_topology.internal_mask)
    edge_vertices = seam_topology.edges[internal_edges]
    relevant = np.any(vertex_labels[edge_vertices] == label, axis=1)
    footprint_edges = internal_edges[relevant]

    global_to_local = np.full(len(mesh.faces), -1, dtype=np.int64)
    global_to_local[face_ids] = np.arange(len(face_ids), dtype=np.int64)
    face_pairs = seam_topology.incident_faces[footprint_edges]
    local_pairs = global_to_local[face_pairs]
    keep = np.all(local_pairs >= 0, axis=1)
    footprint_edges = footprint_edges[keep]
    local_pairs = local_pairs[keep]

    baseline = _UnionFind(len(face_ids))
    cut = _UnionFind(len(face_ids))
    seam_flags = seam_topology.seam_mask[footprint_edges]
    for pair, is_seam in zip(local_pairs, seam_flags, strict=True):
        left, right = (int(pair[0]), int(pair[1]))
        baseline.union(left, right)
        if not is_seam:
            cut.union(left, right)
    return footprint_edges, local_pairs, seam_flags, baseline.roots(), cut.roots()


def _residue_record(
    mesh: trimesh.Trimesh,
    seam_topology: UVSeamTopology,
    vertex_labels: np.ndarray,
    label: str,
    face_ids: np.ndarray,
    face_masses: np.ndarray,
) -> dict[str, object]:
    footprint_edges, local_pairs, seam_flags, baseline_roots, cut_roots = _residue_connectivity(
        mesh,
        seam_topology,
        vertex_labels,
        label,
        face_ids,
    )
    baseline_component_count = int(len(np.unique(baseline_roots)))
    cut_component_count = int(len(np.unique(cut_roots)))
    extra_components = cut_component_count - baseline_component_count
    cycle_rank = int(len(local_pairs) - len(face_ids) + baseline_component_count)

    baseline_mass: dict[int, float] = defaultdict(float)
    piece_mass: dict[tuple[int, int], float] = defaultdict(float)
    for index, mass in enumerate(face_masses):
        baseline_root = int(baseline_roots[index])
        cut_root = int(cut_roots[index])
        baseline_mass[baseline_root] += float(mass)
        piece_mass[(baseline_root, cut_root)] += float(mass)

    fragmentation_mass = 0.0
    squared_piece_mass: dict[int, float] = defaultdict(float)
    for (baseline_root, _piece_root), mass in piece_mass.items():
        squared_piece_mass[baseline_root] += mass * mass
    for baseline_root, total_mass in baseline_mass.items():
        if total_mass <= _EPS:
            continue
        fragmentation_mass += total_mass - squared_piece_mass[baseline_root] / total_mass

    total_mass = float(np.sum(face_masses))
    fragmentation_mass = min(max(float(fragmentation_mass), 0.0), total_mass)
    seam_crossings = int(np.count_nonzero(seam_flags))
    seam_length = float(np.sum(seam_topology.edge_lengths_3d[footprint_edges[seam_flags]]))
    return {
        "residue": label,
        "footprint_area": total_mass,
        "face_count": int(len(face_ids)),
        "dual_edge_count": int(len(local_pairs)),
        "baseline_component_count": baseline_component_count,
        "component_count_after_seams": cut_component_count,
        "extra_component_count": int(extra_components),
        "cycle_rank": max(cycle_rank, 0),
        "seam_crossing_edge_count": seam_crossings,
        "nonseparating_seam_crossing_edge_count": max(seam_crossings - extra_components, 0),
        "seam_crossing_length_3d": seam_length,
        "fragmentation_mass": fragmentation_mass,
    }


def residue_footprint_pieces(
    mesh: trimesh.Trimesh,
    uv: np.ndarray,
    vertex_labels: Sequence[str],
    *,
    atol: float = 1e-9,
) -> dict[str, list[dict[str, object]]]:
    """Locate each post-seam residue-footprint component in UV space.

    Component mass is the same integrated piecewise-linear corner indicator
    used by the TopoPPI objective.  Its UV centroid is the exact first
    moment of that indicator on each affine UV triangle, aggregated with the
    original 3-D face-area mass convention.
    """

    labels = np.asarray(vertex_labels, dtype=object)
    if labels.shape != (len(mesh.vertices),):
        raise ValueError("vertex_labels must contain one residue label per mesh vertex.")
    corners = as_corner_uv(mesh, uv)
    if not np.isfinite(corners).all():
        raise ValueError("Residue-footprint centroids require finite UV coordinates.")
    seams = uv_seam_topology(mesh, corners, atol=atol)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    areas = np.asarray(mesh.area_faces, dtype=np.float64)
    result: dict[str, list[dict[str, object]]] = {}

    for label, (face_ids, face_masses) in _face_residue_masses(mesh, labels).items():
        _edges, _pairs, _flags, _baseline_roots, cut_roots = _residue_connectivity(
            mesh,
            seams,
            labels,
            label,
            face_ids,
        )
        mass_by_root: dict[int, float] = defaultdict(float)
        moment_by_root: dict[int, np.ndarray] = defaultdict(lambda: np.zeros(2, dtype=np.float64))
        faces_by_root: dict[int, list[int]] = defaultdict(list)
        for local_index, face_index in enumerate(face_ids):
            face_index = int(face_index)
            root = int(cut_roots[local_index])
            triangle = corners[face_index]
            triangle_sum = np.sum(triangle, axis=0)
            labelled_corners = np.flatnonzero(labels[faces[face_index]] == label)
            hat_mass = float(areas[face_index]) / 3.0
            for corner_index in labelled_corners:
                hat_centroid = (triangle_sum + triangle[int(corner_index)]) / 4.0
                moment_by_root[root] += hat_mass * hat_centroid
            mass_by_root[root] += float(face_masses[local_index])
            faces_by_root[root].append(face_index)

        total_mass = float(np.sum(face_masses))
        pieces = []
        for root, mass in mass_by_root.items():
            if mass <= _EPS:
                continue
            centroid = moment_by_root[root] / mass
            pieces.append(
                {
                    "footprint_mass": float(mass),
                    "footprint_mass_fraction": float(mass / total_mass) if total_mass > _EPS else 0.0,
                    "uv_centroid": [float(centroid[0]), float(centroid[1])],
                    "face_count": int(len(faces_by_root[root])),
                    "minimum_face_index": int(min(faces_by_root[root])),
                }
            )
        pieces.sort(
            key=lambda piece: (
                -float(piece["footprint_mass"]),
                float(piece["uv_centroid"][0]),
                float(piece["uv_centroid"][1]),
                int(piece["minimum_face_index"]),
            )
        )
        result[label] = pieces
    return result


def analyze_residue_footprints(
    mesh: trimesh.Trimesh,
    uv: np.ndarray,
    vertex_labels: Sequence[str],
    *,
    atol: float = 1e-9,
) -> tuple[list[dict[str, object]], UVSeamTopology]:
    """Analyze seam-induced fragmentation of every residue footprint."""

    labels = np.asarray(vertex_labels, dtype=object)
    if labels.shape != (len(mesh.vertices),):
        raise ValueError("vertex_labels must contain one residue label per mesh vertex.")
    seams = uv_seam_topology(mesh, uv, atol=atol)
    records = [
        _residue_record(mesh, seams, labels, label, face_ids, masses)
        for label, (face_ids, masses) in _face_residue_masses(mesh, labels).items()
    ]
    return records, seams


def write_residue_footprint_sidecar(
    mesh: trimesh.Trimesh,
    uv: np.ndarray,
    source_labels: Sequence[str],
    path: str,
    *,
    residue_weights: Mapping[str, float] | None = None,
    input_source_vertices: Sequence[int] | None = None,
    atol: float = 1e-9,
) -> dict[str, object]:
    """Write the deterministic numeric sidecar consumed by the C++ energy core."""

    vertex_labels = mesh_vertex_residue_labels(mesh, source_labels)
    face_masses = _face_residue_masses(mesh, vertex_labels)
    labels = sorted(face_masses)
    label_ids = {label: index for index, label in enumerate(labels)}
    weights = residue_weights or {}
    objective_weights = [float(weights.get(label, 1.0)) for label in labels]
    if any(not np.isfinite(weight) or weight < 0.0 for weight in objective_weights):
        raise ValueError("Residue footprint weights must be finite and non-negative.")

    face_entries: list[list[tuple[int, float]]] = [[] for _ in range(len(mesh.faces))]
    for label, (face_ids, masses) in face_masses.items():
        label_id = label_ids[label]
        for face, mass in zip(face_ids, masses, strict=True):
            face_entries[int(face)].append((label_id, float(mass)))
    for entries in face_entries:
        entries.sort()

    seams = uv_seam_topology(mesh, uv, atol=atol)
    internal_edges = np.flatnonzero(seams.internal_mask)
    if input_source_vertices is None:
        identity_key = (
            OPTCUTS_GEOMETRY_VERTEX_IDS if OPTCUTS_GEOMETRY_VERTEX_IDS in mesh.metadata else SOURCE_VERTEX_IDS
        )
        source_vertices = np.asarray(
            mesh.metadata.get(identity_key, np.arange(len(mesh.vertices))),
            dtype=np.int64,
        )
        if source_vertices.shape != (len(mesh.vertices),):
            raise ValueError(f"{identity_key} must contain one ID per mesh vertex.")
    else:
        source_vertices = np.asarray(input_source_vertices, dtype=np.int64)
        if source_vertices.ndim != 1 or len(source_vertices) == 0:
            raise ValueError("input_source_vertices must be a non-empty one-dimensional sequence.")
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("TOPOPPI_FOOTPRINT_V2\n")
        handle.write(f"COUNTS {len(mesh.faces)} {len(labels)} {len(internal_edges)} {len(source_vertices)}\n")
        handle.write("SOURCES" + "".join(f" {int(value)}" for value in source_vertices) + "\n")
        handle.write("WEIGHTS" + "".join(f" {weight:.17g}" for weight in objective_weights) + "\n")
        for face, entries in enumerate(face_entries):
            values = "".join(f" {label_id} {mass:.17g}" for label_id, mass in entries)
            handle.write(f"FACE {face} {len(entries)}{values}\n")
        for output_edge, mesh_edge in enumerate(internal_edges):
            vertex0, vertex1 = (int(value) for value in seams.source_edges[mesh_edge])
            face0, face1 = (int(value) for value in seams.incident_faces[mesh_edge])
            representative = seams.edges[mesh_edge]
            edge_labels = sorted(
                {str(vertex_labels[int(representative[0])]), str(vertex_labels[int(representative[1])])}
            )
            label_values = "".join(f" {label_ids[label]}" for label in edge_labels)
            handle.write(
                f"EDGE {output_edge} {vertex0} {vertex1} {face0} {face1} "
                f"{int(seams.seam_mask[mesh_edge])} {len(edge_labels)}{label_values}\n"
            )

    return {
        "schema_version": 2,
        "face_count": int(len(mesh.faces)),
        "mesh_vertex_count": int(len(mesh.vertices)),
        "input_vertex_count": int(len(source_vertices)),
        "residue_count": int(len(labels)),
        "internal_edge_count": int(len(internal_edges)),
        "initial_seam_edge_count": int(np.count_nonzero(seams.seam_mask)),
        "residue_labels": labels,
        "residue_weights": objective_weights,
    }


def residue_fragmentation_report(
    meshes: Sequence[trimesh.Trimesh],
    source_labels: Sequence[str],
    *,
    uv_key: str = "uv",
    interaction_weights: Mapping[str, float] | None = None,
    objective_weights: Mapping[str, float] | None = None,
    atol: float = 1e-9,
) -> dict[str, object]:
    """Aggregate residue-footprint fragmentation across an atlas.

    Natural footprint disconnections are the baseline.  Within each original
    component, the score is ``1 - sum(piece_mass_fraction ** 2)`` and is then
    averaged by that component's share of the residue footprint mass.
    """

    aggregate: dict[str, dict[str, float | int | str]] = {}
    total_seam_edges = 0
    for mesh in meshes:
        vertex_labels = mesh_vertex_residue_labels(mesh, source_labels)
        records, seams = analyze_residue_footprints(
            mesh,
            as_corner_uv(mesh, key=uv_key),
            vertex_labels,
            atol=atol,
        )
        total_seam_edges += int(np.count_nonzero(seams.seam_mask))
        for record in records:
            label = str(record["residue"])
            target = aggregate.setdefault(
                label,
                {
                    "residue": label,
                    "footprint_area": 0.0,
                    "face_count": 0,
                    "dual_edge_count": 0,
                    "baseline_component_count": 0,
                    "component_count_after_seams": 0,
                    "extra_component_count": 0,
                    "cycle_rank": 0,
                    "seam_crossing_edge_count": 0,
                    "nonseparating_seam_crossing_edge_count": 0,
                    "seam_crossing_length_3d": 0.0,
                    "fragmentation_mass": 0.0,
                },
            )
            for key in (
                "footprint_area",
                "face_count",
                "dual_edge_count",
                "baseline_component_count",
                "component_count_after_seams",
                "extra_component_count",
                "cycle_rank",
                "seam_crossing_edge_count",
                "nonseparating_seam_crossing_edge_count",
                "seam_crossing_length_3d",
                "fragmentation_mass",
            ):
                target[key] += record[key]  # type: ignore[operator]

    weights = interaction_weights or {}
    objective = objective_weights or {}
    residue_records = []
    for label in sorted(aggregate):
        record = aggregate[label]
        area = float(record.pop("footprint_area"))
        fragmentation_mass = float(record.pop("fragmentation_mass"))
        residue_records.append(
            {
                **record,
                "footprint_area": area,
                "fragmentation": float(fragmentation_mass / area) if area > _EPS else 0.0,
                "interaction_weight": float(weights.get(label, 0.0)),
                "objective_weight": float(objective.get(label, 1.0 + float(weights.get(label, 0.0)))),
            }
        )

    residue_count = len(residue_records)
    footprint_area = float(sum(float(record["footprint_area"]) for record in residue_records))
    fragmentation_mass = float(
        sum(float(record["footprint_area"]) * float(record["fragmentation"]) for record in residue_records)
    )
    interaction_weight_sum = float(sum(float(record["interaction_weight"]) for record in residue_records))
    objective_weight_sum = float(sum(float(record["objective_weight"]) for record in residue_records))
    interaction_fragmentation = (
        float(
            sum(float(record["interaction_weight"]) * float(record["fragmentation"]) for record in residue_records)
            / interaction_weight_sum
        )
        if interaction_weight_sum > 0.0
        else None
    )
    objective_fragmentation = (
        float(
            sum(float(record["objective_weight"]) * float(record["fragmentation"]) for record in residue_records)
            / objective_weight_sum
        )
        if objective_weight_sum > 0.0
        else 0.0
    )
    cycle_rank = int(sum(int(record["cycle_rank"]) for record in residue_records))
    cyclic_residues = sum(int(record["cycle_rank"]) > 0 for record in residue_records)
    nonseparating = int(sum(int(record["nonseparating_seam_crossing_edge_count"]) for record in residue_records))
    observed_nonlocal_residues = sum(
        int(record["nonseparating_seam_crossing_edge_count"]) > 0 for record in residue_records
    )
    return {
        "definition": "mass_aware_fragmentation_within_original_face_dual_components",
        "footprint_mass": "integrated_piecewise_linear_residue_corner_indicator",
        "residue_count": int(residue_count),
        "footprint_area": footprint_area,
        "baseline_component_count": int(sum(int(record["baseline_component_count"]) for record in residue_records)),
        "component_count_after_seams": int(
            sum(int(record["component_count_after_seams"]) for record in residue_records)
        ),
        "extra_component_count": int(sum(int(record["extra_component_count"]) for record in residue_records)),
        "fragmented_residue_count": int(sum(float(record["fragmentation"]) > _EPS for record in residue_records)),
        "mean_fragmentation": float(np.mean([float(record["fragmentation"]) for record in residue_records]))
        if residue_records
        else 0.0,
        "area_weighted_fragmentation": float(fragmentation_mass / footprint_area) if footprint_area > _EPS else 0.0,
        "interaction_weighted_fragmentation": interaction_fragmentation,
        "interaction_weight_sum": interaction_weight_sum,
        "objective_weighted_fragmentation": objective_fragmentation,
        "objective_weight_sum": objective_weight_sum,
        "seam_edge_count": int(total_seam_edges),
        "residue_seam_crossing_edge_count": int(
            sum(int(record["seam_crossing_edge_count"]) for record in residue_records)
        ),
        "nonlocality_audit": {
            "status": "nonlocal_structure_present" if cycle_rank > 0 else "tree_like_footprints_only",
            "cycle_rank": cycle_rank,
            "cyclic_residue_count": int(cyclic_residues),
            "cyclic_residue_ratio": float(cyclic_residues / residue_count) if residue_count else 0.0,
            "nonseparating_seam_crossing_edge_count": nonseparating,
            "observed_nonseparating_residue_count": int(observed_nonlocal_residues),
        },
        "residues": residue_records,
    }
