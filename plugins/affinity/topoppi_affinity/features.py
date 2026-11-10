"""Label-free native TopoPPI geometry features for affinity models.

Both complete partner groups are surfaced independently, then aggregated with
sum/absolute-difference operators. Native 3-D footprints are integrated before
UV cuts; these descriptors do not measure OptCuts or atlas-layout performance.
Mutations change residue identities on the fixed WT structure. They are not
relaxed mutant structures and are not a physical free-energy calculation.
"""
from __future__ import annotations

import json
import re
from dataclasses import asdict
from itertools import combinations
from pathlib import Path

import numpy as np
from Bio.Data.PDBData import protein_letters_3to1_extended
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from topoppi.atlas.footprints import (
    geometric_contact_partner_map,
    mesh_vertex_residue_labels,
    source_atom_residue_labels,
)
from topoppi.config import SurfaceConfig, TopologyConfig
from topoppi.interface import chain_group, load_partner_data, prepare_interface
from topoppi.mesh.provenance import connected_face_components

AA = "ACDEFGHIKLMNPQRSTVWY"
# Standard residue descriptors: Kyte-Doolittle hydropathy, approximate residue
# volume (A^3), formal charge with neutral His, polarity, aromaticity.
_HYDRO = dict(zip(AA, [1.8, 2.5, -3.5, -3.5, 2.8, -.4, -3.2, 4.5, -3.9, 3.8,
                      1.9, -3.5, -1.6, -3.5, -4.5, -.8, -.7, 4.2, -.9, -1.3], strict=True))
_VOLUME = dict(zip(AA, [88.6, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7,
                       168.6, 166.7, 162.9, 114.1, 112.7, 143.8, 173.4, 89.,
                       116.1, 140., 227.8, 193.6], strict=True))
PROPERTIES = ("hydropathy", "volume", "charge", "polar", "aromatic")
CACHE_PROTOCOL = "native_affinity_v1"


def residue_properties(aa: str) -> dict[str, float]:
    return {"hydropathy": float(_HYDRO.get(aa, 0)), "volume": float(_VOLUME.get(aa, 0)),
            "charge": float((aa in "KR") - (aa in "DE")),
            "polar": float(aa in "CDEHKNQRSTY"), "aromatic": float(aa in "FWY")}


def feature_group(name: str) -> str:
    """Groups retain an honest surface ablation for mutation features too."""
    if name.startswith("mutation.surface_"):
        return "surface"
    if name.startswith("mutation.contact_"):
        return "contact"
    return name.split(".", 1)[0]


def _chains(value) -> tuple[str, ...]:
    return chain_group(value)


def _group_data(loader, chains):
    return load_partner_data(loader, chains)


def _footprint_areas(mesh, source_labels):
    labels = mesh_vertex_residue_labels(mesh, source_labels)
    # Exactly the same corner-indicator mass integral as atlas.footprints:
    # each triangle distributes one third of its physical area to each corner.
    corner_labels = labels[np.asarray(mesh.faces)].reshape(-1)
    masses = np.repeat(np.asarray(mesh.area_faces) / 3., 3)
    unique, inv = np.unique(corner_labels, return_inverse=True)
    return dict(zip(unique, np.bincount(inv, weights=masses), strict=True))


def _side_surface(prepared, residues):
    mesh, atoms = prepared.mesh, prepared.atoms
    labels = source_atom_residue_labels(atoms)
    for label, area in _footprint_areas(mesh, labels).items():
        residues[str(label)]["surface_area"] = float(area)
    patches = prepared.patches
    areas, boundaries, bending, compactness = [], [], [], []
    for patch in patches:
        area = float(patch.area)
        areas.append(area)
        edges, counts = np.unique(np.sort(patch.edges, axis=1), axis=0, return_counts=True)
        boundary_edges = edges[counts == 1]
        boundary = float(np.linalg.norm(patch.vertices[boundary_edges[:, 0]] -
                                        patch.vertices[boundary_edges[:, 1]], axis=1).sum())
        boundaries.append(boundary)
        # Integrated unsigned dihedral bending, divided by area below.
        adjacency = np.asarray(patch.face_adjacency)
        if len(adjacency):
            normal_dot = np.sum(patch.face_normals[adjacency[:, 0]] *
                                patch.face_normals[adjacency[:, 1]], axis=1)
            bending.append(float(np.sum(np.arccos(np.clip(normal_dot, -1., 1.)) *
                                          np.linalg.norm(np.diff(patch.vertices[patch.face_adjacency_edges],
                                                                  axis=1)[:, 0], axis=1))))
        compactness.append(4. * np.pi * area / max(boundary * boundary, 1.))
        for label, mass in _footprint_areas(patch, labels).items():
            residues[str(label)]["footprint_area"] += float(mass)
            residues[str(label)]["footprint_patch_count"] += 1.
    total = float(sum(areas))
    result = {"total_area": float(mesh.area), "interface_area": total,
              "interface_fraction": total / max(float(mesh.area), 1.),
              "patch_count": float(len(patches)), "patch_largest_area": max(areas, default=0.),
              "boundary_length": float(sum(boundaries)),
              "bending_per_area": float(sum(bending)) / max(total, 1.),
              "compactness_mean": float(np.mean(compactness)) if compactness else 0.,
              "surface_component_count": float(len(connected_face_components(mesh)))}
    side_records = [residues[str(k)] for k in set(labels)]
    footprints = [r["footprint_area"] for r in side_records if r["footprint_area"] > 0]
    result["footprint_residue_count"] = float(len(footprints))
    result["footprint_area_mean"] = float(np.mean(footprints)) if footprints else 0.
    result["footprint_area_std"] = float(np.std(footprints)) if footprints else 0.
    for prop in PROPERTIES:
        result[f"interface_{prop}_area"] = sum(r["footprint_area"] * residue_properties(r["aa"])[prop]
                                                for r in side_records)
    report = {"surface_generation": prepared.surface_report,
              "topology_status": prepared.topology_report["status"],
              "selected_face_count": prepared.topology_report.get("selected_face_count", 0),
              "accepted_patch_count": len(patches)}
    return result, report


def extract_complex_features(structure_path, chains_a, chains_b, cache_dir=None, *, prepared=None) -> dict:
    """Extract physical features; no labels, atlas annotations or styles accepted.

    Residue keys preserve the numbering in the supplied PDB/mmCIF (including
    insertion codes). For SKEMPI cleaned mutations use its matching cleaned PDB.
    The cache uses source path/mtime/size to avoid repeating surface generation.
    """
    path = Path(structure_path).resolve()
    # Canonical ordering guarantees bitwise swap invariance as well as symmetric
    # aggregation. It does not select a dominant pair or omit a declared chain.
    group_a, group_b = sorted((_chains(chains_a), _chains(chains_b)))
    if set(group_a) & set(group_b):
        raise ValueError("Partner chain groups overlap.")
    if prepared is None:
        stat = path.stat()
        file_metadata = {"path": str(path), "mtime_ns": stat.st_mtime_ns, "size": stat.st_size}
    else:
        file_metadata = {key: prepared.source.get(key) for key in ("path", "mtime_ns", "size")}
    source = {**file_metadata,
              "groups": [list(group_a), list(group_b)], "protocol": CACHE_PROTOCOL,
              "surface": asdict(SurfaceConfig()), "topology": asdict(TopologyConfig())}
    cache_path = None
    if cache_dir and prepared is None:
        # Preserve full path hierarchy: equal PDB names from different structure
        # sources must not reuse one another's residue numbering or geometry.
        cache_path = Path(cache_dir) / "complexes" / path.parent.relative_to(path.anchor)
        cache_path = cache_path / f"{path.stem}__{'+'.join(group_a)}__{'+'.join(group_b)}.json"
        if cache_path.exists():
            bundle = json.loads(cache_path.read_text())
            if bundle["metadata"]["source"] == source:
                return bundle
    prepared = prepared or prepare_interface(path, group_a, group_b)
    if (prepared.source['path'] != str(path) or
            prepared.surface_config != SurfaceConfig() or prepared.topology_config != TopologyConfig()):
        raise ValueError("Prepared interface differs from the frozen native affinity protocol")
    prepared_a, prepared_b = prepared.partner(group_a), prepared.partner(group_b)
    coords_a, atoms_a = prepared_a.coords, prepared_a.atoms
    coords_b, atoms_b = prepared_b.coords, prepared_b.atoms
    residues = {}
    for side, coords, atoms, other in ((0, coords_a, atoms_a, coords_b), (1, coords_b, atoms_b, coords_a)):
        labels = source_atom_residue_labels(atoms)
        distances = cKDTree(other).query(coords)[0]
        for label in dict.fromkeys(labels):
            indices = np.flatnonzero(labels == label)
            residue = atoms[int(indices[0])].get_parent()
            aa = protein_letters_3to1_extended.get(residue.get_resname(), "X")
            residues[str(label)] = {"aa": aa, "chain": str(residue.get_parent().id),
                "residue_number": int(residue.id[1]), "insertion": str(residue.id[2]).strip(),
                "side": side, "coords": coords[indices].tolist(),
                "centroid": coords[indices].mean(axis=0).tolist(),
                "partner_distance": float(distances[indices].min()), "partners": {},
                "contact_degree": 0., "contact_atom_pairs": 0.,
                "surface_area": 0., "footprint_area": 0., "footprint_patch_count": 0.}
    partner_map = geometric_contact_partner_map(coords_a, atoms_a, coords_b, atoms_b, distance_cutoff=6.)
    for left, partners in partner_map.items():
        for right, count in partners.items():
            residues[left]["partners"][right] = count
            residues[right]["partners"][left] = count
    for record in residues.values():
        record["contact_degree"] = float(len(record["partners"]))
        record["contact_atom_pairs"] = float(sum(record["partners"].values()))
    absolute = {}
    for prefix in ("sequence", "contact"):
        sides = []
        for side, current_group in ((0, group_a), (1, group_b)):
            selected = [r for r in residues.values() if r["side"] == side and
                        (prefix == "sequence" or r["contact_degree"] > 0)]
            side_features = {"residue_count": float(len(selected))}
            if prefix == "sequence":
                side_features["chain_count"] = float(len(current_group))
            for aa in AA:
                side_features[f"fraction_{aa}"] = sum(r["aa"] == aa for r in selected) / max(len(selected), 1)
            for prop in PROPERTIES:
                side_features[f"{prop}_mean"] = sum(residue_properties(r["aa"])[prop]
                                                    for r in selected) / max(len(selected), 1)
            sides.append(side_features)
        for key in sides[0]:
            absolute[f"{prefix}.{key}_sum"] = sides[0][key] + sides[1][key]
            absolute[f"{prefix}.{key}_absdiff"] = abs(sides[0][key] - sides[1][key])
    absolute["contact.residue_pairs"] = float(sum(len(p) for p in partner_map.values()))
    absolute["contact.atom_pairs_6A"] = float(sum(sum(p.values()) for p in partner_map.values()))
    for cutoff in (3., 4.5, 8.):
        absolute[f"contact.atom_pairs_{cutoff:g}A"] = float(cKDTree(coords_a).count_neighbors(cKDTree(coords_b), cutoff))
    absolute["contact.min_distance"] = float(cKDTree(coords_b).query(coords_a)[0].min())
    for prop in ("charge", "hydropathy", "polar", "aromatic"):
        absolute[f"contact.{prop}_products"] = sum(residue_properties(residues[left]["aa"])[prop] *
                                                  residue_properties(residues[right]["aa"])[prop]
                                                  for left, p in partner_map.items() for right in p)
    side_a, report_a = _side_surface(prepared_a, residues)
    side_b, report_b = _side_surface(prepared_b, residues)
    for key in side_a:
        absolute[f"surface.{key}_sum"] = side_a[key] + side_b[key]
        absolute[f"surface.{key}_absdiff"] = abs(side_a[key] - side_b[key])
    keys = list(residues)
    centroids = np.array([residues[k]["centroid"] for k in keys])
    neighborhoods = cKDTree(centroids).query_ball_point(centroids, r=8.)
    for i, key in enumerate(keys):
        record = residues[key]
        neighbors = [residues[keys[j]] for j in neighborhoods[i] if j != i]
        record["local_neighbor_count"] = float(len(neighbors))
        for prop in PROPERTIES:
            record[f"neighbor_{prop}"] = sum(residue_properties(r["aa"])[prop] for r in neighbors)
            record[f"partner_{prop}"] = sum(residue_properties(residues[k]["aa"])[prop]
                                              for k in record["partners"])
    bundle = {"absolute_features": absolute, "feature_groups": {k: feature_group(k) for k in absolute},
              "residues": residues, "metadata": {"source": source, "sides": [report_a, report_b],
                "representation": "native_Gaussian_surface_3D_footprints_no_UV_optimization",
                "geometry_protocol": "fixed_WT_structure_identity_substitution",
                "numbering": "supplied_structure_author_chain_residue_insertion",
                "missing_representation": "hydrogens_waters_ions_cofactors_unresolved_residues"}}
    if not np.isfinite(list(absolute.values())).all():
        raise ValueError("Non-finite structure features.")
    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(bundle))
    return bundle


def parse_mutations(mutations) -> list[dict]:
    """Parse CHAIN:WTnumber[insertion]MT, legacy SKEMPI notation, or records."""
    if isinstance(mutations, str):
        mutations = [m.strip() for m in mutations.split(",") if m.strip()]
    result = []
    for mutation in mutations:
        if isinstance(mutation, dict):
            item = dict(mutation)
            item.setdefault("insertion", "")
        else:
            explicit = re.fullmatch(r"([^:\s,]+):([A-Z])(-?\d+)([A-Z]?)([A-Z])", str(mutation))
            match = explicit or re.fullmatch(r"([A-Z])([^\d\s,]+?)(-?\d+)([A-Z]?)([A-Z])", str(mutation))
            if not match:
                raise ValueError(f"Invalid mutation notation: {mutation!r}; use CHAIN:WTnumberMT, e.g. A:G42V or 1:G42V")
            if explicit:
                chain, wt, position, insertion, mt = match.groups()
            else:
                wt, chain, position, insertion, mt = match.groups()
            item = {"wt": wt, "chain": chain, "position": int(position), "insertion": insertion, "mt": mt}
        if item["wt"] not in AA or item["mt"] not in AA:
            raise ValueError(f"Mutation must use canonical amino acids: {mutation}")
        result.append(item)
    return result


def mutation_features(bundle: dict, mutation_list) -> dict[str, float]:
    """Joint fixed-WT mutation encoding, with explicit distal coverage features.

    All listed WT identities are checked against the supplied structure. Missing
    residues or WT mismatches raise errors and must appear in extraction failures.
    Identity substitutions are removed after validation and encode exactly zero.
    Pair products/distance terms let the joint model represent nonadditivity.
    """
    residues = bundle["residues"]
    lookup = {(r["chain"], r["residue_number"], r["insertion"]): (key, r) for key, r in residues.items()}
    mutations, seen = [], set()
    for mutation in parse_mutations(mutation_list):
        position = (mutation["chain"], int(mutation["position"]), mutation["insertion"])
        if position in seen:
            raise ValueError(f"Duplicate mutation site: {position}")
        seen.add(position)
        if position not in lookup:
            raise ValueError(f"Mutation residue absent from supplied structure: {position}")
        key, residue = lookup[position]
        if residue["aa"] != mutation["wt"]:
            raise ValueError(f"WT identity mismatch for {position}: structure={residue['aa']}, mutation={mutation['wt']}")
        if mutation["wt"] != mutation["mt"]:
            delta = {p: residue_properties(mutation["mt"])[p] - residue_properties(mutation["wt"])[p]
                     for p in PROPERTIES}
            mutations.append((mutation, key, residue, delta))
    out = {"mutation.n_mutations": float(len(mutations))}
    for aa in AA:
        out[f"mutation.wt_{aa}"] = sum(m[0]["wt"] == aa for m in mutations)
        out[f"mutation.mt_{aa}"] = sum(m[0]["mt"] == aa for m in mutations)
        out[f"mutation.delta_{aa}"] = out[f"mutation.mt_{aa}"] - out[f"mutation.wt_{aa}"]
    contexts = {"contact_degree": "contact_degree", "contact_atom_pairs": "contact_atom_pairs",
                "contact_partner_distance": "partner_distance", "surface_footprint_area": "footprint_area",
                "surface_total_area": "surface_area", "surface_patch_count": "footprint_patch_count",
                "contact_local_neighbor_count": "local_neighbor_count"}
    for feature, source in contexts.items():
        out[f"mutation.{feature}_sum"] = sum(m[2][source] for m in mutations)
    out["mutation.surface_covered_count"] = float(sum(m[2]["footprint_area"] > 0 for m in mutations))
    out["mutation.contact_interface_count"] = float(sum(m[2]["partner_distance"] <= 6. for m in mutations))
    out["mutation.contact_near_count"] = float(sum(6. < m[2]["partner_distance"] <= 10. for m in mutations))
    out["mutation.contact_distal_count"] = float(sum(m[2]["partner_distance"] > 10. for m in mutations))
    for prop in PROPERTIES:
        out[f"mutation.delta_{prop}"] = sum(m[3][prop] for m in mutations)
        out[f"mutation.absdelta_{prop}"] = sum(abs(m[3][prop]) for m in mutations)
        out[f"mutation.contact_delta_{prop}_degree"] = sum(m[3][prop] * m[2]["contact_degree"] for m in mutations)
        out[f"mutation.surface_delta_{prop}_footprint"] = sum(m[3][prop] * m[2]["footprint_area"] for m in mutations)
        out[f"mutation.contact_delta_{prop}_neighbor"] = sum(m[3][prop] * m[2][f"neighbor_{prop}"] for m in mutations)
        out[f"mutation.contact_delta_{prop}_partner"] = sum(m[3][prop] * m[2][f"partner_{prop}"] for m in mutations)
    pairs = list(combinations(mutations, 2))
    out["epistasis.pair_count"] = float(len(pairs))
    out["epistasis.cross_partner_pairs"] = float(sum(a[2]["side"] != b[2]["side"] for a, b in pairs))
    distances = [float(cdist(a[2]["coords"], b[2]["coords"]).min()) for a, b in pairs]
    out["epistasis.pair_min_distance"] = min(distances, default=0.)
    out["epistasis.pair_mean_distance"] = float(np.mean(distances)) if distances else 0.
    out["epistasis.close_pair_count"] = float(sum(d <= 8. for d in distances))
    out["epistasis.shared_contact_partners"] = float(sum(len(set(a[2]["partners"]) & set(b[2]["partners"]))
                                                         for a, b in pairs))
    for prop in PROPERTIES:
        out[f"epistasis.delta_{prop}_product_proximity"] = sum(a[3][prop] * b[3][prop] * np.exp(-d / 8.)
                                                               for (a, b), d in zip(pairs, distances, strict=True))
    return {k: float(v) for k, v in out.items()}
