"""Area-conserving quadrature and portable representation inputs (Å, Å²)."""

import io
import json
import struct
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

from .schema import ComplexSource, Patch, Surface


def face_areas(patch):
    tri = patch.xyz[patch.faces]
    return np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) * 0.5


def vertex_dual_mass(patch):
    return np.bincount(patch.faces.ravel(), weights=np.repeat(face_areas(patch) / 3, 3), minlength=len(patch.xyz))


def mesh_edges(faces):
    if not len(faces):
        return np.empty((0, 2), dtype=np.int64)
    pairs = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    pairs = np.sort(pairs, axis=1)
    return np.unique(pairs[pairs[:, 0] != pairs[:, 1]], axis=0)


def _surface(source, xyzs, rows, columns, masses, sides, patch_ids, edges, metadata):
    xyz = np.concatenate(xyzs) if xyzs else np.empty((0, 3), dtype=np.float64)
    row = np.concatenate(rows) if rows else np.empty(0, dtype=np.int64)
    col = np.concatenate(columns) if columns else np.empty(0, dtype=np.int64)
    values = np.concatenate(masses) if masses else np.empty(0, dtype=np.float64)
    mass = csr_matrix((values, (row, col)), shape=(len(xyz), len(source.residue_keys)))
    return Surface(
        xyz,
        np.concatenate(edges) if edges else np.empty((0, 2), dtype=np.int64),
        mass,
        np.concatenate(sides) if sides else np.empty(0, dtype=np.int8),
        np.concatenate(patch_ids) if patch_ids else np.empty(0, dtype=np.int32),
        list(source.residue_keys),
        list(source.amino_acids),
        {k: np.asarray(v).copy() for k, v in source.residue_fields.items()},
        metadata,
    )


def native_surface(source):
    """One barycentric dual cell per identity-welded vertex; edges stay native."""
    xyzs, rows, columns, masses, sides, patch_ids, edges = [], [], [], [], [], [], []
    offset = 0
    for physical_patch, patch in enumerate(source.patches):
        n = len(patch.xyz)
        xyzs.append(patch.xyz.copy())
        rows.append(np.arange(n) + offset)
        columns.append(patch.residue_index)
        masses.append(vertex_dual_mass(patch))
        sides.append(np.full(n, patch.side, dtype=np.int8))
        patch_ids.append(np.full(n, physical_patch, dtype=np.int32))
        edges.append(mesh_edges(patch.faces) + offset)
        offset += n
    result = _surface(
        source,
        xyzs,
        rows,
        columns,
        masses,
        sides,
        patch_ids,
        edges,
        {"carrier": "native_vertex_dual_cells", "mass_units": "angstrom_squared"},
    )
    offsets = np.cumsum([0] + [len(p.xyz) for p in source.patches[:-1]])
    result.faces = (
        np.concatenate([p.faces + off for p, off in zip(source.patches, offsets, strict=False)])
        if source.patches
        else np.empty((0, 3), dtype=np.int64)
    )
    return result


def _subtriangle_centroids(subdivision):
    """Equal-area n² subdivisions in triangle coordinates (1-u-v,u,v)."""
    n = int(subdivision)
    if n < 1 or n != subdivision:
        raise ValueError("subdivision must be a positive integer")
    uv = []
    for i in range(n):
        for j in range(n - i):
            uv.append(((i + 1 / 3) / n, (j + 1 / 3) / n))
            if i + j < n - 1:
                uv.append(((i + 2 / 3) / n, (j + 2 / 3) / n))
    uv = np.asarray(uv)
    return np.column_stack([1 - uv.sum(axis=1), uv])


def _dual_face_barycentric(subdivision):
    bary = _subtriangle_centroids(subdivision)
    eye, center = np.eye(3), np.full(3, 1 / 3)
    dual = []
    for corner in range(3):
        v = eye[corner]
        dual.extend([[v, (v + eye[(corner + 1) % 3]) / 2, center], [v, center, (v + eye[(corner + 2) % 3]) / 2]])
    return np.einsum("qk,tkc->tqc", bary, np.asarray(dual)).reshape(-1, 3)


def surface_quadrature(surface, faces=None, subdivision=1):
    """Spread carried vertex mass over decoded dual geometry, conservatively.

    Per-corner geometric subarea / full vertex dual area multiplies that
    vertex's entire sparse residue-mass row. Thus geometry changes position,
    but does not silently replace the carried source measure. Vertices without
    positive-area support retain their carried mass as explicit point masses.
    """
    faces = surface.faces if faces is None else faces
    if faces is None:
        raise ValueError("surface_quadrature requires native triangle faces")
    faces = np.asarray(faces, dtype=np.int64)
    bary = _dual_face_barycentric(subdivision)
    count = len(bary)
    tri = surface.xyz[faces]
    areas = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2
    dual_mass = np.bincount(faces.ravel(), weights=np.repeat(areas / 3, 3), minlength=len(surface.xyz))
    vertex_rows = np.repeat(faces, count // 3, axis=1).ravel()
    geometric_mass = np.repeat(areas / count, count)
    factors = np.divide(
        geometric_mass, dual_mass[vertex_rows], out=np.zeros_like(geometric_mass), where=dual_mass[vertex_rows] > 0
    )
    xyz = np.einsum("qc,fcd->fqd", bary, tri).reshape(-1, 3)
    positive = factors > 0
    xyz, vertex_rows, factors = xyz[positive], vertex_rows[positive], factors[positive]
    isolated = np.flatnonzero((dual_mass == 0) & (surface.weights > 0))
    if len(isolated):
        xyz = np.concatenate([xyz, surface.xyz[isolated]])
        vertex_rows = np.concatenate([vertex_rows, isolated])
        factors = np.concatenate([factors, np.ones(len(isolated))])
    mass = surface.mass[vertex_rows].multiply(factors[:, None]).tocsr()
    return Surface(
        xyz,
        np.empty((0, 2), dtype=np.int64),
        mass,
        surface.side[vertex_rows],
        surface.patch[vertex_rows],
        list(surface.residue_keys),
        list(surface.amino_acids),
        {k: np.asarray(v).copy() for k, v in surface.residue_fields.items()},
        {
            "carrier": "barycentric_dual_quadrature",
            "subdivision": int(subdivision),
            "points_per_face": count,
            "integration_only": True,
            "mass_units": "angstrom_squared",
            "isolated_mass_point_count": len(isolated),
        },
    )


def reference_quadrature(source, subdivision=1):
    """Integrate piecewise-constant corner dual cells with 6*n² points/face.

    Each quadrilateral corner region is split at the face centroid into two
    equal-area triangles. Their deterministic centroid rules preserve each
    residue's mass exactly up to floating-point summation. Empty edges are
    intentional: use native_surface, not these integration points, for hydrophobic connectivity.
    """
    return surface_quadrature(native_surface(source), subdivision=subdivision)


def canonical_source_bytes(source):
    """Return the actual canonical uncompressed payload; R0 is len(payload).

    v1 uses little-endian float32 XYZ/mass/residue fields, uint32 connectivity
    and indices, UTF-8 length-prefixed keys and one AA byte per residue.
    UV, cut IDs, labels, targets, case/family names and provenance are excluded.
    Dual mass is included explicitly as part of the common source measure.
    """
    out = io.BytesIO()
    out.write(b"RP_SRC_1")

    def integers(*values):
        out.write(struct.pack("<" + "I" * len(values), *values))

    def string(value):
        encoded = str(value).encode("utf-8")
        integers(len(encoded))
        out.write(encoded)

    integers(len(source.residue_keys), len(source.patches), len(source.residue_fields))
    for key, aa in zip(source.residue_keys, source.amino_acids, strict=False):
        string(key)
        out.write(aa.encode("ascii"))
    for patch in sorted(source.patches, key=lambda p: (p.side, p.patch_id)):
        integers(patch.side, patch.patch_id, len(patch.xyz), len(patch.faces))
        out.write(np.asarray(patch.xyz, dtype="<f4").tobytes())
        out.write(np.asarray(patch.faces, dtype="<u4").tobytes())
        out.write(np.asarray(patch.residue_index, dtype="<u4").tobytes())
        out.write(np.asarray(vertex_dual_mass(patch), dtype="<f4").tobytes())
    for name in sorted(source.residue_fields):
        values = np.asarray(source.residue_fields[name], dtype="<f4")
        string(name)
        integers(values.ndim, *values.shape)
        out.write(values.tobytes())
    return out.getvalue()


def save_source(source, path):
    """Save source NPZ plus JSON provenance; archive file size is not R0."""
    path = Path(path).with_suffix(".npz")
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays, patch_records = {}, []
    for i, patch in enumerate(source.patches):
        prefix = f"patch_{i}"
        for name in ("xyz", "faces", "residue_index"):
            arrays[f"{prefix}_{name}"] = getattr(patch, name)
        for name, values in patch.uv.items():
            arrays[f"{prefix}_uv_{name}"] = values
        patch_records.append({"side": patch.side, "patch_id": patch.patch_id, "uv": list(patch.uv)})
    fields = sorted(source.residue_fields)
    for i, name in enumerate(fields):
        arrays[f"residue_field_{i}"] = np.asarray(source.residue_fields[name])
    np.savez(path, **arrays)
    description = {
        "sample_id": source.sample_id,
        "family_id": source.family_id,
        "residue_keys": source.residue_keys,
        "amino_acids": source.amino_acids,
        "patches": patch_records,
        "residue_fields": fields,
        "metadata": source.metadata,
    }
    description_path = path.with_suffix(".json")
    temporary = description_path.with_name(description_path.name + ".tmp")
    temporary.write_text(json.dumps(description, indent=2) + "\n")
    temporary.replace(description_path)
    return path


def load_source(path):
    path = Path(path).with_suffix(".npz")
    desc = json.loads(path.with_suffix(".json").read_text())
    with np.load(path, allow_pickle=False) as arrays:
        patches = [
            Patch(
                p["side"],
                p["patch_id"],
                arrays[f"patch_{i}_xyz"],
                arrays[f"patch_{i}_faces"],
                arrays[f"patch_{i}_residue_index"],
                {k: arrays[f"patch_{i}_uv_{k}"] for k in p["uv"]},
            )
            for i, p in enumerate(desc["patches"])
        ]
        fields = {name: arrays[f"residue_field_{i}"] for i, name in enumerate(desc["residue_fields"])}
    return ComplexSource(
        desc["sample_id"],
        desc["family_id"],
        desc["residue_keys"],
        desc["amino_acids"],
        patches,
        fields,
        desc["metadata"],
    )
