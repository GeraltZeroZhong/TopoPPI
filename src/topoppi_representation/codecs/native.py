"""Mature QEM + Draco mesh baseline with conservative residue-mass transport."""

import numpy as np
from scipy.sparse import csr_matrix, vstack

from ..schema import Encoded, Surface
from .wire import pack, surface_arrays, surface_from, surface_meta, unpack


def mesh_edges(faces):
    f = np.asarray(faces, dtype=np.int64)
    e = np.vstack((f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]))
    e = np.sort(e, axis=1)
    return np.unique(e[e[:, 0] != e[:, 1]], axis=0)


def vertex_mass(patch, n_residues):
    triangles = patch.xyz[patch.faces]
    a = np.linalg.norm(np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]), axis=1) / 2
    w = np.bincount(patch.faces.ravel(), weights=np.repeat(a / 3, 3), minlength=len(patch.xyz))
    return csr_matrix((w, (np.arange(len(w)), patch.residue_index)), shape=(len(w), n_residues))


def coarse_patch(patch, n_residues, fraction=1.0):
    """Collapse lineage is used only by the encoder, not the decoded reader."""
    xyz, faces = patch.xyz.copy(), patch.faces.copy()
    mass = vertex_mass(patch, n_residues)
    if not 0 < fraction <= 1:
        raise ValueError("fraction must lie in (0,1]")
    if fraction < 1 and len(faces) > 4:
        import fast_simplification as fs

        _, _, collapses = fs.simplify(xyz, faces, target_reduction=1 - fraction, return_collapses=True)
        xyz, faces, mapping = fs.replay_simplification(xyz, faces, collapses)
        if len(faces) == 0 or np.any(mapping < 0) or np.any(mapping >= len(xyz)):
            raise ValueError("QEM lineage lost source support; this candidate is infeasible")
        transport = csr_matrix(
            (np.ones(len(mapping)), (mapping, np.arange(len(mapping)))), shape=(len(xyz), len(mapping))
        )
        mass = transport @ mass
    return xyz, faces, mass


def encode_native(source, fraction=1.0, quantization_bits=18, *, geometry_cache=None):
    import DracoPy

    xyzs, edges, masses, sides, patch_ids, parts, payloads = [], [], [], [], [], [], {}
    offset = 0
    for k, p in enumerate(source.patches):
        key = ("N", float(fraction), k)
        if geometry_cache is not None:
            if key not in geometry_cache:
                geometry_cache[key] = coarse_patch(p, len(source.residue_keys), fraction)
            xyz, faces, mass = geometry_cache[key]
        else:
            xyz, faces, mass = coarse_patch(p, len(source.residue_keys), fraction)
        n = len(xyz)
        # Draco reorders vertices. Explicit integer mass-row identity is compressed
        # inside Draco and charged, permitting its efficient mesh traversal.
        payloads[f"mesh{k}"] = DracoPy.encode(
            np.asarray(xyz, float),
            np.asarray(faces, np.uint32),
            quantization_bits=int(quantization_bits),
            compression_level=7,
            generic_attributes={"mass_row": np.arange(n, dtype=np.uint32)[:, None]},
        )
        xyzs.append(xyz)
        edges.append(mesh_edges(faces) + offset)
        masses.append(mass)
        sides.extend([p.side] * n)
        patch_ids.extend([k] * n)
        parts.append({"side": int(p.side), "patch_id": k, "n": n, "payload": f"mesh{k}"})
        offset += n
    s = Surface(
        np.vstack(xyzs),
        np.vstack(edges),
        vstack(masses, format="csr"),
        np.asarray(sides),
        np.asarray(patch_ids),
        source.residue_keys,
        source.amino_acids,
        source.residue_fields,
        {"representation": "native", "mass_source": "conservative source dual cells"},
    )
    meta = surface_meta(s)
    meta.update(
        codec="N",
        fraction=float(fraction),
        quantization_bits=int(quantization_bits),
        parts=parts,
        mass_row_mapping="integer Draco generic attribute; charged",
    )
    arrays = surface_arrays(s)
    for name in ["xyz", "edges", "side", "patch"]:
        del arrays[name]  # reconstructed from the charged mesh/part description
    blob = pack(meta, arrays, payloads)
    return Encoded(
        blob,
        {
            "n_cells": len(s.xyz),
            "n_patches": len(parts),
            "codec": "N",
            "fraction": fraction,
            "quantization_bits": quantization_bits,
        },
    )


def decode_native(blob):
    import DracoPy

    meta, arrays, payloads = unpack(blob)
    xyz, edges, sides, patch_ids, rows, faces = [], [], [], [], [], []
    offset = 0
    areas = []
    for part in meta["parts"]:
        d = DracoPy.decode(payloads[part["payload"]])
        order = np.asarray(d.get_attribute_by_name("mass_row")["data"]).ravel().astype(int)
        identities, first = np.unique(order, return_index=True)
        if not np.array_equal(identities, np.arange(part["n"])):
            raise ValueError("Draco lost encoded mass-row identities")
        # Draco may split a non-manifold coarse vertex into topological copies.
        # The charged integer attribute identifies that single physical cell;
        # restore it without duplicating its area or inventing proximity welds.
        points = np.asarray(d.points)[first]
        if not np.array_equal(np.asarray(d.points), points[order]):
            raise ValueError("Draco copies of one mass identity have inconsistent XYZ")
        local_faces = order[np.asarray(d.faces, dtype=np.int64)]
        xyz.append(points)
        edges.append(mesh_edges(local_faces) + offset)
        rows.extend(identities + offset)
        faces.append(local_faces + offset)
        sides.extend([part["side"]] * len(identities))
        patch_ids.extend([part["patch_id"]] * len(identities))
        t = points[local_faces]
        areas.append(float(np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1).sum() / 2))
        offset += len(identities)
    arrays["side"] = np.asarray(sides)
    arrays["patch"] = np.asarray(patch_ids)
    s = surface_from(meta, arrays, xyz=np.vstack(xyz), edges=np.vstack(edges))
    s.mass = s.mass[np.asarray(rows)]
    s.faces = np.vstack(faces)
    s.metadata.update(
        decoded_geometry_area=float(sum(areas)),
        decoded_geometry_area_by_patch=areas,
        codec="N",
        fraction=meta["fraction"],
        quantization_bits=meta["quantization_bits"],
    )
    return s
