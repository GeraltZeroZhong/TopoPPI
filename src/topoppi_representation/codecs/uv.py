"""Conservative UV carriers and a self-contained Haar/quantization/Zstd code.

Residue dual cells carry source area; XYZ is their area centroid. Disconnected
pieces in a pixel get separate slots. Only positive-length native adjacency is
used to join pieces, including across cuts. Chemical occupancy is derived from
the charged residue-mass table, never separately encoded as a free field.
"""

from collections import defaultdict
from itertools import combinations

import numpy as np
from scipy.sparse import csr_matrix
from shapely.geometry import Polygon, box

from ..schema import Encoded, Surface
from .wire import pack, surface_arrays, surface_from, surface_meta, unpack


class _Union:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, i):
        while self.parent[i] != i:
            self.parent[i] = self.parent[self.parent[i]]
            i = self.parent[i]
        return i

    def join(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.parent[b] = a


def _face_edges(patch, uv):
    edges = defaultdict(list)
    for fi, face in enumerate(patch.faces):
        for k in range(3):
            next_corner = (k + 1) % 3
            a, b = int(face[k]), int(face[next_corner])
            if a > b:
                a, b, k0, k1 = b, a, next_corner, k
            else:
                k0, k1 = k, next_corner
            if np.linalg.norm(patch.xyz[b] - patch.xyz[a]) > 0:
                edges[a, b].append((fi, uv[fi, [k0, k1]]))
    return edges


def _continuous(a, b, scale):
    return np.allclose(a, b, atol=scale * 1e-9, rtol=0)


def _clip_segment(a, b, lower, upper):
    """Canonical native-edge parameter interval inside a grid rectangle."""
    t0, t1 = 0.0, 1.0
    for d in range(2):
        delta = b[d] - a[d]
        if abs(delta) < 1e-14:
            if a[d] < lower[d] - 1e-12 or a[d] > upper[d] + 1e-12:
                return None
        else:
            lo, hi = sorted(((lower[d] - a[d]) / delta, (upper[d] - a[d]) / delta))
            t0, t1 = max(t0, lo), min(t1, hi)
    return (t0, t1) if t1 - t0 > 1e-10 else None


def rasterize(source, layout="aware", target_cells=None):
    """Return (Surface, grid_info), without using experimental labels.

    target_cells is an approximate total occupied-cell budget; thin boundaries
    and disconnected slots can increase the actual count, which is reported.
    """
    target_cells = sum(len(p.xyz) for p in source.patches) if target_cells is None else int(target_cells)
    if target_cells < 1:
        raise ValueError("target_cells must be positive")
    patch_info, total_area = [], 0.0
    for patch in source.patches:
        uv = np.asarray(patch.uv[layout], dtype=float)
        tri = np.asarray(patch.xyz)[patch.faces]
        area = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2
        edge_map = _face_edges(patch, uv)
        scale = max(float(np.ptp(uv.reshape(-1, 2), axis=0).max()), 1e-12)
        groups = _Union(len(patch.faces))
        for entries in edge_map.values():
            for (fa, a), (fb, b) in combinations(entries, 2):
                if _continuous(a, b, scale):
                    groups.join(fa, fb)
        chart_faces = defaultdict(list)
        for fi in range(len(patch.faces)):
            if area[fi] > 0:
                chart_faces[groups.find(fi)].append(fi)
        patch_info.append((patch, uv, area, edge_map, scale, chart_faces))
        total_area += float(area.sum())
    if total_area == 0:
        raise ValueError("No positive-area source faces")

    pieces, charts, joins, seams = [], [], [], []
    for patch_no, (patch, uv, area, edge_map, scale, chart_faces) in enumerate(patch_info):
        # For every native edge, keep intervals in its original 0..1 parameter.
        edge_records = defaultdict(list)
        for face_list in chart_faces.values():
            chart_id = len(charts)
            raw = uv[face_list]
            origin = raw.reshape(-1, 2).min(axis=0)
            normalizer = float(np.ptp(raw.reshape(-1, 2), axis=0).max())
            if normalizer == 0:
                raise ValueError("A positive-area source chart has degenerate UV")
            local = (raw - origin) / normalizer
            extent = local.reshape(-1, 2).max(axis=0)
            e0, e1 = local[:, 1] - local[:, 0], local[:, 2] - local[:, 0]
            uv_area = np.abs(e0[:, 0] * e1[:, 1] - e0[:, 1] * e1[:, 0]) / 2
            if np.any(uv_area <= 1e-16):
                raise ValueError("A positive-area source face has degenerate UV")
            desired = max(1.0, target_cells * float(area[face_list].sum()) / total_area)
            fill = min(1.0, float(uv_area.sum() / np.prod(extent)))
            box_cells = desired / fill
            nx = max(1, int(round(np.sqrt(box_cells * extent[0] / extent[1]))))
            ny = max(1, int(np.ceil(box_cells / nx)))
            step = extent / [nx, ny]
            charts.append(
                {
                    "shape": [ny, nx],
                    "origin": origin.tolist(),
                    "step": (step * normalizer).tolist(),
                    "side": int(patch.side),
                    "patch": patch_no,
                }
            )
            for fi, tuv, auv in zip(face_list, local, uv_area, strict=False):
                face = patch.faces[fi]
                xyz = patch.xyz[face]
                polygon = Polygon(tuv)
                jacobian = area[fi] / auv
                inverse = np.linalg.inv(np.column_stack((tuv[1] - tuv[0], tuv[2] - tuv[0])))
                centroid = tuv.mean(axis=0)
                dual = [
                    Polygon([tuv[k], (tuv[k] + tuv[(k + 1) % 3]) / 2, centroid, (tuv[k] + tuv[(k - 1) % 3]) / 2])
                    for k in range(3)
                ]
                lo = np.maximum(0, np.floor(tuv.min(axis=0) / step).astype(int))
                hi = np.minimum([nx - 1, ny - 1], np.floor(tuv.max(axis=0) / step).astype(int))
                face_pieces = {}
                for row in range(lo[1], hi[1] + 1):
                    for col in range(lo[0], hi[0] + 1):
                        lower, upper = np.array([col, row]) * step, np.array([col + 1, row + 1]) * step
                        cell = box(*lower, *upper)
                        intersection = polygon.intersection(cell)
                        if intersection.area <= 1e-16:
                            continue
                        pmass = defaultdict(float)
                        for k in range(3):
                            mass = dual[k].intersection(cell).area * jacobian
                            if mass > 0:
                                pmass[int(patch.residue_index[face[k]])] += mass
                        uvcenter = np.asarray(intersection.centroid.coords[0])
                        bc = inverse @ (uvcenter - tuv[0])
                        position = xyz[0] + bc[0] * (xyz[1] - xyz[0]) + bc[1] * (xyz[2] - xyz[0])
                        pid = len(pieces)
                        pieces.append(
                            {
                                "chart": chart_id,
                                "row": row,
                                "col": col,
                                "mass": pmass,
                                "xyz": position,
                                "poly": intersection,
                                "side": patch.side,
                                "patch": patch_no,
                            }
                        )
                        face_pieces[row, col] = pid
                        for k in range(3):
                            next_corner = (k + 1) % 3
                            a, b = int(face[k]), int(face[next_corner])
                            p0, p1 = tuv[k], tuv[next_corner]
                            if a > b:
                                a, b, p0, p1 = b, a, p1, p0
                            interval = _clip_segment(p0, p1, lower, upper)
                            if interval is not None:
                                edge_records[a, b].append((fi, pid, *interval))
                # Four-neighbour pixels, only when actual occupied pieces share a segment.
                for (row, col), pa in face_pieces.items():
                    for neighbor in ((row + 1, col), (row, col + 1)):
                        pb = face_pieces.get(neighbor)
                        if pb is not None and pieces[pa]["poly"].intersection(pieces[pb]["poly"]).length > 1e-10:
                            joins.append((pa, pb, False))
        for key, entries in edge_map.items():
            for (fa, a), (fb, b) in combinations(entries, 2):
                continuous = _continuous(a, b, scale)
                ra = [v for v in edge_records[key] if v[0] == fa]
                rb = [v for v in edge_records[key] if v[0] == fb]
                for _, pa, a0, a1 in ra:
                    for _, pb, b0, b1 in rb:
                        if min(a1, b1) - max(a0, b0) > 1e-10:
                            joins.append((pa, pb, not continuous))
                            if not continuous:
                                seams.append((pa, pb))

    union = _Union(len(pieces))
    for a, b, seam in joins:
        pa, pb = pieces[a], pieces[b]
        if not seam and (pa["chart"], pa["row"], pa["col"]) == (pb["chart"], pb["row"], pb["col"]):
            union.join(a, b)
    components = defaultdict(list)
    for i in range(len(pieces)):
        components[union.find(i)].append(i)
    lookup = np.empty(len(pieces), dtype=int)
    xyz, rows, cols, masses, side, patch_ids, node_grid = [], [], [], [], [], [], []
    slots = defaultdict(int)
    for node, members in enumerate(components.values()):
        first = pieces[members[0]]
        grid_key = (first["chart"], first["row"], first["col"])
        slot = slots[grid_key]
        slots[grid_key] += 1
        node_grid.append((*grid_key, slot))
        mass = defaultdict(float)
        numerator = np.zeros(3)
        for pi in members:
            lookup[pi] = node
            piece = pieces[pi]
            numerator += piece["xyz"] * sum(piece["mass"].values())
            for residue, value in piece["mass"].items():
                mass[residue] += value
        xyz.append(numerator / sum(mass.values()))
        for residue, value in sorted(mass.items()):
            rows.append(node)
            cols.append(residue)
            masses.append(value)
        side.append(first["side"])
        patch_ids.append(first["patch"])
    edges = sorted({tuple(sorted((int(lookup[a]), int(lookup[b])))) for a, b, _ in joins if lookup[a] != lookup[b]})
    seam_edges = sorted({tuple(sorted((int(lookup[a]), int(lookup[b])))) for a, b in seams if lookup[a] != lookup[b]})
    node_grid = np.asarray(node_grid, dtype=np.int32)
    for ci, chart in enumerate(charts):
        chart["slots"] = int(node_grid[node_grid[:, 0] == ci, 3].max() + 1)
    metadata = {
        "representation": "uv",
        "layout": layout,
        "target_cells": target_cells,
        "chart_count": len(charts),
        "seam_edge_count": len(seam_edges),
        "extra_slot_count": int(sum(v - 1 for v in slots.values())),
        "adjacency": "positive-length occupied pieces; explicit native seams",
        "mass_measure": "source residue dual-cell area",
    }
    surface = Surface(
        np.asarray(xyz),
        np.asarray(edges, dtype=int).reshape(-1, 2),
        csr_matrix((masses, (rows, cols)), shape=(len(xyz), len(source.residue_keys))),
        np.asarray(side),
        np.asarray(patch_ids),
        list(source.residue_keys),
        list(source.amino_acids),
        dict(source.residue_fields),
        metadata,
    )
    return surface, {
        "charts": charts,
        "node_grid": node_grid,
        "seam_edges": np.asarray(seam_edges, dtype=int).reshape(-1, 2),
    }


def _encoding_arrays(surface, grid, layout, quant_step, transform):
    """Finite representation before entropy coding."""
    import pywt

    if quant_step < 0:
        raise ValueError("quant_step must be nonnegative")
    # The diagnostic raster description must not acquire coding state in place.
    grid = {**grid, "charts": [dict(chart) for chart in grid["charts"]]}
    arrays = surface_arrays(surface)
    source_xyz = arrays.pop("xyz")
    arrays["uv_node_grid"] = grid["node_grid"]  # Also the sparse occupancy mask.
    for ci, chart in enumerate(grid["charts"]):
        h, w = chart["shape"]
        level = pywt.dwt_max_level(min(h, w), "haar") if transform else 0
        chart["level"] = int(level)
        plane = np.zeros((chart["slots"], 3, h, w), dtype=float)
        nodes = np.flatnonzero(grid["node_grid"][:, 0] == ci)
        origin = source_xyz[nodes[0]].astype(float)
        chart["xyz_origin"] = origin.tolist()
        for node in nodes:
            _, row, col, slot = grid["node_grid"][node]
            plane[slot, :, row, col] = source_xyz[node].astype(float) - origin
        bands = []
        for slot in range(chart["slots"]):
            channels = []
            for d in range(3):
                values = plane[slot, d]
                if transform:
                    values, _ = pywt.coeffs_to_array(pywt.wavedec2(values, "haar", mode="symmetric", level=level))
                channels.append(values)
            bands.append(channels)
        values = np.asarray(bands)
        if quant_step:
            values = np.rint(values / quant_step)
            dtype = np.int32 if np.max(np.abs(values)) <= np.iinfo(np.int32).max else np.int64
            arrays[f"uv_xyz:{ci}"] = values.astype(dtype)
        else:
            # q=0 is an explicitly high-precision transport diagnostic. Keeping
            # residuals double avoids an extra rounding after changing origin.
            arrays[f"uv_xyz:{ci}"] = values.astype(np.float64)
    meta = surface_meta(surface)
    name = {"ordinary": "Uo", "aware": "Ua", "aware_uniform": "Ua0"}[layout]
    meta["codec"] = name + ("" if transform else "-id")
    meta["uv_codec"] = {
        "version": 1,
        "charts": grid["charts"],
        "layout": layout,
        "transform": bool(transform),
        "wavelet": "haar",
        "mode": "symmetric",
        "quant_step": float(quant_step),
        "occupancy": "uv_node_grid",
        "zstd_level": 3,
        "xyz_source_precision": "float32",
        "xyz_frame": "charged chart translation in common XYZ frame",
        "unquantized_storage_precision": "float64",
    }
    return meta, arrays


def encode_uv(source, layout="aware", target_cells=None, quant_step=0.001, transform=True):
    """The envelope includes XYZ coefficients, occupancy, topology and all sidecars."""
    surface, grid = rasterize(source, layout, target_cells)
    meta, arrays = _encoding_arrays(surface, grid, layout, quant_step, transform)
    blob = pack(meta, arrays)
    return Encoded(
        blob,
        {
            "n_cells": len(surface.xyz),
            "charts": len(grid["charts"]),
            "extra_slot_count": surface.metadata["extra_slot_count"],
            "seam_edges": len(grid["seam_edges"]),
        },
    )


def _decode_arrays(meta, arrays):
    """Recover the carrier from a finite representation, before or after Zstd."""
    import pywt

    spec = meta["uv_codec"]
    grid = arrays["uv_node_grid"]
    xyz = np.empty((meta["n_cells"], 3), dtype=float)
    for ci, chart in enumerate(spec["charts"]):
        h, w = chart["shape"]
        values = arrays[f"uv_xyz:{ci}"].astype(float)
        if spec["quant_step"]:
            values *= spec["quant_step"]
        nodes = np.flatnonzero(grid[:, 0] == ci)
        if spec["transform"]:
            template = pywt.wavedec2(np.zeros((h, w)), "haar", mode="symmetric", level=chart["level"])
            _, slices = pywt.coeffs_to_array(template)
        for slot in range(chart["slots"]):
            plane = []
            for d in range(3):
                field = values[slot, d]
                if spec["transform"]:
                    coeffs = pywt.array_to_coeffs(field, slices, output_format="wavedec2")
                    field = pywt.waverec2(coeffs, "haar", mode="symmetric")[:h, :w]
                plane.append(field)
            plane = np.asarray(plane)
            for node in nodes[grid[nodes, 3] == slot]:
                _, row, col, _ = grid[node]
                xyz[node] = plane[:, row, col] + np.asarray(chart["xyz_origin"])
    return surface_from(meta, arrays, xyz=xyz)


def decode_uv(blob):
    """Decode without source geometry, source UVs, source lineage or labels."""
    meta, arrays, _ = unpack(blob)
    return _decode_arrays(meta, arrays)
