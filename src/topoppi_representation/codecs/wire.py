"""Small self-contained binary envelope; every header, array and payload is charged."""

import json
import struct

import numpy as np
from scipy.sparse import csr_matrix

from ..schema import Surface

MAGIC = b"RPPI1\0"


def pack(meta, arrays, payloads=None):
    import zstandard as zstd

    header = {"meta": meta, "arrays": {}, "payloads": {}}
    parts, offset = [], 0
    for key, value in sorted(arrays.items()):
        a = np.ascontiguousarray(value)
        if a.dtype.hasobject:
            raise TypeError("Object arrays are not a scientific bitstream")
        a = a.astype(a.dtype.newbyteorder("<"), copy=False)
        raw = a.tobytes()
        header["arrays"][key] = [a.dtype.str, list(a.shape), offset, len(raw)]
        parts.append(raw)
        offset += len(raw)
    for key, raw in sorted((payloads or {}).items()):
        header["payloads"][key] = [offset, len(raw)]
        parts.append(raw)
        offset += len(raw)
    description = json.dumps(header, separators=(",", ":"), sort_keys=True, allow_nan=False).encode()
    raw = struct.pack("<I", len(description)) + description + b"".join(parts)
    return MAGIC + zstd.ZstdCompressor(level=3).compress(raw)


def unpack(blob):
    import zstandard as zstd

    if not blob.startswith(MAGIC):
        raise ValueError("Not an RPPI version 1 code")
    raw = zstd.ZstdDecompressor().decompress(blob[len(MAGIC) :])
    n = struct.unpack("<I", raw[:4])[0]
    header = json.loads(raw[4 : 4 + n])
    data = memoryview(raw)[4 + n :]
    arrays = {
        k: np.frombuffer(data[o : o + length], dtype=d).reshape(s).copy()
        for k, (d, s, o, length) in header["arrays"].items()
    }
    payloads = {k: bytes(data[o : o + length]) for k, (o, length) in header["payloads"].items()}
    return header["meta"], arrays, payloads


def surface_arrays(surface):
    m = surface.mass.tocsr(copy=True)
    m.sum_duplicates()
    m.eliminate_zeros()
    m.sort_indices()
    out = {
        "xyz": np.asarray(surface.xyz, dtype="<f4"),
        "edges": np.asarray(surface.edges, dtype="<u4").reshape(-1, 2),
        "mass_data": m.data.astype("<f4"),
        "mass_indices": m.indices.astype("<u4"),
        "mass_indptr": m.indptr.astype("<u4"),
        "side": np.asarray(surface.side, dtype="u1"),
        "patch": np.asarray(surface.patch, dtype="<i4"),
    }
    for name, value in surface.residue_fields.items():
        out["field:" + name] = np.asarray(value, dtype="<f4")
    if surface.faces is not None:
        out["faces"] = np.asarray(surface.faces, dtype="<u4")
    return out


def surface_meta(surface):
    return {
        "residue_keys": list(surface.residue_keys),
        "amino_acids": list(surface.amino_acids),
        "n_cells": len(surface.xyz),
        "mass_unit": "angstrom^2",
        "xyz_unit": "angstrom",
        "mass_precision": "float32",
        "field_names": sorted(surface.residue_fields),
        "surface_metadata": surface.metadata,
    }


def surface_from(meta, arrays, xyz=None, edges=None):
    n = meta["n_cells"]
    r = len(meta["residue_keys"])
    m = csr_matrix((arrays["mass_data"].astype(float), arrays["mass_indices"], arrays["mass_indptr"]), shape=(n, r))
    return Surface(
        np.asarray(arrays["xyz"] if xyz is None else xyz, dtype=float),
        np.asarray(arrays["edges"] if edges is None else edges, dtype=np.int64).reshape(-1, 2),
        m,
        arrays["side"].astype(int),
        arrays["patch"].astype(int),
        meta["residue_keys"],
        meta["amino_acids"],
        {k: arrays["field:" + k].astype(float) for k in meta.get("field_names", [])},
        meta.get("surface_metadata", {}),
        arrays.get("faces"),
    )


def save_surface(surface, path):
    from pathlib import Path

    Path(path).write_bytes(pack(surface_meta(surface), surface_arrays(surface)))


def load_surface(path):
    from pathlib import Path

    meta, arrays, _ = unpack(Path(path).read_bytes())
    return surface_from(meta, arrays)
