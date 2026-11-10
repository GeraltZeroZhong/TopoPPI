"""P-direct: an explicitly limited, charged table of requested quantities.

No geometry is implied by this table. Unstored queries are unavailable; the only
optional extrapolation is the declared 4/12 Å linear interpolation at 8 Å.
"""

import numpy as np

from ..schema import Encoded
from .wire import pack, unpack


def encode_direct(values, *, units=None):
    """values maps declared names/query keys to scalar or missing (None)."""
    keys = sorted(values)
    valid = np.array([values[k] is not None and np.isfinite(values[k]) for k in keys], dtype=np.uint8)
    array = np.array([values[k] if valid[i] else 0 for i, k in enumerate(keys)], dtype=np.float32)
    blob = pack(
        {"codec": "P-direct", "keys": keys, "units": units or {}, "value_precision": "float32"},
        {"values": array, "valid": valid},
    )
    return Encoded(blob, {"codec": "P-direct", "n_queries": len(keys)})


def decode_direct(blob):
    meta, arrays, _ = unpack(blob)
    return {k: float(arrays["values"][i]) if arrays["valid"][i] else None for i, k in enumerate(meta["keys"])}


def query_direct(table, key):
    return {"value": table.get(key), "available": key in table and table[key] is not None}


def interpolate_scale(table, residue_key, scale=8.0):
    if scale != 8:
        return {"value": None, "available": False, "method": "unsupported scale"}
    a = table.get(f"{residue_key}@4")
    b = table.get(f"{residue_key}@12")
    if a is None or b is None:
        return {"value": None, "available": False, "method": "linear_4_12"}
    return {"value": (a + b) / 2, "available": True, "method": "linear_4_12"}
