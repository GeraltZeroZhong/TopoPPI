"""Load residue graph arrays from explicit local manifests."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def load_graphs(directory):
    """Read graph manifests and their named NPZ array files."""
    directory = Path(directory)
    manifest = json.loads((directory / "graphs.json").read_text())
    result = {}
    for complex_id, filename in manifest.items():
        with np.load(directory / filename, allow_pickle=False) as arrays:
            result[complex_id] = {key: arrays[key] for key in arrays.files if key != "metadata_json"}
    return result
