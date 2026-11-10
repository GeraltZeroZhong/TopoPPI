"""Self-contained native interface records, with no UV or model requirement."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import trimesh

from topoppi.config import SurfaceConfig, TopologyConfig
from topoppi.interface import PartnerInterface, PreparedInterface
from topoppi.io.atom_records import atom_records, restore_atoms


def _json(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def save_interface(prepared: PreparedInterface, path: str | Path) -> Path:
    """Save native geometry and selected atoms to a self-contained NPZ record."""
    arrays = {}

    def mesh_record(mesh, prefix):
        if mesh is None:
            return None
        arrays[prefix + "_vertices"] = np.asarray(mesh.vertices)
        arrays[prefix + "_faces"] = np.asarray(mesh.faces)
        return {"prefix": prefix, "metadata": mesh.metadata}

    partners = []
    for side, partner in enumerate(prepared.partners):
        arrays[f"coords_{side}"] = partner.coords
        partners.append({"chains": partner.chains, "atoms": atom_records(partner.atoms),
            "mesh": mesh_record(partner.mesh, f"mesh_{side}"),
            "patches": [mesh_record(p, f"patch_{side}_{i}") for i, p in enumerate(partner.patches)],
            "surface_report": partner.surface_report, "topology_report": partner.topology_report})
    record = {"format": "topoppi-native-interface", "version": 1,
              "interface": prepared.metadata(), "partners": partners}
    arrays["document_json"] = np.asarray(json.dumps(record, default=_json))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    return path


def load_interface(path: str | Path) -> PreparedInterface:
    """Load native geometry, atom identities and preparation metadata from NPZ."""
    with np.load(path, allow_pickle=False) as arrays:
        record = json.loads(str(arrays["document_json"]))
        if record["format"] != "topoppi-native-interface" or record["version"] != 1:
            raise ValueError("Unsupported native interface format")

        def restore_mesh(item):
            if item is None:
                return None
            mesh = trimesh.Trimesh(arrays[item["prefix"] + "_vertices"].copy(),
                                   arrays[item["prefix"] + "_faces"].copy(), process=False)
            mesh.metadata.update(item["metadata"])
            # Provenance arrays are consumed by geometry and footprint code.
            for key, value in list(mesh.metadata.items()):
                if isinstance(value, list):
                    mesh.metadata[key] = np.asarray(value)
            return mesh

        partners = []
        for side, item in enumerate(record["partners"]):
            coords = arrays[f"coords_{side}"].copy()
            partners.append(PartnerInterface(tuple(item["chains"]), coords, restore_atoms(item["atoms"], coords),
                restore_mesh(item["mesh"]), [restore_mesh(p) for p in item["patches"]],
                item["surface_report"], item["topology_report"]))
    metadata = record["interface"]
    prepared = PreparedInterface(metadata["analysis_id"], metadata["source"], tuple(partners),
        SurfaceConfig(**metadata["surface"]), TopologyConfig(**metadata["topology"]), metadata["numbering"])
    if "contact_distance_angstrom" in metadata:
        prepared.interface.contact_distance_angstrom = float(metadata["contact_distance_angstrom"])
    return prepared
