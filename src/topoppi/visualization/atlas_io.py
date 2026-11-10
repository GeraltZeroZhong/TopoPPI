"""Save an interface atlas and its annotations for computation-free rendering."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from zipfile import BadZipFile

import numpy as np
import trimesh
from scipy.spatial import cKDTree

from topoppi.atlas.uv import as_corner_uv
from topoppi.config import VisualizationConfig
from topoppi.io.atom_records import atom_records as _atom_records
from topoppi.io.atom_records import restore_atoms as _restore_atoms
from topoppi.visualization.visualizer import InterfaceVisualizer


@dataclass
class AtlasDocument:
    patches: list[trimesh.Trimesh]
    visualizer: InterfaceVisualizer
    style: dict
    metadata: dict


def _json_value(value):
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_json_value(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def save_atlas(path, patches, visualizer, *, style_config=None, run_metadata=None) -> Path:
    """Write one NPZ with exact meshes, atom identities, resolved annotations and style.

    The file stores numeric/string arrays and JSON. Supply all retained patches,
    including those hidden by the marker view's display threshold, so the atlas
    can be reopened in either display mode.
    """
    if not patches:
        raise ValueError("Cannot save an empty atlas.")
    capture = getattr(visualizer, "capture_surface_camera", None)
    if capture is not None:
        visualizer.last_style["surface_camera"] = capture()
    arrays, patch_records = {}, []
    for index, patch in enumerate(patches):
        prefix = f"patch_{index}"
        arrays[f"{prefix}_vertices"] = np.asarray(patch.vertices, dtype=float)
        arrays[f"{prefix}_faces"] = np.asarray(patch.faces, dtype=np.int64)
        metadata, metadata_arrays = {}, {}
        for key, value in patch.metadata.items():
            if isinstance(value, np.ndarray) and value.dtype.kind != "O":
                array_key = f"{prefix}_meta_{key}"
                arrays[array_key] = value
                metadata_arrays[key] = array_key
            else:
                metadata[key] = _json_value(value)
        # A saved atlas always has a complete display layout, including seams.
        uv_key = "uv_global" if "uv_global" in patch.metadata else "uv"
        uv = as_corner_uv(patch, key=uv_key)
        if not np.isfinite(uv).all():
            raise ValueError("Cannot save an atlas with non-finite UV coordinates.")
        arrays[f"{prefix}_display_uv"] = uv
        patch_records.append({"prefix": prefix, "metadata": metadata, "metadata_arrays": metadata_arrays})

    arrays["coords_A"] = np.asarray(visualizer.coords_A, dtype=float)
    arrays["coords_B"] = np.asarray(visualizer.coords_B, dtype=float)
    last_style = dict(getattr(visualizer, "last_style", {}) or {})
    style = {**asdict(visualizer.config), **last_style}
    if style_config:
        style.update(style_config)
    if style.get("annotation_file") and (style.get("annotation_values") is None or
                                         style["annotation_file"] != last_style.get("annotation_file")):
        from topoppi.visualization.footprint_rendering import read_residue_annotations

        style["annotation_values"] = read_residue_annotations(style["annotation_file"], visualizer.residue_metadata_A)
    state = {
        "config": asdict(visualizer.config),
        "chain_a_id": visualizer.chain_a_id, "chain_b_id": visualizer.chain_b_id,
        "structure_label": visualizer.structure_label,
        "contact_distance_angstrom": visualizer.contact_distance_angstrom,
        "atoms_A": _atom_records(visualizer.atoms_A), "atoms_B": _atom_records(visualizer.atoms_B),
        "interaction_partner_map": visualizer.interaction_partner_map,
        "interaction_type_source": visualizer.interaction_type_source,
        "interaction_residue_source": visualizer.interaction_residue_source,
        "prolif_data": visualizer.prolif_data, "prolif_partners": visualizer.prolif_partners,
        "geometric_types_cache": visualizer._geometric_types_cache,
    }
    document = _json_value({"format": "topoppi-atlas", "version": 1, "patches": patch_records,
                            "visualizer": state, "style": style, "metadata": run_metadata or {}})
    arrays["document_json"] = np.asarray(json.dumps(document, allow_nan=False))
    target = Path(path).expanduser()
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("wb") as stream:
        np.savez_compressed(stream, **arrays)
    return target


def load_atlas(path) -> AtlasDocument:
    """Load a self-contained atlas without structure, annotation or solver files."""
    invalid_file_message = (
        "This file is not a valid TopoPPI atlas archive. "
        "Choose a .npz file created with Save Atlas or --export-atlas. "
        "If the file is incomplete or damaged, export it again or obtain a complete copy."
    )
    try:
        arrays = np.load(Path(path).expanduser(), allow_pickle=False)
    except (ValueError, EOFError, BadZipFile) as exc:
        raise ValueError(invalid_file_message) from exc
    if not isinstance(arrays, np.lib.npyio.NpzFile):
        raise ValueError(invalid_file_message)
    try:
        with arrays:
            if "document_json" not in arrays:
                raise ValueError("This file is not a TopoPPI atlas. Use an atlas saved by --export-atlas or Save atlas.")
            document = json.loads(str(arrays["document_json"].item()))
            if document.get("format") == "topoppi-native-interface":
                raise ValueError("This file contains native interface geometry. For rendering, create an atlas with "
                                 "'topoppi STRUCTURE -A CHAIN -B PARTNER --export-atlas atlas.npz -o map.png', "
                                 "then use 'topoppi render atlas.npz -o figure.svg'.")
            if document.get("format") != "topoppi-atlas" or document.get("version") != 1:
                raise ValueError("Unsupported TopoPPI atlas format.")
            patches = []
            for record in document["patches"]:
                prefix = record["prefix"]
                patch = trimesh.Trimesh(vertices=arrays[f"{prefix}_vertices"].copy(),
                                        faces=arrays[f"{prefix}_faces"].copy(), process=False)
                patch.metadata.update(record["metadata"])
                for key, array_key in record["metadata_arrays"].items():
                    patch.metadata[key] = arrays[array_key].copy()
                display_uv = arrays[f"{prefix}_display_uv"].copy()
                patch.metadata.setdefault("uv_global", display_uv)
                patch.metadata.setdefault("uv", display_uv.copy())
                as_corner_uv(patch, key="uv_global")
                patches.append(patch)
            if not patches:
                raise ValueError("The saved atlas has no patches.")
            state = document["visualizer"]
            coords_a, coords_b = arrays["coords_A"].copy(), arrays["coords_B"].copy()
    except (EOFError, BadZipFile) as exc:
        raise ValueError(invalid_file_message) from exc
    atoms_a = _restore_atoms(state["atoms_A"], coords_a)
    atoms_b = _restore_atoms(state["atoms_B"], coords_b)
    config_fields = dict(state["config"])
    config_fields["highlight_residues"] = tuple(config_fields.get("highlight_residues", ()))
    config_fields["annotation_file"] = ""
    config = VisualizationConfig(**config_fields)
    config.validate()
    viz = InterfaceVisualizer(
        atoms_a, coords_a, coords_b, atoms_b,
        chain_a_id=state["chain_a_id"], chain_b_id=state["chain_b_id"], structure_label=state["structure_label"],
        config=replace(config, use_geometric_interaction_fallback=False),
        contact_distance_angstrom=state["contact_distance_angstrom"],
    )
    viz.config = config
    viz.interaction_partner_map = state["interaction_partner_map"]
    viz.prolif_data = (None if state["prolif_data"] is None else
                      {key: set(values) for key, values in state["prolif_data"].items()})
    viz.prolif_partners = state["prolif_partners"]
    viz.interaction_type_source = state["interaction_type_source"]
    viz.interaction_residue_source = state["interaction_residue_source"]
    viz._interaction_engine = viz.interaction_type_source if viz.prolif_data is not None else None
    cache = state["geometric_types_cache"]
    viz._geometric_types_cache = None if cache is None else {key: set(values) for key, values in cache.items()}
    viz.tree_B = cKDTree(coords_b) if viz.interaction_type_source == "geometric_fallback" and len(coords_b) else None
    style = dict(document["style"])
    style["annotation_file"] = ""
    viz.last_style = dict(style)
    return AtlasDocument(patches=patches, visualizer=viz, style=style, metadata=document["metadata"])
