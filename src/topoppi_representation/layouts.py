"""Uniform-residue OptCuts layout from a canonical native source."""

import time
from copy import deepcopy

import numpy as np

from .source import canonical_source_bytes

LAYOUT = "aware_uniform"
PROTOCOL = "canonical_native_uniform_residue_optcuts20_v1"


def augment_uniform_layout(source, *, optcuts_bin, timeout=120.0, cancel_event=None, log=None):
    """Return a separate source with uniform layouts; failures remain explicit.

    Constants are lambda=20, distortion bound=4.1, and 120 seconds per patch.
    The label-index metadata required by OptCuts points to native vertices, not
    to atoms. No Parameterizer cleanup, remeshing, or coordinate fitting occurs.
    """
    import trimesh
    from threadpoolctl import threadpool_limits

    from topoppi.atlas.uv import as_corner_uv
    from topoppi.config import OptCutsConfig
    from topoppi.mesh.provenance import (
        OPTCUTS_GEOMETRY_VERTEX_IDS,
        SOURCE_ATOM_IDS,
        SOURCE_FACE_IDS,
        SOURCE_VERTEX_IDS,
    )
    from topoppi.optimization.optcuts import OptCutsUVOptimizer

    before = canonical_source_bytes(source)
    result = deepcopy(source)
    started, reports = time.perf_counter(), []
    with threadpool_limits(limits=1):
        for patch in result.patches:
            if cancel_event is not None and cancel_event.is_set():
                raise InterruptedError("Uniform UV layout cancelled")
            if log:
                log(f"Laying out partner {patch.side + 1}, patch {patch.patch_id + 1}")
            tick = time.perf_counter()
            # The common source contract has already normalized coordinates to f32.
            xyz = np.asarray(patch.xyz, dtype=np.float32).astype(np.float64)
            if not np.array_equal(xyz, patch.xyz):
                raise ValueError("Input coordinates do not satisfy the canonical float32 source contract")
            mesh = trimesh.Trimesh(vertices=xyz, faces=patch.faces, process=False)
            identity = np.arange(len(xyz), dtype=np.int64)
            mesh.metadata.update(
                {
                    SOURCE_ATOM_IDS: identity.copy(),
                    SOURCE_VERTEX_IDS: identity.copy(),
                    OPTCUTS_GEOMETRY_VERTEX_IDS: identity.copy(),
                    SOURCE_FACE_IDS: np.arange(len(patch.faces), dtype=np.int64),
                }
            )
            labels = [result.residue_keys[int(i)] for i in patch.residue_index]
            weights = {key: 1.0 for key in labels}
            entry = {
                "side": int(patch.side),
                "patch": int(patch.patch_id),
                "vertices": len(xyz),
                "faces": len(patch.faces),
                "residues": len(weights),
                "residue_weight_min": 1.0,
                "residue_weight_max": 1.0,
            }
            config = OptCutsConfig(
                optcuts_bin=str(optcuts_bin),
                optcuts_mode=100,
                timeout_sec=float(timeout),
                optcuts_distortion_bound=4.1,
                residue_fragmentation_weight=20.0,
            )
            try:
                OptCutsUVOptimizer(config, cancel_event=cancel_event).optimize_patches(
                    [mesh],
                    initialization="automatic",
                    pack=False,
                    build_report=False,
                    source_residue_labels=labels,
                    residue_weights=weights,
                )
                if not np.array_equal(mesh.faces, patch.faces) or not np.array_equal(mesh.vertices, xyz):
                    raise ValueError("OptCuts modified the frozen native geometry or face order")
                uv = np.asarray(as_corner_uv(mesh, key="uv_optcuts"), dtype=np.float64)
                if uv.shape != (len(patch.faces), 3, 2) or not np.isfinite(uv).all():
                    raise ValueError("OptCuts returned an invalid canonical face-corner layout")
                patch.uv[LAYOUT] = uv.copy()
                entry.update(
                    status="ok",
                    execution=mesh.metadata.get("optcuts_execution"),
                    native_face_corner_order_preserved=True,
                )
            except InterruptedError:
                raise
            except (OSError, RuntimeError, ValueError) as exc:
                if cancel_event is not None and cancel_event.is_set():
                    raise InterruptedError("Uniform UV layout cancelled") from exc
                patch.uv.pop(LAYOUT, None)
                entry.update(status="layout_failed", error=f"{type(exc).__name__}: {exc}")
            entry["seconds"] = time.perf_counter() - tick
            reports.append(entry)
    if canonical_source_bytes(result) != before:
        raise RuntimeError("Uniform layout augmentation changed canonical source bytes")
    for old, new in zip(source.patches, result.patches, strict=False):
        for key, values in old.uv.items():
            if key != LAYOUT and not np.array_equal(values, new.uv[key]):
                raise RuntimeError(f"Uniform layout augmentation changed existing {key} coordinates")
    result.metadata["uniform_layout"] = {
        "protocol": PROTOCOL,
        "layout_key": LAYOUT,
        "input_fields": ["canonical_float32_xyz", "faces", "vertex_residue_index", "residue_keys"],
        "source_contract": "canonical_native_source_only",
        "residue_weighting": "uniform_one",
        "fragmentation_lambda": 20.0,
        "distortion_bound": 4.1,
        "timeout_seconds_per_patch": float(timeout),
        "optcuts_mode": 100,
        "initialization": "automatic",
        "atomic_coordinates_read": False,
        "atomic_contacts_read": False,
        "previous_uv_used_for_initialization": False,
        "parameterizer_used": False,
        "label_index_kind": "native_vertex_index_not_atom_index",
        "canonical_source_bytes_unchanged": True,
        "legacy_uv_arrays_unchanged": True,
        "r0_bytes": len(before),
        "patch_reports": reports,
        "complete": all(r["status"] == "ok" for r in reports),
        "seconds": time.perf_counter() - started,
    }
    return result
