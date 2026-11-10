"""Prepare finite residue-aware UV carriers for spatial descriptor queries."""
from __future__ import annotations

from dataclasses import replace
from types import SimpleNamespace

import numpy as np

from .adapters import from_geometry
from .descriptors.spatial import SpatialDescriptorCache
from .source import canonical_source_bytes

SASA_FIELDS = ("sasa.isolated_partner", "sasa.bound_complex", "sasa.buried")
UV_RESOLUTIONS = (1.0, 0.25, 0.0625)
UV_QUANTIZATION_STEPS = (0.001, 0.02, 0.1)


def _check_cancel(cancel_event):
    if cancel_event is not None and cancel_event.is_set():
        raise InterruptedError("Spatial representation preparation cancelled")


def attach_residue_sasa(source, residues):
    """Add the complete partner AA/SASA table, retaining zero-footprint rows.

    ``residues`` is a DataFrame with side, residue_key, aa and the three
    ``sasa.*`` fields. Its side numbering follows sorted partner chain groups.
    """
    rows = residues.to_dict("records")
    table = {f"side{int(row['side'])}:{row['residue_key']}": row for row in rows}
    keys = list(source.residue_keys)
    missing = [key for key in keys if key not in table]
    if missing:
        raise ValueError(f"SASA table lacks surface residue {missing[0]}")
    keys.extend(sorted(set(table) - set(keys)))
    fields = {name: np.pad(np.asarray(values, float), (0, len(keys) - len(values)), constant_values=np.nan)
              for name, values in source.residue_fields.items()}
    fields.update({name: np.array([float(table[key][name]) for key in keys]) for name in SASA_FIELDS})
    metadata = dict(source.metadata)
    metadata["sasa_transport"] = {
        "unit": "angstrom^2", "domain": "complete_partner_table",
        "zero_footprint_rows": len(keys) - len(source.residue_keys),
        "n_points": sorted({int(row["n_points"]) for row in rows if row.get("n_points")}),
    }
    return replace(source, residue_keys=keys, amino_acids=[table[key]["aa"] for key in keys],
                   residue_fields=fields, metadata=metadata)


def prepare_spatial_source(geometry, *, sasa_residues=None, cancel_event=None):
    """Build a canonical source from both Core partner surfaces and WT SASA.

    Partner groups are sorted. Topology preparation preserves native triangle
    coordinates; identity-welded XYZ uses float32. All selected partner residue
    identities and SASA fields contribute to the source byte count.
    """
    from topoppi.config import ParameterizationConfig
    from topoppi.mesh.parameterization import Parameterizer

    partners = []
    for partner in sorted(geometry.partners, key=lambda p: p.chains):
        patches = []
        for patch in partner.layout_patches():
            _check_cancel(cancel_event)
            triangles = np.asarray(patch.triangles).copy()
            carrier = Parameterizer(ParameterizationConfig()).prepare_patch(patch)
            if carrier is None or not np.array_equal(triangles, carrier.triangles):
                raise ValueError("Spatial source preparation requires retained native triangle geometry")
            patches.append(carrier)
        if not patches:
            raise ValueError("Spatial descriptors require a retained surface on each partner")
        partners.append(SimpleNamespace(chains=partner.chains, atoms=partner.atoms, patches=patches))
    if len(partners) != 2:
        raise ValueError("Spatial descriptors require both binding partners")
    source = from_geometry(SimpleNamespace(partners=partners, analysis_id=geometry.analysis_id))
    if sasa_residues is None:
        from .descriptors.physical import residue_physics
        sasa_residues, _ = residue_physics(geometry, n_points=960, probe_radius=1.4)
    return attach_residue_sasa(source, sasa_residues)


def encode_spatial_source(source, *, budget_ratio=0.25, optcuts_bin=None, timeout=120.0,
                          cancel_event=None, log=None):
    """Return a decoded UV surface, its complete stream, and encoding details.

    Uniform residue-aware OptCuts uses lambda 20 and distortion bound 4.1.
    A resolution-first search selects the first stream fitting ``budget_ratio``
    times the canonical source size. The search uses cell fractions 1, 1/4,
    1/16, each with quantization steps 0.001, 0.02 and 0.1 Å. An infeasible
    budget raises ValueError. Existing uniform layouts may be reused.
    """
    from .codecs.uv import _encoding_arrays, decode_uv, rasterize
    from .codecs.wire import pack
    from .layouts import augment_uniform_layout
    from .schema import Encoded

    if not np.isfinite(budget_ratio) or budget_ratio <= 0:
        raise ValueError("budget_ratio must be a positive finite number")
    _check_cancel(cancel_event)
    if not all("aware_uniform" in patch.uv for patch in source.patches):
        if log:
            log("Computing residue-aware UV layouts")
        source = augment_uniform_layout(source, optcuts_bin=optcuts_bin or "OptCuts_bin",
                                        timeout=timeout, cancel_event=cancel_event, log=log)
        failures = [row for row in source.metadata["uniform_layout"]["patch_reports"] if row["status"] != "ok"]
        if failures:
            raise RuntimeError("Residue-aware UV layout failed: " + failures[0].get("error", "missing layout"))
    source_bytes = len(canonical_source_bytes(source))
    budget = int(source_bytes * budget_ratio)
    n_vertices = sum(len(patch.xyz) for patch in source.patches)
    attempts = []
    for fraction in UV_RESOLUTIONS:
        _check_cancel(cancel_event)
        target_cells = max(1, round(n_vertices * fraction))
        if log:
            log(f"Encoding UV surface with {target_cells} target cells")
        surface, grid = rasterize(source, "aware_uniform", target_cells)
        for quant_step in UV_QUANTIZATION_STEPS:
            _check_cancel(cancel_event)
            meta, arrays = _encoding_arrays(surface, grid, "aware_uniform", quant_step, True)
            encoded = Encoded(pack(meta, arrays))
            setting = {"target_cells": target_cells, "quant_step": quant_step}
            attempts.append({**setting, "encoded_bytes": encoded.nbytes})
            if encoded.nbytes <= budget:
                decoded = decode_uv(encoded.blob)
                report = {
                    "representation": "residue_aware_uv", "layout": "aware_uniform",
                    "source_bytes": source_bytes, "budget_ratio": float(budget_ratio), "budget_bytes": budget,
                    "encoded_bytes": encoded.nbytes, "setting": setting, "attempts": attempts,
                    "decoded_cells": len(decoded.xyz), "sasa_fields": list(SASA_FIELDS),
                    "layout_settings": {"fragmentation_weight": 20.0, "distortion_bound": 4.1,
                                        "timeout_seconds_per_patch": float(timeout)},
                }
                return decoded, encoded, report
    minimum = min(row["encoded_bytes"] for row in attempts)
    raise ValueError(f"UV encoding budget is infeasible: smallest stream {minimum} bytes exceeds {budget} bytes")


def prepare_spatial_descriptors(geometry, *, sasa_residues=None, budget_ratio=0.25,
                                optcuts_bin=None, timeout=120.0, cancel_event=None, log=None):
    """Prepare one decoded carrier and query cache for a WT interface."""
    if log:
        log("Preparing the shared surface and complete-partner SASA table")
    source = prepare_spatial_source(geometry, sasa_residues=sasa_residues, cancel_event=cancel_event)
    surface, encoded, report = encode_spatial_source(source, budget_ratio=budget_ratio,
                                                    optcuts_bin=optcuts_bin, timeout=timeout,
                                                    cancel_event=cancel_event, log=log)
    _check_cancel(cancel_event)
    return SpatialDescriptorCache(surface), encoded, report
