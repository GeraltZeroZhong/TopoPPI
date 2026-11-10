"""Analytical checks for public spatial queries and finite UV preparation."""
from __future__ import annotations

import copy
import threading

import numpy as np
import pandas as pd
import pytest

from topoppi_representation import SpatialDescriptorCache
from topoppi_representation.schema import ComplexSource, Patch
from topoppi_representation.source import native_surface
from topoppi_representation.spatial import attach_residue_sasa, encode_spatial_source


def square_source():
    xy = np.array([[0., 0.], [2., 0.], [2., 2.], [0., 2.]])
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    patches = [Patch(side, 0, np.c_[xy, np.full(4, side * 2.)], faces.copy(),
                     np.full(4, side), {"aware_uniform": xy[faces].copy()}) for side in (0, 1)]
    return ComplexSource(residue_keys=["side0:AA:VAL:7B", "side1:BB:ASP:8"],
                         amino_acids=["V", "D"], patches=patches)


def test_global_and_local_descriptors_have_physical_values_and_reuse_queries():
    cache = SpatialDescriptorCache(native_surface(square_source()))
    values, names, _ = cache.describe()
    result = dict(zip(names, values, strict=True))
    assert len(names) == 60
    assert result["partner_sum.mass"] == pytest.approx(8.)
    assert result["partner_absdiff.mass"] == pytest.approx(0.)
    assert result["partner_sum.radius_squared"] == pytest.approx(4.)
    assert result["partner_sum.L_H"] == pytest.approx(1.)
    assert result["partner_sum.AA_fraction_V"] == pytest.approx(1.)
    assert result["partner_sum.AA_fraction_D"] == pytest.approx(1.)
    values, names, quality = cache.describe([("V", "AA", 7, "B", "A")])
    result = dict(zip(names, values, strict=True))
    assert len(names) == 121
    assert result["query.delta_hydropathy"] == pytest.approx(1.8 - 4.2)
    assert result["query.mean.footprint_mass"] == pytest.approx(4.)
    assert result["query.mean.hydropathy_4A"] == pytest.approx(4.2)
    assert result["query.mean.hydropathy_12A"] == pytest.approx(4.2)
    assert quality["covered_sites"] == ["side0:AA:VAL:7B"]
    local = cache.local["side0:AA:VAL:7B"]
    cache.describe([("V", "AA", 7, "B", "G")])
    assert cache.local["side0:AA:VAL:7B"] is local


def test_missing_footprint_retains_identity_and_nan_local_information():
    source = square_source()
    source.residue_keys.append("side0:AA:LYS:99")
    source.amino_acids.append("K")
    cache = SpatialDescriptorCache(native_surface(source))
    values, names, quality = cache.describe([("K", "AA", 99, "", "A")])
    result = dict(zip(names, values, strict=True))
    assert result["query.n_mutations"] == result["query.missing_sites"] == 1
    assert result["query.wt_K"] == result["query.mt_A"] == 1
    assert result["query.covered_sites"] == 0
    assert np.isnan(result["query.mean.footprint_mass"])
    assert quality["missing_sites"] == ["side0:AA:LYS:99"]
    with pytest.raises(ValueError, match="WT residue mismatch"):
        cache.describe([("A", "AA", 7, "B", "V")])


def test_joint_mutation_order_and_rigid_translation_preserve_descriptors():
    surface = native_surface(square_source())
    moved = copy.deepcopy(surface)
    moved.xyz += [20., -7., 13.]
    first = [("V", "AA", 7, "B", "A"), ("D", "BB", 8, "", "A")]
    a, names, _ = SpatialDescriptorCache(surface).describe(first)
    b, other_names, _ = SpatialDescriptorCache(moved).describe(first[::-1])
    assert names == other_names
    np.testing.assert_allclose(a, b, atol=1e-12)


def test_complete_sasa_transport_sorts_extra_residues_and_retains_parent_aa():
    source = square_source()
    rows = []
    for side, key, aa in [(0, "AA:VAL:7B", "V"), (1, "BB:ASP:8", "D"),
                          (0, "AA:CGU:10", "E"), (0, "AA:LYS:9", "K")]:
        rows.append(dict(side=side, residue_key=key, aa=aa, n_points=960,
                         **{"sasa.isolated_partner": 10., "sasa.bound_complex": 8., "sasa.buried": 2.}))
    result = attach_residue_sasa(source, pd.DataFrame(rows))
    assert result.residue_keys == source.residue_keys + ["side0:AA:CGU:10", "side0:AA:LYS:9"]
    assert result.amino_acids == ["V", "D", "E", "K"]
    assert result.metadata["sasa_transport"]["zero_footprint_rows"] == 2
    assert result.residue_fields["sasa.buried"].tolist() == [2.] * 4
    assert source.residue_fields == {}


def test_uv_stream_is_selected_under_budget_and_decoded_without_native_fallback():
    pytest.importorskip("zstandard")
    pytest.importorskip("pywt")
    from topoppi_representation import decode_uv

    source = square_source()
    decoded, code, report = encode_spatial_source(source, budget_ratio=10.)
    assert code.nbytes <= report["budget_bytes"]
    assert report["setting"] == {"target_cells": 8, "quant_step": .001}
    np.testing.assert_array_equal(decoded.xyz, decode_uv(code.blob).xyz)
    assert decoded.faces is None
    with pytest.raises(ValueError, match="budget is infeasible"):
        encode_spatial_source(source, budget_ratio=.25)


def test_spatial_preparation_honors_cancellation_before_encoding():
    pytest.importorskip("zstandard")
    pytest.importorskip("pywt")
    stop = threading.Event()
    stop.set()
    with pytest.raises(InterruptedError):
        encode_spatial_source(square_source(), cancel_event=stop)


def test_core_adapter_canonicalizes_partner_order_without_changing_meshes():
    from dataclasses import replace
    from pathlib import Path

    import topoppi
    from topoppi_representation import prepare_spatial_source
    from topoppi_representation.descriptors.physical import residue_physics
    from topoppi_representation.source import canonical_source_bytes

    structure = topoppi.load_structure(Path(__file__).parent / "fixtures/synthetic_complex.pdb")
    geometry = topoppi.compute_geometry(topoppi.get_interface(structure, "B", "A"))
    triangles = [[p.triangles.copy() for p in partner.patches] for partner in geometry.partners]
    residues, _ = residue_physics(geometry)
    source = prepare_spatial_source(geometry, sasa_residues=residues)
    default_source = prepare_spatial_source(geometry)
    reordered = replace(geometry, partners=geometry.partners[::-1])
    other = prepare_spatial_source(reordered, sasa_residues=residues)
    assert canonical_source_bytes(source) == canonical_source_bytes(other)
    assert canonical_source_bytes(source) == canonical_source_bytes(default_source)
    assert source.metadata["chain_groups"] == [["A"], ["B"]]
    assert sum(source.residue_fields["sasa.buried"]) == pytest.approx(residues["sasa.buried"].sum())
    for partner, before in zip(geometry.partners, triangles, strict=True):
        for patch, expected in zip(partner.patches, before, strict=True):
            np.testing.assert_array_equal(patch.triangles, expected)
