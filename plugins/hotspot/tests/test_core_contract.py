"""Object-based inference must reuse Core geometry and residue selection."""
import json
import shutil
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from Bio.PDB import MMCIFIO, PDBParser
from topoppi_hotspot.inference import rank_complex, rank_interface
from topoppi_hotspot.sequence_embeddings import (
    prepare_chain_embedding_registry,
    prepare_interface_embedding_registry,
)
from topoppi_hotspot.solvent_accessibility import residue_sasa

import topoppi

FIXTURE = Path(__file__).parent / "fixtures/tiny_complex.pdb"


def test_shared_geometry_ranking_sasa_and_sequences_survive_source_removal(tmp_path, monkeypatch):
    path = tmp_path / "input.pdb"
    shutil.copy2(FIXTURE, path)
    model = tmp_path / "model"
    model.mkdir()
    (model / "model.json").write_text(json.dumps({
        "kind": "linear", "columns": ["surface.area", "sasa.buried"], "seeds": [0],
    }))
    np.savez(model / "linear_model.npz", mean=np.zeros(2), scale=np.ones(2), coefficient=[1., .1])
    expected, old_metadata = rank_complex(path, "A", "B", "A", model, tmp_path / "legacy")
    sasa = residue_sasa(path, "A", "B")
    geometry = topoppi.compute_geometry(topoppi.get_interface(topoppi.load_structure(path), "B", "A"))
    table = pd.DataFrame([{"complex_id": "input:A:B", "structure_path": str(path),
                           "chains_a": "A", "chains_b": "B"}])
    prepare_chain_embedding_registry(table, tmp_path / "sequences_before")
    path.unlink()

    def unexpected(*args, **kwargs):
        raise AssertionError("A supplied Core object must not be parsed or meshed again")

    monkeypatch.setattr("topoppi.io.io_loader.PDBLoader.__init__", unexpected)
    monkeypatch.setattr("topoppi.interface.prepare_partner_surface", unexpected)
    monkeypatch.setattr("topoppi.geometry.prepare_partner_surface", unexpected)
    result, metadata = rank_interface(geometry, "A", model, tmp_path / "shared")
    pd.testing.assert_frame_equal(result, expected)
    context = topoppi.ExecutionContext(plugin_config={"hotspot": {
        "partner": "B", "model_dir": "unused", "output_dir": "unused",
    }})
    overridden, _ = rank_interface(geometry, "A", model, tmp_path / "overridden", context=context)
    pd.testing.assert_frame_equal(overridden, expected)
    configured, _ = rank_interface(geometry, context=topoppi.ExecutionContext(plugin_config={"hotspot": {
        "partner": "A", "model_dir": model, "output_dir": tmp_path / "configured",
    }, "uninstalled": "opaque"}))
    pd.testing.assert_frame_equal(configured, expected)
    assert metadata["candidate_coverage"] == old_metadata["candidate_coverage"]
    assert metadata["feature_protocol"] == old_metadata["feature_protocol"]
    pd.testing.assert_frame_equal(residue_sasa(path, "B", "A", prepared=geometry), sasa)
    prepare_interface_embedding_registry(geometry, tmp_path / "sequences_after", complex_id="input:A:B")
    for filename in ("sequences.csv", "chain_sequences.csv"):
        assert (tmp_path / "sequences_before" / filename).read_text() == (
            tmp_path / "sequences_after" / filename).read_text()
    geometry.surface_config = replace(geometry.surface_config, grid_resolution=1.23)
    with pytest.raises(ValueError, match="surface/topology protocol"):
        rank_interface(geometry, "A", model, tmp_path / "mismatch")


def test_pdb_and_mmcif_sequence_registries_use_core_residue_keys(tmp_path):
    model = PDBParser(QUIET=True).get_structure("synthetic", FIXTURE)
    writer = MMCIFIO()
    writer.set_structure(model)
    cif = tmp_path / "input.cif"
    writer.save(str(cif))
    for name, path in (("pdb", FIXTURE), ("cif", cif)):
        table = pd.DataFrame([{"complex_id": "synthetic", "structure_path": str(path),
                               "chains_a": "B", "chains_b": "A"}])
        prepare_chain_embedding_registry(table, tmp_path / name)
    for filename in ("sequences.csv", "chain_sequences.csv"):
        assert (tmp_path / "pdb" / filename).read_text() == (tmp_path / "cif" / filename).read_text()
