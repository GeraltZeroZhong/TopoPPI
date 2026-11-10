"""Backend-neutral caller routing and saved annotations use synthetic inputs."""
import json
import logging
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import matplotlib.pyplot as plt
import numpy as np
import pytest
import trimesh

from topoppi import cli, plugins
from topoppi.atlas.uv import set_uv_layout
from topoppi.config import DEFAULT_RUN_CONFIG
from topoppi.gui_app.forms import parse_single_run_form
from topoppi.gui_app.ui_mixin import interaction_sources
from topoppi.gui_app.workflow_mixin import WorkflowMixin
from topoppi.io.io_loader import PDBLoader
from topoppi.pipeline import _build_interaction_partner_map, _resolve_prolif_file
from topoppi.visualization.atlas_io import load_atlas, save_atlas
from topoppi.visualization.visualizer import InterfaceVisualizer

FIXTURE = Path(__file__).parent / "fixtures/tiny_complex.pdb"


@pytest.mark.parametrize("flag", ["--interactions", "--prolif"])
def test_cli_accepts_backend_names_and_both_document_flags(monkeypatch, flag):
    run = Mock()
    monkeypatch.setattr(cli, "run_interface_mapping", run)
    assert cli.main([str(FIXTURE), "--interaction-source", "toy", flag, "toy.json"]) == 0
    config = run.call_args.args[0]
    assert config.interaction_source == "toy" and config.prolif_file == "toy.json"


def test_config_and_form_accept_future_backend_before_discovery():
    config = replace(DEFAULT_RUN_CONFIG, pdb_file=str(FIXTURE), interaction_source="toy")
    config.validate()
    form = parse_single_run_form({"path": str(FIXTURE), "chain_a": "A", "chain_b": "B",
                                  "interaction_source": "toy", "res": "1", "sigma": "1"})
    assert form.to_config().interaction_source == "toy"


def test_gui_discovers_names_without_loading_backend_implementations(monkeypatch):
    from topoppi_chemistry import api
    from topoppi_chemistry.plugin import plugin

    entry = SimpleNamespace(name="toy", load=Mock(side_effect=AssertionError("Eager backend import")))
    monkeypatch.setattr(api, "_backend_entry_points", lambda: [entry])
    registry = plugins.PluginRegistry(entry_points=[])
    registry.register_plugin(plugin)
    monkeypatch.setattr(plugins, "_default_registry", registry)
    assert interaction_sources() == ("toy", "geometric")
    entry.load.assert_not_called()


@pytest.mark.parametrize("backend", ["prolif", "toy"])
def test_pipeline_and_gui_forward_selected_atoms_to_backend(monkeypatch, backend):
    selected = (np.zeros((0, 3)), [], np.zeros((0, 3)), [])
    symbol = "generate_prolif_interactions" if backend == "prolif" else "generate_interactions"
    generate = Mock(return_value="annotations.json")
    monkeypatch.setattr(f"topoppi.pipeline.{symbol}", generate)
    config = replace(DEFAULT_RUN_CONFIG, interaction_source=backend)
    assert _resolve_prolif_file(config, logging.getLogger(), input_sha256="fixture", chain_data=selected) == "annotations.json"
    assert generate.call_args.kwargs["chain_data"] is selected
    assert generate.call_args.kwargs.get("backend", "prolif") == backend

    monkeypatch.setattr(f"topoppi.gui_app.workflow_mixin.{symbol}", generate)
    host = SimpleNamespace(log=lambda _: None)
    result = WorkflowMixin.resolve_prolif_interactions(host, FIXTURE, "A", "B", backend=backend, chain_data=selected)
    assert result == "annotations.json"
    assert generate.call_args.kwargs["chain_data"] is selected
    assert generate.call_args.kwargs.get("backend", "prolif") == backend


def test_recorded_backend_and_canonical_types_render_without_plugin(tmp_path, monkeypatch):
    document = tmp_path / "toy.json"
    document.write_text(json.dumps({"engine": "toy", "chain_a": "A", "chain_b": "B",
        "interactions": [{"res_a_seq": "1", "res_b_seq": "2", "interaction": "PolarContact"}]}))
    loader = PDBLoader(FIXTURE)
    coords_a, atoms_a = loader.get_chain_data("A")
    coords_b, atoms_b = loader.get_chain_data("B")
    config = replace(DEFAULT_RUN_CONFIG, pdb_file=str(FIXTURE), interaction_source="prolif")
    partners, source = _build_interaction_partner_map((coords_a, atoms_a, coords_b, atoms_b), str(document),
        config, logging.getLogger(), input_sha256="")
    assert source == "toy" and partners == {"A:GLY:1": {"B:ALA:2": 1}}
    monkeypatch.setattr(plugins, "run", Mock(side_effect=AssertionError("Recorded data must render without inference")))
    viz = InterfaceVisualizer(atoms_a, coords_a, coords_b, atoms_b, chain_a_id="A", chain_b_id="B",
                              prolif_file=str(document))
    assert viz.interaction_type_source == "toy"
    assert viz.prolif_data == {"1": {"PolarContact"}}
    patch = trimesh.Trimesh(vertices=[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]], faces=[[0, 1, 2]], process=False)
    patch.metadata["source_atom_indices"] = np.asarray([0, 1, 2])
    set_uv_layout(patch, np.asarray(patch.vertices[:, :2]), key="uv_global")
    target = save_atlas(tmp_path / "saved.npz", [patch], viz)
    document.unlink()
    restored = load_atlas(target)
    assert restored.visualizer.interaction_type_source == "toy"
    assert restored.visualizer.prolif_data == {"1": {"PolarContact"}}
    figure = restored.visualizer.plot_patches(restored.patches, show=False)
    try:
        assert restored.visualizer.last_report["interaction_residue_source"] == "toy"
        assert "toy interaction record" in restored.visualizer.last_report["interaction_residue_definition"]
    finally:
        plt.close(figure)
