"""The GUI uses real frozen inference and preserves independently saved results."""
import json
import shutil
from pathlib import Path
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pytest
from topoppi_hotspot.gui_workflow import RankingDocument, RankingSettings, load_ranking, run_ranking
from topoppi_hotspot.plugin import plugin

import topoppi
from topoppi import plugins

FIXTURE = Path(__file__).parent / "fixtures/tiny_complex.pdb"


def test_surface_highlights_all_selected_covered_residues_without_changing_scores(tmp_path):
    visualizer = object()
    geometry = SimpleNamespace(partners=[SimpleNamespace(chains=("A",), patches=["native patch"])],
                               visualizer=lambda _: visualizer)
    rows = [{"residue_key": "A:GLY:1", "rank_score": .9, "footprint_area_A2": 2.0},
            {"residue_key": "A:TYR:2", "rank_score": .5, "footprint_area_A2": 3.0},
            {"residue_key": "A:TRP:3", "rank_score": -.1, "footprint_area_A2": 0.0}]
    document = RankingDocument(tmp_path, {"selected_partner": ["A"]}, rows, geometry)
    _patches, actual_visualizer, style = document.view([row["residue_key"] for row in rows])
    assert actual_visualizer is visualizer
    assert style["highlight_residues"] == ("A:GLY:1", "A:TYR:2")
    assert style["annotation_values"] == {row["residue_key"]: row["rank_score"] for row in rows}
    assert document.view("A:TYR:2")[2]["highlight_residues"] == ("A:TYR:2",)


@pytest.fixture
def ranking_inputs(tmp_path, monkeypatch):
    registry = plugins.PluginRegistry(entry_points=[])
    registry.register_plugin(plugin)
    monkeypatch.setattr(plugins, "_default_registry", registry)
    source = tmp_path / "synthetic.pdb"
    shutil.copyfile(FIXTURE, source)
    model = tmp_path / "fixed-synthetic-model"
    model.mkdir()
    (model / "model.json").write_text(json.dumps({"kind": "linear", "columns": ["surface.area"], "seeds": [0],
                                                 "representation": "Synthetic GUI test"}))
    np.savez(model / "linear_model.npz", mean=[0.], scale=[1.], coefficient=[1.])
    return RankingSettings(str(source), "A", "B", "A", str(model), str(tmp_path / "ranking"))


def test_reuses_geometry_and_reopens_results_without_source_or_model(ranking_inputs, monkeypatch, tmp_path):
    settings = ranking_inputs
    geometry = topoppi.compute_geometry(topoppi.get_interface(topoppi.load_structure(settings.structure), "A", "B"))
    Path(settings.structure).unlink()

    def repeated(*args, **kwargs):
        raise AssertionError("The supplied compatible geometry must be reused")

    monkeypatch.setattr(topoppi, "load_structure", repeated)
    monkeypatch.setattr(topoppi, "compute_geometry", repeated)
    document = run_ranking(settings, prepared=geometry)
    assert document.geometry is geometry
    assert document.metadata["geometry_reused"] is True
    assert len(document.rows) == 1 and document.rows[0]["residue_key"] == "A:GLY:1"
    assert document.rows[0]["rank_score"] > 0
    shutil.rmtree(settings.model_dir)
    reopened = load_ranking(settings.output_dir)
    assert reopened.geometry.analysis_id == geometry.analysis_id
    assert float(reopened.rows[0]["rank_score"]) == document.rows[0]["rank_score"]
    patches, visualizer, style = reopened.view("A:GLY:1")
    assert patches and visualizer.chain_a_id == "A"
    assert style["annotation_values"]["A:GLY:1"] == document.rows[0]["rank_score"]
    figure = visualizer.plot_patches(patches, show=False, style_config=style)
    try:
        figure.savefig(tmp_path / "ranked-surface.png")
        assert (tmp_path / "ranked-surface.png").stat().st_size > 0
    finally:
        plt.close(figure)
    full_csv, core_csv = tmp_path / "full.csv", tmp_path / "scores.csv"
    reopened.export(full_csv)
    reopened.export(core_csv, annotations=True)
    assert "rank_score" in full_csv.read_text().splitlines()[0]
    assert core_csv.read_text().splitlines()[0] == "residue,value"
    assert "A:GLY:1" in core_csv.read_text()


def test_prepares_missing_partner_surface_from_existing_atoms(ranking_inputs, monkeypatch):
    settings = ranking_inputs
    geometry = topoppi.compute_geometry(topoppi.get_interface(topoppi.load_structure(settings.structure), "A", "B"), sides=(0,))
    Path(settings.structure).unlink()
    monkeypatch.setattr(topoppi, "load_structure", lambda *_: pytest.fail("Shared atoms should be reused"))
    document = run_ranking(settings, prepared=geometry)
    assert all(partner.mesh is not None for partner in document.geometry.partners)
    assert document.metadata["geometry_reused"] is False


def test_missing_model_and_existing_results_report_actionable_errors(ranking_inputs):
    from dataclasses import replace

    with pytest.raises(ValueError, match="Choose a model folder"):
        run_ranking(replace(ranking_inputs, model_dir=ranking_inputs.model_dir + "-missing"))
    output = Path(ranking_inputs.output_dir)
    output.mkdir()
    previous = output / "rankings.csv"
    previous.write_text("saved result")
    with pytest.raises(ValueError, match="Choose a new results folder"):
        run_ranking(ranking_inputs)
    assert previous.read_text() == "saved result"


def test_embedding_model_requires_explicit_local_registry(ranking_inputs):
    (Path(ranking_inputs.model_dir) / "embedding_pca.npz").touch()
    with pytest.raises(ValueError, match="Select its prepared embedding registry"):
        run_ranking(ranking_inputs)
    assert not Path(ranking_inputs.output_dir).exists()


@pytest.mark.parametrize("filename", ["features.csv", "candidate_coverage.csv", "representation.bin"])
def test_partial_outputs_are_preserved_before_recomputing(ranking_inputs, monkeypatch, filename):
    output = Path(ranking_inputs.output_dir)
    output.mkdir()
    previous = output / filename
    previous.write_bytes(b"Saved partial ranking data")
    monkeypatch.setattr(topoppi, "compute_geometry", lambda *_args, **_kwargs: pytest.fail("Output conflict reached geometry"))
    with pytest.raises(ValueError, match="Choose a new results folder"):
        run_ranking(ranking_inputs)
    assert previous.read_bytes() == b"Saved partial ranking data"
    assert sorted(path.name for path in output.iterdir()) == [filename]


def test_output_file_is_rejected_before_computing_geometry(ranking_inputs, monkeypatch):
    from dataclasses import replace

    source = Path(ranking_inputs.structure).read_text()
    monkeypatch.setattr(topoppi, "compute_geometry", lambda *_args, **_kwargs: pytest.fail("Validate the destination before geometry"))
    with pytest.raises(ValueError, match="Choose a writable results folder"):
        run_ranking(replace(ranking_inputs, output_dir=ranking_inputs.structure))
    assert Path(ranking_inputs.structure).read_text() == source


def test_cancelled_request_creates_no_output(ranking_inputs):
    from topoppi.errors import ExecutionCancelled

    context = topoppi.ExecutionContext()
    context.cancel_event.set()
    with pytest.raises(ExecutionCancelled):
        run_ranking(ranking_inputs, context=context)
    assert not Path(ranking_inputs.output_dir).exists()


@pytest.mark.parametrize("saved_rows", [0, 1])
def test_truncated_ranking_table_is_rejected_and_complete_copy_recovers(tmp_path, saved_rows):
    metadata = {"selected_partner": ["A"], "candidate_coverage": {"n_ranked_candidates": 2}}
    (tmp_path / "metadata.json").write_text(json.dumps(metadata))
    lines = ["residue_key,rank_score,rank_min,rank_max\n", "A:TYR:2,1.0,1,1\n", "A:GLY:3,0.5,2,2\n"]
    table = tmp_path / "rankings.csv"
    table.write_text(''.join(lines[:saved_rows + 1]))
    with pytest.raises(ValueError, match="Saved rankings table is incomplete"):
        load_ranking(tmp_path)
    table.write_text(''.join(lines))
    assert [row["residue_key"] for row in load_ranking(tmp_path).rows] == ["A:TYR:2", "A:GLY:3"]
