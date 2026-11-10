"""Structure-to-result checks for the bundled mutation-derived ranking model."""
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from topoppi_hotspot.inference import model_requirements, rank_complex, rank_interface
from topoppi_hotspot.mutation_model import MutationSurfaceRegressor
from topoppi_hotspot.resources import default_model_dir, resolve_model_dir

FIXTURE = Path(__file__).parent / "fixtures/tiny_complex.pdb"
MODEL_ID = "topoppi-hotspot-mutation-surface"


@pytest.fixture
def surface_runtime():
    """Real parameterization needs the UV dependencies and residue-aware OptCuts."""
    pytest.importorskip("pywt", reason="Install the representation dependencies for UV inference")
    pytest.importorskip("zstandard", reason="Install the representation dependencies for UV inference")
    from topoppi.config import OptCutsConfig
    from topoppi.optimization.optcuts.joint_optimizer import (
        resolve_optcuts_binary,
        supports_residue_footprint_energy,
    )

    binary = resolve_optcuts_binary(OptCutsConfig(residue_fragmentation_weight=20.))
    if binary is None or not os.access(binary, os.X_OK) or not supports_residue_footprint_energy(binary):
        pytest.skip("Structure inference requires the residue-aware OptCuts executable")
    return binary


def test_official_model_has_complete_portable_feature_schema():
    folder = resolve_model_dir()
    assert folder == default_model_dir().resolve()
    metadata = json.loads((folder / "model.json").read_text())
    assert metadata["model_id"] == MODEL_ID
    assert metadata["columns"] == metadata["feature_names"]
    assert len(set(metadata["columns"])) == 552
    assert metadata["target"] == "ddg" and metadata["units"] == "kcal/mol"
    requirements = model_requirements()
    assert requirements["needs_uv"]
    assert not requirements["needs_embeddings"]
    assert not requirements["needs_graph"]
    with np.load(folder / "weights.npz", allow_pickle=False) as model:
        assert np.isfinite(model["dual_coefficient"]).all()
        assert model["training_vectors"].shape == (len(model["dual_coefficient"]), len(model["mean"]))
        assert np.all(model["scale"] > 0)


def test_default_structure_prediction_replays_from_saved_features(tmp_path, monkeypatch, surface_runtime):
    import socket

    from topoppi_affinity.data import parse_mutations
    from topoppi_affinity.features import extract_complex_features
    from topoppi_affinity.production import PhysicalRegressor, bundled_model, physical_feature_rows

    from topoppi import load_interface
    from topoppi_representation import SpatialDescriptorCache, decode_uv

    def no_network(*args, **kwargs):
        raise AssertionError("The bundled model must run without network access")

    monkeypatch.setattr(socket.socket, "connect", no_network)
    output = tmp_path / "ranking"
    table, metadata = rank_complex(FIXTURE, "A", "B", "A", output_dir=output, optcuts_bin=surface_runtime)
    assert len(table) == 1
    assert metadata["model_id"] == MODEL_ID and metadata["model_version"] == "3.0.0"
    assert metadata["target"] == "ddg" and metadata["units"] == "kcal/mol"
    assert metadata["embeddings"] is None
    assert metadata["candidate_coverage"]["n_native_alanine_context_nodes"] == 1
    np.testing.assert_array_equal(table.rank_score, table.predicted_ddg_kcal_mol)
    assert all((output / filename).is_file() for filename in metadata["files"].values())
    assert metadata["encoding"]["encoded_bytes"] <= metadata["encoding"]["budget_bytes"]

    features = pd.read_csv(output / "features.csv", index_col="residue_key")
    assert list(features.index) == table.residue_key.tolist()
    predictor = MutationSurfaceRegressor(resolve_model_dir() / "model.json")
    replay = predictor.predict(features.to_numpy(float), features.columns)
    np.testing.assert_allclose(replay, table.rank_score, rtol=1e-11, atol=1e-11)

    geometry = load_interface(output / "interface.npz")
    bundle = extract_complex_features(FIXTURE, "A", "B", prepared=geometry)
    mutation = PhysicalRegressor(bundled_model("mutation"))
    queries = pd.DataFrame({"mutations": table.mutations, "temperature_K": 298.15})
    physical, names = physical_feature_rows(geometry, bundle, queries, mutation)
    np.testing.assert_allclose(features[names].to_numpy(float), physical, rtol=1e-12, atol=1e-12, equal_nan=True)
    decoded = decode_uv((output / "representation.bin").read_bytes())
    spatial = SpatialDescriptorCache(decoded)
    values, names, coverage = spatial.describe(parse_mutations(table.mutations.iloc[0]))
    np.testing.assert_allclose(features[names].iloc[0], values, rtol=1e-12, atol=1e-12, equal_nan=True)
    expected_coverage = "missing" if coverage["missing_sites"] else "covered"
    assert table.spatial_coverage.iloc[0] == expected_coverage


def test_default_model_cli_json_matches_saved_ranking(tmp_path, capsys, surface_runtime):
    from topoppi_hotspot import cli

    assert cli.main(["rank-complex", "--structure", str(FIXTURE), "--chains-a", "A", "--chains-b", "B",
                     "--partner", "A", "--output-dir", str(tmp_path), "--json",
                     "--optcuts-bin", surface_runtime]) == 0
    summary = json.loads(capsys.readouterr().out)
    saved = pd.read_csv(summary["files"]["rankings"])
    assert summary["model_id"] == MODEL_ID
    assert summary["target"] == "ddg" and summary["units"] == "kcal/mol"
    assert summary["n_ranked_candidates"] == len(saved) == 1
    np.testing.assert_allclose([row["rank_score"] for row in summary["rankings"]], saved.rank_score)
    assert {"features", "representation", "interface"}.issubset(summary["files"])
    assert all(Path(path).is_file() for path in summary["files"].values())


@pytest.mark.parametrize("failure", ["partner", "missing_structure", "missing_optcuts"])
def test_cli_reports_input_and_uv_failures_without_tracebacks(failure, tmp_path, monkeypatch, capsys):
    from topoppi_hotspot import cli

    structure = tmp_path / "absent.pdb" if failure == "missing_structure" else FIXTURE
    partner = "C" if failure == "partner" else "A"
    args = ["rank-complex", "--structure", str(structure), "--chains-a", "A", "--chains-b", "B",
            "--partner", partner, "--output", str(tmp_path / "result"), "--json"]
    if failure == "missing_optcuts":
        pytest.importorskip("pywt")
        pytest.importorskip("zstandard")
        monkeypatch.delenv("TOPOPPI_OPTCUTS_BIN", raising=False)
        args += ["--optcuts-bin", str(tmp_path / "missing_OptCuts")]
    assert cli.main(args) == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.strip() and "Traceback" not in captured.err
    if failure == "partner":
        assert "exactly match one complete" in captured.err
        with pytest.raises(ValueError, match="exactly match one complete"):
            rank_complex(FIXTURE, "A", "B", "C", output_dir=tmp_path / "api")
    elif failure == "missing_structure":
        assert "absent.pdb" in captured.err
    else:
        assert "OptCuts" in captured.err
    assert not (tmp_path / "result" / "rankings.csv").exists()


def test_desktop_blank_model_saves_reopenable_geometry(tmp_path, monkeypatch, surface_runtime):
    from topoppi_hotspot.gui_workflow import RankingSettings, load_ranking, run_ranking
    from topoppi_hotspot.plugin import plugin

    from topoppi import plugins

    monkeypatch.setenv("TOPOPPI_OPTCUTS_BIN", surface_runtime)
    registry = plugins.PluginRegistry(entry_points=[])
    registry.register_plugin(plugin)
    monkeypatch.setattr(plugins, "_default_registry", registry)
    settings = RankingSettings(str(FIXTURE), "A", "B", "A", "", str(tmp_path / "saved"))
    result = run_ranking(settings)
    reopened = load_ranking(settings.output_dir)
    assert reopened.geometry is not None
    assert reopened.metadata["model_id"] == MODEL_ID
    np.testing.assert_allclose([float(x["rank_score"]) for x in reopened.rows],
                               [float(x["rank_score"]) for x in result.rows])
    _, _, style = reopened.view(reopened.rows[0]["residue_key"])
    assert style["annotation_label"] == "Predicted Ala ΔΔG (kcal/mol)"


def test_cli_preserves_an_existing_ranking(tmp_path, capsys, surface_runtime):
    from topoppi_hotspot import cli

    rank_complex(FIXTURE, "A", "B", "A", output_dir=tmp_path, optcuts_bin=surface_runtime)
    saved = {path.name: path.read_bytes() for path in tmp_path.iterdir() if path.is_file()}
    assert cli.main(["rank-complex", "--structure", str(FIXTURE), "--chains-a", "A", "--chains-b", "B",
                     "--partner", "B", "--output", str(tmp_path), "--json"]) == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "Choose a new output directory" in captured.err
    assert saved == {path.name: path.read_bytes() for path in tmp_path.iterdir() if path.is_file()}


def test_module_cli_returns_failure_status_for_invalid_partner(tmp_path):
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [str(Path(__file__).parents[1]), env.get("PYTHONPATH")]))
    completed = subprocess.run([sys.executable, "-m", "topoppi_hotspot", "rank-complex", "--structure", str(FIXTURE),
                                "--chains-a", "A", "--chains-b", "B", "--partner", "C",
                                "--output", str(tmp_path / "results"), "--json"],
                               cwd=tmp_path, env=env, capture_output=True, text=True, check=False)
    assert completed.returncode == 1
    assert completed.stdout == ""
    assert "exactly match one complete" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_cli_explains_incomplete_custom_model(tmp_path, capsys):
    from topoppi_hotspot import cli

    model = tmp_path / "unfinished-model"
    model.mkdir()
    (model / "model.json").write_text("{}")
    assert cli.main(["rank-complex", "--structure", str(FIXTURE), "--chains-a", "A", "--chains-b", "B",
                     "--partner", "A", "--model-dir", str(model), "--output", str(tmp_path / "results")]) == 1
    assert "must declare a model kind" in capsys.readouterr().err


@pytest.mark.parametrize("metadata, message", [
    (None, "must be a JSON object"),
    ([], "must be a JSON object"),
    ({"kind": "linear"}, "must declare its feature columns"),
    ({"kind": "linear", "columns": [], "feature_protocol": []}, "feature_protocol must be a JSON object"),
    ({"kind": "linear", "columns": [], "feature_protocol": {"surface": {"resoluton": .5}}},
     "Invalid model feature_protocol"),
])
def test_cli_explains_malformed_custom_model_metadata(metadata, message, tmp_path, capsys):
    from topoppi_hotspot import cli

    model = tmp_path / "custom-model"
    model.mkdir()
    (model / "model.json").write_text(json.dumps(metadata))
    output = tmp_path / "results"
    assert cli.main(["rank-complex", "--structure", str(FIXTURE), "--chains-a", "A", "--chains-b", "B",
                     "--partner", "A", "--model-dir", str(model), "--output", str(output), "--json"]) == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert message in captured.err
    assert "Traceback" not in captured.err
    assert not output.exists()


def test_cli_reports_keyboard_cancellation(tmp_path, monkeypatch, capsys):
    from topoppi_hotspot import cli, inference

    def interrupt(*args, **kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr(inference, "rank_complex", interrupt)
    assert cli.main(["rank-complex", "--structure", str(FIXTURE), "--chains-a", "A", "--chains-b", "B",
                     "--partner", "A", "--output", str(tmp_path)]) == 130
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == "Ranking cancelled.\n"


@pytest.mark.parametrize("partner", ["A", "B"])
def test_swapping_input_groups_preserves_each_partner_prediction(tmp_path, surface_runtime, partner):
    original, _ = rank_complex(FIXTURE, "A", "B", partner, output_dir=tmp_path / "original",
                               optcuts_bin=surface_runtime)
    swapped, metadata = rank_complex(FIXTURE, "B", "A", partner, output_dir=tmp_path / "swapped",
                                     optcuts_bin=surface_runtime)
    assert metadata["selected_partner"] == [partner]
    assert original.residue_key.tolist() == swapped.residue_key.tolist()
    np.testing.assert_allclose(original.rank_score, swapped.rank_score, rtol=1e-12, atol=1e-12)
    first = pd.read_csv(tmp_path / "original/features.csv", index_col="residue_key")
    second = pd.read_csv(tmp_path / "swapped/features.csv", index_col="residue_key")
    np.testing.assert_allclose(first, second, rtol=1e-12, atol=1e-12, equal_nan=True)


def test_shared_geometry_predicts_after_original_structure_is_removed(tmp_path, surface_runtime):
    from topoppi_hotspot.gui_workflow import load_ranking

    import topoppi

    source = tmp_path / "complex.pdb"
    source.write_bytes(FIXTURE.read_bytes())
    geometry = topoppi.compute_geometry(topoppi.get_interface(topoppi.load_structure(source), "A", "B"))
    source.unlink()
    output = tmp_path / "saved"
    table, _ = rank_interface(geometry, partner="A", output_dir=output, optcuts_bin=surface_runtime)
    assert np.isfinite(table.rank_score).all()
    document = load_ranking(output)
    assert document.geometry is not None
    assert document.rows[0]["residue_key"] == table.residue_key.iloc[0]
