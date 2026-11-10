import json
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from topoppi_hotspot import cli, inference
from topoppi_hotspot.block_ranking import BlockPairwiseRanker
from topoppi_hotspot.features import FEATURE_SETS
from topoppi_hotspot.graph_operators import build_graph_operators
from topoppi_hotspot.graph_ranking import GraphRanker
from topoppi_hotspot.partner_features import feature_specifications, prepare_partner_features
from topoppi_hotspot.rank_fusion import aggregate_ranks, standard_nonalanine_candidates
from topoppi_hotspot.representations import center_within_queries

from topoppi.config import SurfaceConfig, TopologyConfig


def _toy_features():
    keys = ["A:LYS:1", "A:ALA:2", "A:VAL:3", "B:TYR:1", "B:ALA:2", "C:ARG:1", "C:MSE:2", "C:GLY:3"]
    wt = list("KAVYARMG")
    frame = pd.DataFrame({"residue_key": keys, "wt": wt, "side": [0] * 3 + [1] * 5,
                          "chain": [key.split(":")[0] for key in keys],
                          "residue_number": [int(key.split(":")[2]) for key in keys],
                          "insertion_code": "", "ca_resolved": [True] * 7 + [False]})
    for column in set(column for columns in FEATURE_SETS.values() for column in columns):
        frame[column] = 1.
    for letter in "ACDEFGHIKLMNPQRSTVWYX":
        frame[f"aa.{letter}"] = (frame.wt == letter).astype(float)
    frame["surface.area"] = [2., 10., 6., 2., 100., 6., 4., 8.]
    frame["surface.footprint_area"] = [1., 2., 3., 4., 5., 0., 3., 4.]
    frame["contact.degree"] = [3., 20., 1., 3., 100., 1., 4., 2.]
    local = np.array([[0, 1], [1, 2], [3, 4], [4, 5], [5, 6], [6, 7]])
    graph = {"residue_keys": np.asarray(keys), "footprint_edges": local, "footprint_weights": np.ones(len(local)),
             "spatial_edges": local, "cross_edges": np.array([[0, 3], [1, 4], [2, 5], [2, 6], [2, 7]]),
             "cross_weights": np.ones(5)}
    return frame, graph


def _linear_bundle(directory, columns, coefficient=None, feature_protocol=None):
    directory.mkdir(parents=True, exist_ok=True)
    metadata = {"kind": "linear", "representation": "toy", "columns": columns, "seeds": [0]}
    if feature_protocol is not None:
        metadata["feature_protocol"] = feature_protocol
    (directory / "model.json").write_text(json.dumps(metadata))
    coefficient = np.ones(len(columns)) if coefficient is None else np.asarray(coefficient)
    np.savez(directory / "linear_model.npz", mean=np.zeros(len(coefficient)), scale=np.ones(len(coefficient)),
             coefficient=coefficient)
    return directory


@pytest.fixture
def mocked_structure(monkeypatch):
    features, graph = _toy_features()
    calls = []
    prepared = SimpleNamespace(metadata=lambda: {"source": "synthetic_test"})

    def prepare(structure, first, second, **kwargs):
        calls.append((first, second, kwargs))
        return prepared

    def extract(structure, first, second, *, prepared):
        assert prepared is mocked.prepared
        return features.copy(), graph

    mocked = SimpleNamespace(features=features, graph=graph, calls=calls, prepared=prepared)
    monkeypatch.setattr(inference, "prepare_interface", prepare)
    monkeypatch.setattr(inference, "extract_features", extract)
    return mocked


def test_standalone_full_context_partner_swap_configs_and_candidate_flags(tmp_path, mocked_structure):
    protocol = {"surface": asdict(SurfaceConfig(grid_resolution=1.2)), "topology": asdict(TopologyConfig())}
    model = _linear_bundle(tmp_path / "model", ["relative.partner.zscore.surface.area"], feature_protocol=protocol)
    first, metadata = inference.rank_complex("toy.pdb", "A", "B,C", "C,B", model, tmp_path / "first")
    second, _ = inference.rank_complex("toy.pdb", "C,B", "A", "B,C", model, tmp_path / "second")
    pd.testing.assert_frame_equal(first, second)
    assert len(mocked_structure.calls) == 2  # Exactly one preparation per call.
    assert mocked_structure.calls[0][:2] == (("A",), ("B", "C"))
    assert mocked_structure.calls[0][2]["surface"].grid_resolution == 1.2
    assert set(first.residue_key) == {"B:TYR:1", "C:ARG:1", "C:GLY:3"}
    part = mocked_structure.features[mocked_structure.features.side == 1]
    expected = (part["surface.area"] - part["surface.area"].mean()) / part["surface.area"].std(ddof=0)
    expected.index = part.residue_key
    np.testing.assert_allclose(first.rank_score, expected.loc[first.residue_key])
    assert first.set_index("residue_key").loc["C:ARG:1", "flags"] == "missing_native_footprint"
    assert first.set_index("residue_key").loc["C:GLY:3", "flags"] == "missing_CA"
    np.testing.assert_allclose(first.full_candidate_percentile, [5 / 6, 3 / 6, 1 / 6])
    assert metadata["candidate_coverage"]["n_WT_interface_context_nodes"] == 5
    assert metadata["candidate_coverage"]["n_native_alanine_context_nodes"] == 1
    assert metadata["candidate_coverage"]["n_other_excluded_interface_nodes"] == 1
    assert metadata["labels_used"] is False
    assert metadata["affinity_plugin_used"] is False
    assert (tmp_path / "first" / "rankings.csv").exists()
    assert (tmp_path / "first" / "candidate_coverage.csv").exists()
    with pytest.raises(ValueError, match="exactly match one complete"):
        inference.rank_complex("toy.pdb", "A", "B,C", "B", model, tmp_path / "invalid")
    assert len(mocked_structure.calls) == 2


def test_only_required_auxiliary_block_is_computed_with_reused_preparation(tmp_path, mocked_structure, monkeypatch):
    from topoppi_hotspot import (
        footprint_diffusion,
        mesh_diffusion,
        paired_surface_features,
        solvent_accessibility,
        surface_rim,
    )

    calls = []

    def paired(structure, first, second, *, prepared):
        assert prepared is mocked_structure.prepared
        calls.append(prepared)
        frame = mocked_structure.features[["residue_key", "side"]].copy()
        frame["surface_pair.mean.gap"] = np.arange(len(frame), dtype=float)
        return frame, {"synthetic": True}

    def unexpected(*args, **kwargs):
        raise AssertionError("Unrequested auxiliary extractor ran")

    monkeypatch.setattr(paired_surface_features, "extract_paired_surface", paired)
    monkeypatch.setattr(mesh_diffusion, "extract_mesh_diffusion", unexpected)
    monkeypatch.setattr(surface_rim, "extract_surface_rim", unexpected)
    monkeypatch.setattr(solvent_accessibility, "residue_sasa", unexpected)
    monkeypatch.setattr(footprint_diffusion, "add_graph_diffusion", unexpected)
    model = _linear_bundle(tmp_path / "model", ["new_relative.z.surface_pair.mean.gap"])
    result, metadata = inference.rank_complex("toy.pdb", "A", "B,C", "A", model, tmp_path / "output")
    assert len(calls) == 1
    assert set(metadata["auxiliary_features"]) == {"paired_surface_features"}
    assert metadata["feature_columns"] == ["new_relative.z.surface_pair.mean.gap"]
    assert len(result) == 2


def test_fusion_uses_all_standard_nonalanine_candidates_without_labels(tmp_path):
    all_features, _ = _toy_features()
    all_features["complex_id"] = "toy"
    frame = prepare_partner_features(all_features)
    root = tmp_path / "fusion"
    _linear_bundle(root / "components" / "area", ["surface.area"])
    _linear_bundle(root / "components" / "degree", ["contact.degree"])
    (root / "model.json").write_text(json.dumps({"kind": "fusion", "components": [
        {"name": "area", "path": "components/area"}, {"name": "degree", "path": "components/degree"},
    ]}))
    scores, components = inference.score_model(root, frame)
    eligible = standard_nonalanine_candidates(frame)
    expected = aggregate_ranks([frame["surface.area"].to_numpy(), frame["contact.degree"].to_numpy()],
                               frame.query_id.to_numpy(), eligible)
    np.testing.assert_allclose(scores, expected, equal_nan=True)
    assert set(components) == {"component_rank.area", "component_rank.degree"}
    np.testing.assert_allclose(scores, np.mean(list(components.values()), axis=0), equal_nan=True)
    assert np.isnan(scores[~eligible]).all()
    assert np.isfinite(scores[eligible]).all()
    # The same function accepts numeric features without an outcome table.
    poisoned = frame.assign(hotspot="unused", Er="unused")
    np.testing.assert_allclose(inference.score_model(root, poisoned)[0], scores, equal_nan=True)
    assert inference.model_requirements(root)["columns"] == ["contact.degree", "surface.area"]


def test_auxiliary_prefixes_dispatch_once_and_normalize_the_full_interface(mocked_structure, monkeypatch):
    from topoppi_hotspot import footprint_diffusion, mesh_diffusion, solvent_accessibility, surface_rim

    calls = []

    def auxiliary(name, column):
        def extract(structure, first, second, *, prepared):
            assert prepared is mocked_structure.prepared
            calls.append(name)
            frame = mocked_structure.features[["residue_key", "side"]].copy()
            frame[column] = np.arange(len(frame), dtype=float)
            if name == "mesh":
                frame["mesh_coverage.missing_footprint"] = 0.
            return frame, {"synthetic": name}
        return extract

    def diffusion(features, graphs, *, kinds):
        calls.append("diffusion")
        assert kinds == ("native",)
        assert len(features) == 8 and set(graphs) == {"toy"}
        result = features.copy()
        result["diffusion.native.t0p25.chemical.hydropathy"] = np.arange(len(result), dtype=float)
        return result, {}, {"synthetic": True}

    def sasa(structure, first, second, *, prepared):
        assert prepared is mocked_structure.prepared
        calls.append("sasa")
        result = mocked_structure.features[["residue_key", "side"]].copy()
        result["sasa.buried"] = np.arange(len(result), dtype=float)
        return result

    monkeypatch.setattr(mesh_diffusion, "extract_mesh_diffusion", auxiliary("mesh", "mesh_diffusion.t0.mean.test"))
    monkeypatch.setattr(surface_rim, "extract_surface_rim", auxiliary("rim", "surface_rim.mean.geodesic"))
    monkeypatch.setattr(footprint_diffusion, "add_graph_diffusion", diffusion)
    monkeypatch.setattr(solvent_accessibility, "residue_sasa", sasa)
    columns = ["new_relative.z.mesh_diffusion.t0.mean.test", "mesh_coverage.missing_footprint",
               "new_relative.pct.surface_rim.mean.geodesic", "new_relative.pct.sasa.buried",
               "new_relative.z.diffusion.native.t0p25.chemical.hydropathy",
               "new_relative.z.contact_context.mean.chemical.hydropathy"]
    _, frame, _, reports = inference._prepare_features("toy.pdb", (("A",), ("B", "C")), "toy",
                                                      mocked_structure.prepared, {"columns": columns})
    assert calls == ["mesh", "rim", "diffusion", "sasa"]
    assert len(frame) == 8
    assert np.isfinite(frame[columns].to_numpy()).all()
    assert set(reports) == {"contact_context", "mesh_diffusion", "surface_rim", "footprint_diffusion", "sasa"}
    # Selected partner's percentile denominator includes native Ala and MSE,
    # because feature context is broader than the exported candidate set.
    np.testing.assert_allclose(np.sort(frame.loc[frame.side == 1, "new_relative.pct.sasa.buried"]),
                               [.1, .3, .5, .7, .9])


def test_graph_portable_model_scores_full_context_before_candidate_filtering(tmp_path):
    pytest.importorskip("torch")
    features, graph = _toy_features()
    features["complex_id"] = "toy"
    frame = prepare_partner_features(features)
    operators, _ = build_graph_operators(frame, {"toy": graph})
    rows = np.flatnonzero(standard_nonalanine_candidates(frame))
    y = np.full(len(frame), np.nan)
    y[rows] = (frame.iloc[rows]["surface.area"] > 3.).astype(int)
    X = frame[["surface.area"]].to_numpy()
    model = GraphRanker(epochs=5, mode="local").fit(X, y, frame.query_id.to_numpy(), rows,
                                                   operators["native"], operators["context"])
    model.save(tmp_path / "seed_0")
    (tmp_path / "model.json").write_text(json.dumps({"kind": "graph", "columns": ["surface.area"],
                                                    "operator_kind": "native", "mode": "local", "seeds": [0]}))
    scores, components = inference.score_model(tmp_path, frame, {"toy": graph})
    np.testing.assert_array_equal(scores, model.decision_function(X, operators["native"], operators["context"]))
    assert scores.shape == (len(frame),)
    assert components == {}


@pytest.mark.parametrize("control", ["native_source_shuffle", "native_append_shuffle"])
def test_frozen_native_control_scores_full_WT_context_before_candidate_filtering(
    tmp_path, mocked_structure, monkeypatch, control,
):
    from topoppi_hotspot import feature_controls
    from topoppi_hotspot.features import NEIGHBOR_FEATURES, neighbor_means
    from topoppi_hotspot.model_io import predict_frozen

    features, graph = mocked_structure.features, mocked_structure.graph
    graph["native_neighbor_edges"] = graph["footprint_edges"]
    means, counts = neighbor_means(features[NEIGHBOR_FEATURES], graph["native_neighbor_edges"])
    features["topology.neighbor_count"] = counts
    features.loc[:, feature_controls.native_raw_columns()[1:]] = means
    native = feature_controls.native_raw_columns()
    columns = feature_specifications()["relative_spatial"] + native + [
        f"relative.partner.{role}.{column}" for column in native for role in ("zscore", "percentile")]
    coefficient = np.zeros(len(columns))
    coefficient[columns.index("topology.mean.surface.area")] = 1.
    coefficient[columns.index("relative.partner.percentile.topology.mean.surface.area")] = .2
    model = _linear_bundle(tmp_path / "model", columns, coefficient)
    metadata = json.loads((model / "model.json").read_text())
    metadata.update(incremental_control=control, incremental_control_seed=3)
    (model / "model.json").write_text(json.dumps(metadata))
    apply_control = feature_controls.control_features
    populations = []

    def checked_control(frame, graphs, kind, seed):
        populations.append(len(frame))
        assert set(frame.side) == {0, 1}
        assert frame.wt.eq("A").sum() == 2
        assert "C:MSE:2" in set(frame.residue_key)
        assert "hotspot" not in frame
        return apply_control(frame, graphs, kind, seed)

    def unexpected_operators(*args, **kwargs):
        raise AssertionError("Linear controls do not require graph-ranker operators")

    monkeypatch.setattr(feature_controls, "control_features", checked_control)
    monkeypatch.setattr(inference, "build_graph_operators", unexpected_operators)
    result, report = inference.rank_complex("toy.pdb", "A", "B,C", "B,C", model, tmp_path / "output")
    frame = prepare_partner_features(features.assign(complex_id=report["complex_id"]))
    transformed = apply_control(frame, {report["complex_id"]: graph}, control, 3)
    expected = pd.Series(predict_frozen(model, transformed), index=frame.residue_key)
    np.testing.assert_allclose(result.rank_score, expected.loc[result.residue_key])
    original = pd.Series(predict_frozen(model, frame), index=frame.residue_key)
    assert not np.allclose(result.rank_score, original.loc[result.residue_key])
    assert populations == [8]
    assert set(result.residue_key) == {"B:TYR:1", "C:ARG:1", "C:GLY:3"}
    assert report["candidate_coverage"]["n_WT_interface_context_nodes"] == 5
    assert inference.model_requirements(model)["needs_graph"] == (control == "native_source_shuffle")
    # A preceding complex must not consume RNG state for this WT query.
    extra = prepare_partner_features(features.assign(complex_id="aaa-extra"))
    batch = pd.concat([extra, frame], ignore_index=True)
    monkeypatch.setattr(feature_controls, "control_features", apply_control)
    batch_scores, _ = inference.score_model(model, batch, {"aaa-extra": graph, report["complex_id"]: graph})
    np.testing.assert_array_equal(batch_scores[len(extra):], expected.to_numpy())


def test_block_model_query_centers_full_raw_embeddings_before_scoring(tmp_path):
    features, _ = _toy_features()
    features["complex_id"] = "toy"
    frame = prepare_partner_features(features)
    raw = np.column_stack([np.arange(len(frame), dtype=float), np.arange(len(frame), dtype=float) ** 2])
    centered = center_within_queries(raw, frame.query_id.to_numpy())
    rows = np.flatnonzero(standard_nonalanine_candidates(frame))
    y = (frame.iloc[rows]["surface.area"] > 3.).astype(int)
    geometry = frame[["surface.area"]].to_numpy()
    model = BlockPairwiseRanker().fit(geometry[rows], centered[rows], y, frame.iloc[rows].query_id.to_numpy())
    model.save(tmp_path / "estimator")
    (tmp_path / "model.json").write_text(json.dumps({"kind": "block_linear", "columns": ["surface.area"],
                                                    "model_directory": "estimator", "embedding_dimension": 2}))
    actual, _ = inference.score_model(tmp_path, frame, embeddings=raw)
    np.testing.assert_array_equal(actual, model.decision_function(geometry, centered))
    assert inference.model_requirements(tmp_path)["needs_embeddings"] is True


def _write_registry(directory, features):
    directory = Path(directory)
    (directory / "sequence_embeddings").mkdir(parents=True, exist_ok=True)
    chains, sequences = [], []
    for index, (chain, part) in enumerate(features.groupby("chain", sort=False)):
        sequence_id = f"chain_{index}"
        sequence = "".join(part.wt)
        keys = part.residue_key.tolist()
        chains.append({"complex_id": "registry_name", "chain": chain, "sequence_id": sequence_id,
                       "residue_keys": json.dumps(keys)})
        sequences.append({"sequence_id": sequence_id, "sequence": sequence})
        values = np.column_stack([part["surface.area"].to_numpy(), part["contact.degree"].to_numpy()])
        np.savez(directory / "sequence_embeddings" / f"{sequence_id}.npz", sequence=sequence, embedding=values)
    pd.DataFrame(chains).to_csv(directory / "chain_sequences.csv", index=False)
    pd.DataFrame(sequences).to_csv(directory / "sequences.csv", index=False)


def test_registry_matches_complete_WT_sequences_and_keys_independent_of_complex_name(tmp_path):
    features, _ = _toy_features()
    features["ca_resolved"] = True
    features["complex_id"] = "current_name"
    frame = prepare_partner_features(features)
    _write_registry(tmp_path, features)
    embeddings, report = inference._registry_embeddings(tmp_path, frame, features)
    np.testing.assert_array_equal(embeddings, frame[["surface.area", "contact.degree"]].to_numpy())
    assert report["registry_complex_id"] == "registry_name"
    changed = features.copy()
    changed.loc[0, "wt"] = "V"
    with pytest.raises(ValueError, match="complete observed WT chain"):
        inference._registry_embeddings(tmp_path, frame, changed)


def test_fresh_embedding_path_requires_checkpoint_and_forwards_device(tmp_path, mocked_structure, monkeypatch):
    from topoppi_hotspot import sequence_embeddings

    mocked_structure.features["ca_resolved"] = True
    model = _linear_bundle(tmp_path / "model", ["surface.area"], [0., 1.])
    np.savez(model / "embedding_pca.npz", mean=np.zeros(2), components=np.array([[1., 0.]]))
    with pytest.raises(ValueError, match="explicit --esmc-checkpoint"):
        inference.rank_complex("toy.pdb", "A", "B,C", "A", model, tmp_path / "missing")
    assert not mocked_structure.calls
    calls = []

    def prepare_registry(prepared, directory, *, complex_id):
        assert prepared is mocked_structure.prepared
        assert complex_id == "toy:A:B,C"

    def embed(directory, checkpoint, *, device):
        calls.append((checkpoint, device))
        _write_registry(directory, mocked_structure.features)

    monkeypatch.setattr(sequence_embeddings, "prepare_interface_embedding_registry", prepare_registry)
    monkeypatch.setattr(sequence_embeddings, "embed_sequences", embed)
    result, metadata = inference.rank_complex("toy.pdb", "A", "B,C", "A", model, tmp_path / "fresh",
                                              esmc_checkpoint="explicit.pt", device="cpu")
    assert calls == [("explicit.pt", "cpu")]
    assert len(result) == 2
    assert metadata["embeddings"]["registry_complex_id"] == "registry_name"


def test_fusion_rejects_incompatible_feature_protocols(tmp_path):
    root = tmp_path / "fusion"
    _linear_bundle(root / "child", ["surface.area"],
                   feature_protocol={"surface": asdict(SurfaceConfig(grid_resolution=1.2))})
    (root / "model.json").write_text(json.dumps({"kind": "fusion", "components": [{"name": "child", "path": "child"}]}))
    with pytest.raises(ValueError, match="different structure feature protocols"):
        inference.model_requirements(root)


@pytest.mark.parametrize("output_flag", ["--output-dir", "--output"])
def test_rank_complex_cli_routes_explicit_inputs(tmp_path, monkeypatch, output_flag):
    calls = []

    def rank(*args, **kwargs):
        calls.append((args, kwargs))
        return pd.DataFrame({"residue_key": ["A:LYS:1"]}), {}

    monkeypatch.setattr(inference, "rank_complex", rank)
    cli.main(["rank-complex", "--structure", "new.pdb", "--chains-a", "A,B", "--chains-b", "C",
              "--partner", "A,B", "--model-dir", "frozen", output_flag, str(tmp_path),
              "--embedding-registry", "registry", "--device", "cpu"])
    assert calls[0][0] == ("new.pdb", "A,B", "C", "A,B", "frozen", str(tmp_path))
    assert calls[0][1] == {"embedding_registry": "registry", "esmc_checkpoint": None, "device": "cpu", "optcuts_bin": None}


def test_rank_complex_json_keeps_progress_on_stderr(tmp_path, monkeypatch, capsys):
    def rank(*args, **kwargs):
        print("Computing interface features")
        return pd.DataFrame({"residue_key": ["A:LYS:1"], "rank_score": [1.5], "optional": [np.nan]}), {
            "model_id": "test-model", "model_version": "3.0.0", "model_scope": "Test scope",
            "partner_groups": ["A", "B"], "selected_partner": "A", "score_type": "within-partner ranking score",
        }

    monkeypatch.setattr(inference, "rank_complex", rank)
    assert cli.main(["rank-complex", "--structure", "new.pdb", "--chains-a", "A", "--chains-b", "B",
                     "--partner", "A", "--output", str(tmp_path), "--json"]) == 0
    captured = capsys.readouterr()
    summary = json.loads(captured.out)
    assert "Computing interface features" in captured.err
    assert summary["rankings"] == [{"residue_key": "A:LYS:1", "rank_score": 1.5, "optional": None}]
    assert summary["model_id"] == "test-model"
    assert summary["n_ranked_candidates"] == 1
    assert summary["prediction_status"] == "ok"
    assert summary["score_type"] == "within-partner ranking score"
    assert summary["files"]["rankings"] == str(tmp_path / "rankings.csv")
