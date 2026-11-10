"""Structure inference must reproduce the feature protocol used to fit custom models."""

from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import pytest
from topoppi.interface import prepare_interface

from topoppi_affinity.embedding_baselines import fit_frozen_regressor, pooled_features, zero_surface_slots
from topoppi_affinity.features import AA, extract_complex_features
from topoppi_affinity.graph_dataset import GraphDataset
from topoppi_affinity.inference import build_features, model_info, predict_features


@pytest.mark.parametrize("target", ["pkd", "ddg"])
@pytest.mark.parametrize("representation", ["plm_only", "plm_contact", "plm_surface", "plm_analytic"])
def test_custom_model_structure_features_match_training_representation(tmp_path, monkeypatch, target, representation):
    structure = Path(__file__).parent / "fixtures/tiny_complex.pdb"
    prepared = prepare_interface(structure, "A", "B", source="user", complex_id="pair", log=lambda *_: None)
    descriptors = extract_complex_features(structure, "A", "B", prepared=prepared)
    frame = pd.DataFrame([dict(sample_id="input", source="user", complex_id="pair", structure_path=str(structure),
        chains_a="A", chains_b="B", mutations="A:G1V" if target == "ddg" else "", **{target: 1.})])
    frame.to_csv(tmp_path / "inputs.csv", index=False)
    bundle = {"model": None, "feature_names": np.array([]), "report": {"target": target, "representation": representation}}

    def synthetic_embeddings(sequence_dir, **_):
        # Deterministic sequence-dependent values exercise WT/MT pooling without
        # loading a language model or contacting its upstream host.
        sequence_dir = Path(sequence_dir)
        destination = sequence_dir / "embeddings"
        destination.mkdir(exist_ok=True)
        for row in pd.read_csv(sequence_dir / "sequences.csv").itertuples():
            identity = np.array([AA.index(letter) + 1 for letter in row.sequence], dtype=np.float32)
            values = identity[:, None] * np.linspace(.01, 1., 1152, dtype=np.float32)[None, :]
            np.savez_compressed(destination / f"{row.sequence_id}.npz", embedding=values,
                                sequence=row.sequence, model_name="esmc_600m")
        return {"model": "synthetic test values"}

    monkeypatch.setattr("topoppi_affinity.sequence_embeddings.extract_embeddings", synthetic_embeddings)
    paths, protocol = build_features(prepared, descriptors, frame, bundle, tmp_path / "features", log=lambda *_: None)
    with np.load(paths["features"], allow_pickle=False) as source:
        actual, names = source["X"], source["feature_names"]
        curvature = str(source["curvature_protocol"])
    assert curvature == ("analytic_gaussian" if representation == "plm_analytic" else "angle_defect")
    assert protocol["curvature"] == curvature

    # This is the existing training-data route, independent of _pooled_vectors.
    dataset = GraphDataset(frame, tmp_path / "features/graphs.csv", tmp_path / "features/sequences",
                           use_surface=True, use_plm=True, max_nodes=None, target=target)
    vector, expected_names = pooled_features(dataset[0], target=target, use_surface=representation != "plm_only")
    expected = vector[None, :]
    if representation == "plm_contact":
        expected = zero_surface_slots(expected, expected_names)
        assert np.count_nonzero(actual[:, np.char.startswith(names, "surface.")]) == 0
    if representation == "plm_only":
        assert not np.char.startswith(names, "surface.").any()
    np.testing.assert_array_equal(names, expected_names)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=0)

    # A fitted RBF is sensitive to accidental nonzero surface slots even when
    # names/dimensions match. Verify the complete saved-model prediction route.
    training = np.vstack([expected * factor for factor in (.8, 1., 1.2, 1.4)])
    labels = [5., 6., 7., 8.] if target == "pkd" else [-2., -1., 1., 2.]
    model = fit_frozen_regressor(training, labels, ["a", "b", "c", "d"], family="rbf",
        params={"alpha": .1, "gamma": 1. / training.shape[1]}, antisymmetric=target == "ddg")
    bundle.update(model=model, feature_names=expected_names)
    path = tmp_path / "custom.joblib"
    joblib.dump(bundle, path)
    predicted = predict_features(path, **paths)
    np.testing.assert_allclose(predicted.predicted, model.predict(expected), rtol=1e-12, atol=1e-12)
    assert predicted.status.tolist() == ["ok"]


def test_public_analytic_representation_keeps_its_curvature_protocol():
    bundle = {"model": None, "feature_names": [], "report": {
        "target": "ddg", "representation": "sequence_and_analytic_surface"}}
    assert model_info(bundle)["curvature"] == "analytic_gaussian"
