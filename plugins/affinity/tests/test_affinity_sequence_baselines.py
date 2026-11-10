"""Strict sequence isolation, joint mutations and train-only regression behavior."""
import json

import numpy as np
import pandas as pd
import pytest

from topoppi_affinity.sequence_baselines import (
    SequenceFeatureStore,
    fit_sequence_regressor,
    mutation_context_features,
    prepare_sequence_features,
)


def sequence_fixture(tmp_path):
    folder = tmp_path / "sequences"
    (folder / "embeddings").mkdir(parents=True)
    sequences = ["ACDE", "GG", "KWV", "AYDE", "ACNE", "AYNE"]
    rng = np.random.default_rng(7)
    for index, sequence in enumerate(sequences):
        np.savez(folder / "embeddings" / f"seq{index}.npz", sequence=sequence,
                 embedding=rng.normal(size=(len(sequence), 5)), model_name="test_frozen")
    pd.DataFrame({"sequence_id": [f"seq{i}" for i in range(len(sequences))], "sequence": sequences}).to_csv(
        folder / "sequences.csv", index=False)
    records = []
    for chain, index, keys in [("A", 0, ["A:ALA:1", "A:CYS:2", "A:ASP:3A", "A:GLU:4"]),
                               ("B", 1, ["B:GLY:1", "B:GLY:2"]),
                               ("C", 2, ["C:LYS:1", "C:TRP:2", "C:VAL:3"])]:
        records.append({"source": "test", "complex_id": "complex", "chain": chain,
                        "sequence_id": f"seq{index}", "residue_keys": json.dumps(keys)})
    pd.DataFrame(records).to_csv(folder / "wild_type_chains.csv", index=False)
    pd.DataFrame([{"sample_id": sample, "complex_id": "complex", "chain": "A", "sequence_id": f"seq{index}"}
                  for sample, index in [("single1", 3), ("single2", 4), ("double", 5)]]).to_csv(
        folder / "mutant_chains.csv", index=False)
    frame = pd.DataFrame([{"sample_id": sample, "source": "test", "complex_id": "complex",
                           "chains_a": "A,C", "chains_b": "B", "mutations": mutations,
                           "structure_path": "/never/open/this.pdb", "ddg": "not_a_label"}
                          for sample, mutations in [("single1", "CA2Y"), ("single2", "DA3AN"),
                                                    ("double", "CA2Y,DA3AN"), ("identity", "CA2C")]])
    return folder, frame


def test_sequence_readout_needs_no_structure_masks_or_labels_and_is_partner_symmetric(tmp_path):
    folder, frame = sequence_fixture(tmp_path)
    store = SequenceFeatureStore(folder)
    row = frame.iloc[0].to_dict()
    original = store.features(row, task="absolute")["global"]
    swapped = store.features({**row, "chains_a": "B", "chains_b": "C,A"}, task="absolute")["global"]
    np.testing.assert_array_equal(original[0], swapped[0])
    assert all(name.startswith("sequence.") for name in original[1])
    report = prepare_sequence_features(frame, folder, tmp_path / "features", task="absolute")
    assert report["n_success"] == 4
    assert report["observed_affinity_labels_loaded"] is False
    assert report["coordinate_files_read"] is False
    assert "not complete experimental construct" in report["input_scope"]
    with np.load(report["features"]["global"]["path"]) as result:
        assert set(result.files) == {"X", "sample_ids", "feature_names"}
        assert len(result["X"]) == 4


def test_joint_mutant_uses_joint_chain_embedding_and_identity_encodes_zero(tmp_path):
    folder, frame = sequence_fixture(tmp_path)
    store = SequenceFeatureStore(folder)
    values = [store.features(row, task="mutation") for row in frame.to_dict("records")]
    for view in ("global", "context"):
        assert not np.allclose(values[2][view][0], values[0][view][0] + values[1][view][0])
        np.testing.assert_array_equal(values[3][view][0], 0.)
    # Changed author insertion code is mapped to its observed-sequence position.
    assert np.linalg.norm(values[1]["context"][0]) > 0


def test_mutation_context_is_odd_under_joint_wt_mt_reversal(tmp_path):
    folder, _ = sequence_fixture(tmp_path)
    store = SequenceFeatureStore(folder)
    wild, mutant = {"A": store._embedding("seq0")}, {"A": store._embedding("seq5")}
    forward, names = mutation_context_features(wild, mutant, [("A", 1), ("A", 2)])
    reverse, reverse_names = mutation_context_features(mutant, wild, [("A", 1), ("A", 2)])
    np.testing.assert_array_equal(forward, -reverse)
    np.testing.assert_array_equal(names, reverse_names)


def test_missing_joint_mutant_reports_failure_instead_of_wt_fallback(tmp_path):
    folder, frame = sequence_fixture(tmp_path)
    store = SequenceFeatureStore(folder)
    del store.mutants[("double", "A")]
    with pytest.raises(KeyError):
        store.features(frame.iloc[2].to_dict(), task="mutation")
    mutants = pd.read_csv(folder / "mutant_chains.csv")
    mutants[mutants.sample_id.ne("double")].to_csv(folder / "mutant_chains.csv", index=False)
    report = prepare_sequence_features(frame, folder, tmp_path / "output", task="mutation")
    assert report["n_success"] == 3 and report["n_failed"] == 1
    failures = pd.read_csv(tmp_path / "output/sequence_coverage.csv")
    assert failures.loc[failures.sample_id.eq("double"), "status"].item() == "error"


def test_wrong_joint_mutant_identity_fails_even_when_embedding_shapes_match(tmp_path):
    folder, frame = sequence_fixture(tmp_path)
    store = SequenceFeatureStore(folder)
    store.mutants[("double", "A")] = "seq3"
    with pytest.raises(ValueError, match="joint mutant sequence"):
        store.features(frame.iloc[2].to_dict(), task="mutation")


@pytest.mark.parametrize("family", ["ridge", "rbf"])
@pytest.mark.parametrize("target", ["pkd", "ddg"])
def test_regression_fits_only_supplied_training_and_mutation_predictions_are_odd(family, target):
    rng = np.random.default_rng(12)
    matrix, labels = rng.normal(size=(16, 7)), rng.normal(size=16)
    groups = np.array([str(index // 3) for index in range(16)])
    model = fit_sequence_regressor(matrix[:12], labels[:12], groups[:12], target=target, family=family)
    predicted = model.predict(matrix[12:])
    assert np.isfinite(predicted).all()
    if target == "ddg":
        np.testing.assert_allclose(model.predict(-matrix[12:]), -predicted, atol=1e-12)
        np.testing.assert_allclose(model.predict(np.zeros((1, 7))), 0., atol=1e-12)
    else:
        np.testing.assert_allclose(model.scaler.mean_, matrix[:12].mean(axis=0), atol=1e-12)
