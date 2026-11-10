"""Fold isolation, actual old/new solver equivalence and coverage behavior."""
import joblib
import numpy as np
import pandas as pd
import pytest

from topoppi_affinity import benchmark_models as benchmark
from topoppi_affinity.embedding_baselines import fit_frozen_regressor
from topoppi_affinity.residual_baselines import fit_origin_slope


def fixture_data(n=24):
    rng = np.random.default_rng(19)
    frame = pd.DataFrame({
        "sample_id": [f"sample_{i}" for i in range(n)],
        "complex_id": [f"complex_{i // 3}" for i in range(n)],
        "family_group": [f"family_{i // 6}" for i in range(n)],
        "initial_interaction_group": [f"interaction_{i // 3}" for i in range(n)],
        "variant_group": [f"variant_{i // 2}" for i in range(n)],
        "mutations": [f"AA{i + 1}V" for i in range(n)],
        "n_mutations": np.where(np.arange(n) % 3 == 0, 2, 1),
        "pkd": rng.normal(6., 1., n), "ddg": rng.normal(size=n),
    })
    ids = frame.sample_id.to_numpy(dtype=str)
    plm = rng.normal(size=(n, 4)).astype(np.float32)
    legacy = rng.normal(size=(n, 3)).astype(np.float32)
    pair = np.column_stack([plm, rng.normal(size=(n, 2))]).astype(np.float32)
    features = {
        "plm_surface": benchmark.LoadedFeatures(ids[::-1], plm[::-1], np.array(["plm.a", "contact.b", "surface.c", "surface.d"])),
        "plm_analytic": benchmark.LoadedFeatures(ids[::-1], plm[::-1], np.array(["plm.a", "contact.b", "surface.c", "surface.d"])),
        "legacy_v1": benchmark.LoadedFeatures(ids, legacy, np.array(["sequence.a", "hidden_area", "contact.c"]),
                                              np.array(["sequence", "surface", "contact"])),
        "bilinear_pair_surface": benchmark.LoadedFeatures(ids, pair, np.array([
            "plm.a", "contact.b", "surface.c", "surface.d", "pair.contact.0", "surface.pair.0"])),
    }
    priors = {}
    for key, method in (("stab", "StaB-ddG-stability-only"), ("ba", "BA-Cycle-partner-MC8")):
        priors[key] = pd.DataFrame({"sample_id": ids[::-1], "predicted": rng.normal(size=n),
                                   "target": "ddg", "method": method, "status": "ok", "error": ""})
    return frame, features, priors


@pytest.mark.parametrize("family", ["ridge", "rbf"])
def test_odd_reduction_matches_old_augmentation_with_unequal_complex_weights(family):
    rng = np.random.default_rng(77)
    x = rng.normal(size=(13, 7)) + .7
    groups = np.array(["a"] * 8 + ["b"] * 3 + ["c"] * 2)
    labels, prior = rng.normal(size=13), rng.normal(size=13)
    slope = fit_origin_slope(prior, labels, groups)
    residual = labels - slope * prior
    params = {"alpha": .1, "gamma": .1 / x.shape[1]}
    old = fit_frozen_regressor(x, residual, groups, family=family, params=params, antisymmetric=True)
    new = benchmark.fit_odd_regressor(x, residual, groups, family=family, **params)
    query = rng.normal(size=(5, 7))
    np.testing.assert_allclose(new.predict(query), old.predict(query), rtol=1e-10, atol=1e-11)
    np.testing.assert_allclose(new.predict(-query), -new.predict(query), atol=1e-12)
    np.testing.assert_allclose(new.predict(np.zeros((1, 7))), 0., atol=1e-12)
    np.testing.assert_allclose(new.scaler.scale_, old.scaler.scale_, rtol=1e-13)
    if family == "rbf":
        assert len(new.training) == 13
        assert len(old.estimator.X_fit_) == 26


def test_absolute_fold_ignores_validation_labels_and_contact_uses_both_masks(tmp_path):
    frame, features, priors = fixture_data()
    train, valid = frame.iloc[:18], frame.iloc[18:]
    first = benchmark.fit_predict_benchmark_fold(train, valid, features, priors,
                                                task="absolute_family", output=tmp_path / "first")
    # Poison labels rather than merely permute them: prediction must not read them.
    changed = valid.assign(pkd="not an affinity label", ddg="not an affinity label")
    second = benchmark.fit_predict_benchmark_fold(train, changed, features, priors,
                                                 task="absolute_family", output=tmp_path / "second")
    np.testing.assert_array_equal(first["predictions"].predicted, second["predictions"].predicted)
    assert "pkd" not in first["predictions"] and "ddg" not in first["predictions"]
    memory = first["report"]["memory_estimate"]
    assert memory["n_absolute"] == 18 and memory["n_mutation"] == 0
    assert memory["absolute_one_float64_gram_MiB"] == 18**2 * 8 / 2**20
    model = joblib.load(tmp_path / "first/topoppi_surface/model.joblib")["model"]
    train_x = features["plm_surface"].aligned(train.sample_id)
    np.testing.assert_allclose(model.blocks[0].scaler.mean_, train_x.mean(axis=0), atol=1e-7)
    assert first["report"]["method_reports"]["legacy_original_recipe"]["parameters"] == {
        "alpha": 10., "gamma": 1./3., "gamma_factor": 1.}
    # Real surface values in both independent kernels must be absent from contact.
    for key in ("plm_surface", "bilinear_pair_surface"):
        block = features[key]
        mask = np.char.startswith(block.feature_names, "surface.")
        block.X[:, mask] = block.X[:, mask] * 1000 + 10000
    features["legacy_v1"].X[:, 1] += 10000
    third = benchmark.fit_predict_benchmark_fold(train, valid, features, priors,
                                                task="absolute_family", output=tmp_path / "third")
    for method in ("topoppi_contact", "bilinear_contact"):
        left = first["predictions"].loc[first["predictions"].method.eq(method), "predicted"]
        right = third["predictions"].loc[third["predictions"].method.eq(method), "predicted"]
        np.testing.assert_array_equal(left, right)


def test_family_coverage_matched_contacts_prior_fallback_and_npz_alignment(tmp_path):
    frame, features, priors = fixture_data()
    train, valid = frame.iloc[:18], frame.iloc[18:]
    block = features["plm_analytic"]
    # Only a surface coordinate is invalid: contact still uses the same cohort.
    for sample in [train.sample_id.iloc[0], valid.sample_id.iloc[0]]:
        block.X[np.flatnonzero(block.sample_ids == sample), 2] = np.nan
    path = tmp_path / "analytic.npz"
    np.savez(path, sample_ids=block.sample_ids, X=block.X, feature_names=block.feature_names)
    features["plm_analytic"] = path
    failed_prior = valid.sample_id.iloc[1]
    priors["ba"].loc[priors["ba"].sample_id.eq(failed_prior), ["predicted", "status", "error"]] = [np.nan, "error", "missing backbone"]
    result = benchmark.fit_predict_benchmark_fold(train, valid.drop(columns=["ddg", "pkd"]), features, priors,
                                                 task="mutation_family", output=tmp_path / "fold")
    reports = result["report"]["method_reports"]
    assert reports["ridge_surface"]["train_sample_ids"] == reports["ridge_contact"]["train_sample_ids"]
    assert reports["ridge_surface"]["n_train"] == 17
    predicted = result["predictions"].set_index(["method", "sample_id"])
    absent = valid.sample_id.iloc[0]
    for method in ("ridge_surface", "ridge_contact", "generalization_surface_ensemble"):
        assert predicted.loc[(method, absent), "status"] == "missing_features"
        assert np.isnan(predicted.loc[(method, absent), "predicted"])
    assert predicted.loc[("legacy_surface", absent), "status"] == "ok"
    assert predicted.loc[("ba_raw", failed_prior), "status"] == "prior_failed"
    combined = predicted.loc[("generalization_surface_ensemble", failed_prior)]
    assert combined.components_used == "base_only"
    assert combined.predicted == predicted.loc[("ridge_surface", failed_prior), "predicted"]
    assert reports["ridge_surface"]["curvature"] == "analytic_Gaussian_Goldman"
    memory = result["report"]["memory_estimate"]
    assert memory["n_absolute"] == 0 and memory["n_mutation"] == 18
    assert memory["mutation_dimension"] == 4
    assert memory["mutation_float64_features_MiB"] == 18 * 4 * 8 / 2**20


def test_known_nested_residual_only_training_and_failed_prior_keeps_base(tmp_path, monkeypatch):
    frame, features, priors = fixture_data()
    train, valid = frame.iloc[:18], frame.iloc[18:]
    monkeypatch.setattr(benchmark, "candidate_parameters", lambda *_: [
        {"alpha": .1, "gamma": .1/4, "gamma_factor": .1},
        {"alpha": 10., "gamma": 1./4, "gamma_factor": 1.},
    ])
    failed = valid.sample_id.iloc[-1]
    priors["stab"].loc[priors["stab"].sample_id.eq(failed), ["predicted", "status", "error"]] = [np.nan, "error", "missing backbone"]
    first = benchmark.fit_predict_benchmark_fold(train, valid, features, priors,
                                                task="mutation_variant", output=tmp_path / "first")
    second = benchmark.fit_predict_benchmark_fold(train, valid.assign(ddg=object()), features, priors,
                                                 task="mutation_variant", output=tmp_path / "second")
    for method in ("stab_surface_residual", "stab_contact_residual"):
        before = first["report"]["method_reports"][method]
        after = second["report"]["method_reports"][method]
        assert before["parameters"] == after["parameters"]
        assert set(before["train_sample_ids"]) == set(train.sample_id)
        assert before["parameters"]["inner_group"] == "variant_group"
    np.testing.assert_allclose(first["predictions"].predicted, second["predictions"].predicted, equal_nan=True)
    rows = first["predictions"].set_index(["method", "sample_id"])
    assert rows.loc[("known_surface_ensemble", failed), "components_used"] == "base_only"
    assert rows.loc[("known_surface_ensemble", failed), "predicted"] == rows.loc[("legacy_surface", failed), "predicted"]
    assert set(first["predictions"].method) == set(benchmark.VARIANT_METHODS)


def test_insufficient_priors_do_not_cancel_other_methods(tmp_path):
    frame, features, priors = fixture_data()
    priors["stab"].loc[:, "status"] = "error"
    result = benchmark.fit_predict_benchmark_fold(frame.iloc[:18], frame.iloc[18:], features, priors,
                                                 task="mutation_variant", output=tmp_path)
    rows = result["predictions"]
    assert rows.query("method == 'stab_surface_residual'").status.eq("not_fitted").all()
    assert rows.query("method == 'known_surface_ensemble'").components_used.eq("base_only").all()
    assert rows.query("method == 'legacy_surface'").status.eq("ok").all()
    estimate = benchmark.benchmark_memory_estimate(7, 11, 13)
    assert estimate["mutation_old_augmented_gram_MiB"] == 4 * estimate["mutation_one_float64_gram_MiB"]
    assert estimate["mutation_float64_features_MiB"] == 11 * 13 * 8 / 2**20


def test_no_training_validation_id_overlap(tmp_path):
    frame, features, priors = fixture_data()
    with pytest.raises(ValueError, match="overlap"):
        benchmark.fit_predict_benchmark_fold(frame.iloc[:18], frame.iloc[17:], features, priors,
                                             task="absolute_family", output=tmp_path)
