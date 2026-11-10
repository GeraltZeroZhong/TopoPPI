"""The distributed default models run from user structures without downloads."""
import json
from pathlib import Path

import numpy as np
import pytest

from topoppi_affinity.production import PhysicalRegressor, bundled_model
from topoppi_affinity.resources import resource_status
from topoppi_affinity.workflow import analyze, load_analysis


@pytest.mark.parametrize("task", ["absolute", "mutation"])
def test_bundled_models_are_ready_without_a_cache(task, tmp_path, monkeypatch):
    import urllib.request
    monkeypatch.setattr(urllib.request, "urlopen", lambda *a, **kw: pytest.fail("Default model attempted a download"))
    state = resource_status(cache_dir=tmp_path)["models"][task]
    assert state["installed"] and state["ready"] and state["bundled"]
    model = PhysicalRegressor(bundled_model(task))
    X = model.arrays["median"][None, :].copy()
    expected = model.predict(X, model.feature_names)
    assert np.isfinite(expected).all()
    np.testing.assert_allclose(model.predict(X[:, ::-1], model.feature_names[::-1]), expected)


@pytest.mark.parametrize("task", ["absolute", "mutation"])
def test_default_structure_prediction_and_reopening(task, tmp_path, monkeypatch):
    import urllib.request
    from importlib.resources import files
    monkeypatch.setattr(urllib.request, "urlopen", lambda *a, **kw: pytest.fail("Default workflow attempted a download"))
    structure = Path(str(files("topoppi.assets").joinpath("examples", "synthetic_complex.pdb")))
    options = {"mutations": ["A:A1V", "A:A1V,B:F2Y"]} if task == "mutation" else {}
    output = tmp_path / "result"
    result = analyze(structure, "A", "B", output=output, task=task,
                     cache_dir=tmp_path / "empty-cache", log=lambda text: None, **options)
    assert result["prediction_status"] == "ok", result.get("prediction_error")
    assert result["device"] == "cpu"
    saved = load_analysis(output)
    assert np.isfinite(saved.predictions.predicted).all()
    assert saved.interface is not None
    assert saved.report["model"]["feature_scheme"]
    assert (output / "predictions.csv").is_file()
    assert (output / "residues.csv").is_file()


def test_wrong_target_is_reported_before_computation(tmp_path):
    from importlib.resources import files

    from topoppi_affinity.workflow import validate_inputs
    structure = Path(str(files("topoppi.assets").joinpath("examples", "synthetic_complex.pdb")))
    with pytest.raises(ValueError, match="conflicts"):
        validate_inputs(structure, "A", "B", task="mutation", model=bundled_model("absolute"),
                        mutations=["A:A1V"], output=tmp_path / "result")


def test_portable_kernel_matches_weighted_sklearn_with_missing_inputs(tmp_path):
    from sklearn.impute import SimpleImputer
    from sklearn.kernel_ridge import KernelRidge
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    rng = np.random.default_rng(19)
    X = rng.normal(size=(31, 5))
    X[::4, 1] = np.nan
    X[:, 4] = np.nan
    target = rng.normal(size=len(X))
    weights = np.linspace(.2, 2., len(X))
    alpha, gamma = 1., 1. / X.shape[1]
    model = make_pipeline(SimpleImputer(strategy="median", add_indicator=True, keep_empty_features=True),
                          StandardScaler(), KernelRidge(kernel="rbf", alpha=alpha, gamma=gamma))
    model.fit(X, target, kernelridge__sample_weight=weights)
    imputer, scaler, kernel = [step for _, step in model.steps]
    np.savez_compressed(tmp_path / "weights.npz", median=imputer.statistics_,
                        indicator_indices=imputer.indicator_.features_, mean=scaler.mean_, scale=scaler.scale_,
                        training_vectors=kernel.X_fit_, dual_coefficient=kernel.dual_coef_)
    names = [f"feature_{i}" for i in range(X.shape[1])]
    path = tmp_path / "model.json"
    path.write_text(json.dumps({"target": "ddg", "feature_names": names, "weights": "weights.npz",
                               "training": {"family": "rbf", "parameters": {"alpha": alpha, "gamma": gamma}}}))
    portable = PhysicalRegressor(path)
    # Include fitted vectors and unseen queries spanning more than two prediction blocks.
    queries = np.vstack([X, rng.normal(size=(513, X.shape[1]))])
    queries[::3, 1] = np.nan
    queries[:, 4] = np.nan
    np.testing.assert_allclose(portable.predict(queries[:, ::-1], names[::-1]),
                               model.predict(queries), rtol=1e-11, atol=1e-11)


def test_default_mutation_is_kernel_ridge_and_preserves_zero_substitution():
    model = PhysicalRegressor(bundled_model("mutation"))
    assert model.family == "rbf"
    assert "coefficient" not in model.arrays
    assert len(model.feature_names) == 431
    values = model.arrays["median"][None, :].copy()
    values[0, model.feature_names.index("mutation.n_mutations")] = 0.
    np.testing.assert_array_equal(model.predict(values, model.feature_names), [0.])
