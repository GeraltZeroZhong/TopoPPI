"""The portable mutation readout preserves fitted preprocessing and RBF scores."""
import json

import numpy as np
import pytest
from sklearn.impute import SimpleImputer
from sklearn.kernel_ridge import KernelRidge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from topoppi_hotspot.mutation_model import MutationSurfaceRegressor


def test_portable_regression_matches_fitted_pipeline_with_missing_values(tmp_path):
    rng = np.random.default_rng(18)
    X = rng.normal(size=(24, 5))
    X[::3, 1] = np.nan
    X[:, 4] = np.nan
    target = X[:, 0] - X[:, 2] * .4
    pipeline = make_pipeline(
        SimpleImputer(strategy="median", add_indicator=True, keep_empty_features=True),
        StandardScaler(), KernelRidge(kernel="rbf", alpha=1., gamma=.2),
    )
    pipeline.fit(X, target, kernelridge__sample_weight=np.linspace(.5, 1.5, len(X)))
    imputer, scaler, regressor = (step for _, step in pipeline.steps)
    np.savez_compressed(
        tmp_path / "weights.npz", median=imputer.statistics_,
        indicator_indices=imputer.indicator_.features_, mean=scaler.mean_, scale=scaler.scale_,
        training_vectors=regressor.X_fit_, dual_coefficient=regressor.dual_coef_,
    )
    names = ["a", "b", "c", "d", "e"]
    (tmp_path / "model.json").write_text(json.dumps({
        "feature_names": names, "target": "ddg", "weights": "weights.npz",
        "training": {"parameters": {"gamma": .2}},
    }))
    model = MutationSurfaceRegressor(tmp_path / "model.json")
    queries = rng.normal(size=(301, 5))
    queries[::4, 1] = np.nan
    queries[:, 4] = np.nan
    np.testing.assert_allclose(model.predict(queries, names), pipeline.predict(queries), atol=1e-12)
    np.testing.assert_allclose(model.predict(queries[:, ::-1], names[::-1]), pipeline.predict(queries), atol=1e-12)
    assert model.predict(np.empty((0, 5)), names).shape == (0,)
    with pytest.raises(ValueError, match="Missing model features"):
        model.predict(queries[:, :-1], names[:-1])


def test_bundled_mutation_model_has_physical_and_surface_inputs():
    from importlib.resources import files

    path = files("topoppi_hotspot").joinpath("data", "mutation_surface", "model.json")
    model = MutationSurfaceRegressor(path)
    assert len(model.feature_names) == 552
    assert model.metadata["kind"] == "mutation_surface"
    assert model.metadata["encoding"]["budget_ratio"] == .25
    assert model.arrays["training_vectors"].shape[1] == len(model.arrays["mean"])
    assert np.isfinite(model.predict(model.arrays["median"][None, :], model.feature_names)).all()
