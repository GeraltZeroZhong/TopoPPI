"""Regression models for affinity descriptors and grouped evaluation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import joblib
import numpy as np
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from topoppi_affinity.thermodynamics import pkd_to_dg, pkd_to_kd


def complex_weights(groups):
    _, inverse, counts = np.unique(np.asarray(groups).astype(str), return_inverse=True, return_counts=True)
    weights = 1.0 / counts[inverse]
    return weights / np.mean(weights)


@dataclass
class AffinityRegressor:
    estimator: object
    feature_names: tuple[str, ...]
    target: str
    model_name: str
    training_metadata: dict

    def predict(self, X, feature_names):
        names = list(map(str, feature_names))
        if len(names) != len(set(names)):
            raise ValueError("Duplicate feature names")
        lookup = {name: index for index, name in enumerate(names)}
        missing = set(self.feature_names) - lookup.keys()
        if missing:
            raise ValueError(f"Missing model features: {sorted(missing)[:5]}")
        matrix = np.asarray(X, dtype=float)
        if matrix.ndim != 2 or matrix.shape[1] != len(names):
            raise ValueError("Feature matrix does not match feature names")
        selected = matrix[:, [lookup[name] for name in self.feature_names]]
        if np.isinf(selected).any():
            raise ValueError("Infinite feature values")
        predictions = np.asarray(self.estimator.predict(selected), dtype=float)
        # WT complex context remains nonzero for a validated identity mutation.
        if self.target == "ddg":
            no_change = (matrix[:, lookup["mutation.n_mutations"]] == 0
                         if "mutation.n_mutations" in lookup else np.all(matrix == 0, axis=1))
            predictions[no_change] = 0.0
        return predictions

    def predict_affinity(self, X, feature_names, *, temperature_K=298.15):
        if self.target != "pkd":
            raise ValueError("Absolute affinity conversion requires a pKd model")
        pkd = self.predict(X, feature_names)
        return {"pkd": pkd, "kd_M": pkd_to_kd(pkd),
                "dg_reference_kcal_mol": pkd_to_dg(pkd, temperature_K)}

    def save(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        joblib.dump(self, path)


def fit_regressor(X, y, feature_names, *, target, model_name, groups, seed=42, n_jobs=4):
    X, y = np.asarray(X, dtype=float), np.asarray(y, dtype=float)
    if X.ndim != 2 or X.shape != (len(y), len(feature_names)) or len(y) < 2:
        raise ValueError("Training rows and feature schema must align")
    if not np.isfinite(y).all() or np.isinf(X).any():
        raise ValueError("Nonfinite training targets or infinite features")
    if target not in {"pkd", "ddg"}:
        raise ValueError("Target must be pkd or ddg")
    if model_name.startswith("ridge_"):
        alpha = float(model_name.split("_", 1)[1])
        regressor = Ridge(alpha=alpha)
        steps = [("impute", SimpleImputer(strategy="median", keep_empty_features=True)),
                 ("scale", StandardScaler()), ("regressor", regressor)]
    elif model_name.startswith("trees_"):
        leaf = int(model_name.split("_", 1)[1])
        regressor = ExtraTreesRegressor(n_estimators=192, min_samples_leaf=leaf,
                                       max_features=0.75, random_state=seed, n_jobs=n_jobs)
        steps = [("impute", SimpleImputer(strategy="median", keep_empty_features=True)),
                 ("regressor", regressor)]
    else:
        raise ValueError(f"Unknown model: {model_name}")
    pipeline = Pipeline(steps)
    pipeline.fit(X, y, regressor__sample_weight=complex_weights(groups))
    return AffinityRegressor(pipeline, tuple(map(str, feature_names)), target, model_name,
                             {"seed": seed, "n_train": len(y), "n_train_complexes": len(set(groups)),
                              "weighting": "equal total weight per training complex",
                              "ddg_sign": "positive means weaker mutant binding",
                              "geometry": "wild-type native Gaussian interface; explicit mutation identity"})


def feature_sets(feature_groups, *, target):
    groups = np.asarray(feature_groups).astype(str)
    base = {"sequence", "contact"}
    if target == "ddg":
        base |= {"mutation", "epistasis"}
    masks = {"contact_baseline": np.isin(groups, sorted(base)),
             "topoppi": np.ones(len(groups), dtype=bool)}
    if target == "pkd":
        masks["sequence_composition"] = groups == "sequence"
    return {key: value for key, value in masks.items() if np.any(value)}

