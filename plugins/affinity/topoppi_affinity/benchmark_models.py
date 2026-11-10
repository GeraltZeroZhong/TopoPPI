"""Train the frozen benchmark recipes on one supplied fold, without scoring it.

The caller owns splits, input preparation and evaluation. This module never
loads a development artifact or reads a validation affinity label.
"""
from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.kernel_ridge import KernelRidge
from sklearn.linear_model import Ridge
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from topoppi_affinity.bilinear_kernel import WeightedKernelRegressor
from topoppi_affinity.embedding_baselines import (
    candidate_parameters,
    fit_frozen_regressor,
    fit_kernel_block,
    zero_surface_slots,
)
from topoppi_affinity.evaluation import regression_metrics
from topoppi_affinity.features import feature_group
from topoppi_affinity.legacy_kernel import concatenate_blocks, fit_legacy_bundle
from topoppi_affinity.models import complex_weights
from topoppi_affinity.residual_baselines import ResidualRegressor, fit_origin_slope

ABSOLUTE_METHODS = (
    "topoppi_surface", "topoppi_contact", "bilinear_surface", "bilinear_contact",
    "legacy_original_recipe", "legacy_alternative_recipe", "null_train_mean",
)
VARIANT_METHODS = (
    "legacy_surface", "legacy_contact", "stab_raw", "ba_raw",
    "stab_surface_residual", "stab_contact_residual",
    "known_surface_ensemble", "known_contact_ensemble", "null_zero", "null_train_mean",
)
GENERALIZATION_METHODS = (
    "ridge_surface", "ridge_contact", "ba_raw",
    "generalization_surface_ensemble", "generalization_contact_ensemble",
    "legacy_surface", "legacy_contact", "null_zero", "null_train_mean",
)
TASK_METHODS = {
    "absolute_interaction": ABSOLUTE_METHODS, "absolute_family": ABSOLUTE_METHODS,
    "mutation_variant": VARIANT_METHODS,
    "mutation_interaction": GENERALIZATION_METHODS, "mutation_family": GENERALIZATION_METHODS,
}


@dataclass
class LoadedFeatures:
    sample_ids: np.ndarray
    X: np.ndarray
    feature_names: np.ndarray
    feature_groups: np.ndarray | None = None

    @classmethod
    def read(cls, value):
        if isinstance(value, cls):
            return value
        if isinstance(value, dict):
            return cls(**value)
        with np.load(value, allow_pickle=False) as source:
            return cls(source["sample_ids"].astype(str), source["X"], source["feature_names"].astype(str),
                       source["feature_groups"].astype(str) if "feature_groups" in source else None)

    def aligned(self, ids):
        lookup = {str(sample): index for index, sample in enumerate(self.sample_ids)}
        if len(lookup) != len(self.sample_ids):
            raise ValueError("Feature sample IDs must be unique")
        matrix = np.full((len(ids), len(self.feature_names)), np.nan, dtype=np.result_type(self.X.dtype, np.float32))
        for row, sample in enumerate(ids):
            if str(sample) in lookup:
                matrix[row] = self.X[lookup[str(sample)]]
        return matrix


def _odd_scaler(matrix, weights):
    """Exactly the scaler of weighted [X,-X], without materializing both."""
    scaler = StandardScaler().partial_fit(matrix, sample_weight=weights / 2.)
    scaler.partial_fit(-matrix, sample_weight=weights / 2.)
    return scaler


@dataclass
class OddRegressor:
    """Equivalent to reverse augmentation with half weights and odd prediction."""
    scaler: object
    estimator: object
    training: np.ndarray | None = None
    gamma: float | None = None

    def predict(self, matrix):
        values = self.scaler.transform(np.asarray(matrix, dtype=float))
        if self.training is None:
            return self.estimator.predict(values)
        gram = rbf_kernel(values, self.training, gamma=self.gamma)
        gram -= rbf_kernel(-values, self.training, gamma=self.gamma)
        gram *= .5
        return self.estimator.predict(gram)


def fit_odd_regressor(matrix, labels, groups, *, family, alpha, gamma=None):
    """Reduce the existing +/- augmented Ridge/RBF problem to N observations.

For RBF, K_odd=(K(X,X')-K(X,-X'))/2 with original weights and alpha is
algebraically identical to the 2N problem with half weights. No new tuning.
"""
    matrix = np.asarray(matrix, dtype=float)
    weights = complex_weights(groups)
    scaler = _odd_scaler(matrix, weights)
    values = scaler.transform(matrix)
    if family == "ridge":
        estimator = Ridge(alpha=alpha, fit_intercept=False).fit(values, labels, sample_weight=weights)
        return OddRegressor(scaler, estimator)
    if family != "rbf":
        raise ValueError("Odd reduction supports ridge and rbf")
    gram = rbf_kernel(values, gamma=gamma)
    gram -= rbf_kernel(-values, values, gamma=gamma)
    gram *= .5
    estimator = KernelRidge(alpha=alpha, kernel="precomputed").fit(gram, labels, sample_weight=weights)
    return OddRegressor(scaler, estimator, values, gamma)


def fit_benchmark_residual(matrix, prior, frame, *, seed=11):
    """The existing nine RBF candidates, two canonical-variant inner folds."""
    if frame.variant_group.nunique() < 2:
        return None, {}, [], "Fewer than two successful prior variant groups"
    splits = list(GroupKFold(n_splits=2, shuffle=True, random_state=seed).split(
        frame, groups=frame.variant_group))
    if any(not np.any(prior[train] != 0.) for train, _ in splits):
        return None, {}, [], "An inner training fold has no nonzero successful prior"
    labels, groups = frame.ddg.to_numpy(float), frame.complex_id.to_numpy()
    candidates, best_params, best_score = [], None, float("inf")
    for index, params in enumerate(candidate_parameters("rbf", matrix.shape[1])):
        predicted = np.full(len(frame), np.nan)
        for train, valid in splits:
            slope = fit_origin_slope(prior[train], labels[train], groups[train])
            model = fit_odd_regressor(matrix[train], labels[train] - slope * prior[train], groups[train],
                                      family="rbf", alpha=params["alpha"], gamma=params["gamma"])
            predicted[valid] = slope * prior[valid] + model.predict(matrix[valid])
            del model
        metrics = regression_metrics(labels, predicted, groups)
        candidates.append({"candidate": index, **params, **metrics})
        if metrics["macro_mae"] < best_score:
            best_params, best_score = params, metrics["macro_mae"]
    slope = fit_origin_slope(prior, labels, groups)
    model = fit_odd_regressor(matrix, labels - slope * prior, groups, family="rbf",
                              alpha=best_params["alpha"], gamma=best_params["gamma"])
    return ResidualRegressor(slope, model), best_params, candidates, ""


def benchmark_memory_estimate(n_absolute, n_mutation, mutation_dimension):
    """Analytical dense-array sizes for explicitly supplied counts and dimensions."""
    return {
        "n_absolute": n_absolute,
        "n_mutation": n_mutation,
        "mutation_dimension": mutation_dimension,
        "absolute_one_float64_gram_MiB": n_absolute**2 * 8 / 2**20,
        "mutation_one_float64_gram_MiB": n_mutation**2 * 8 / 2**20,
        "mutation_old_augmented_gram_MiB": (2*n_mutation)**2 * 8 / 2**20,
        "mutation_float64_features_MiB": n_mutation * mutation_dimension * 8 / 2**20,
        "scope": "Declared training-row counts and largest loaded mutation feature block; complete-case subsets may be smaller. Array sizes exclude solver workspaces and process overhead.",
        "odd_kernel": "N by N equivalent odd kernel, not a materialized 2N by 2N augmented kernel",
    }


def _prior_frame(value, frame, method):
    if value is None:
        result = pd.DataFrame({"sample_id": frame.sample_id.to_numpy(), "predicted": np.nan,
                               "status": "missing", "error": "Prior input not supplied"})
    else:
        source = value.copy() if isinstance(value, pd.DataFrame) else pd.read_csv(value, keep_default_na=False)
        if source.sample_id.duplicated().any():
            raise ValueError("Prior sample IDs must be unique")
        if not source.target.eq("ddg").all() or not source.method.eq(method).all():
            raise ValueError(f"Prior must use ddg and method {method}")
        result = source.set_index("sample_id").reindex(frame.sample_id).reset_index()
        result["status"] = result.status.fillna("missing")
        result["error"] = result.error.fillna("Missing prior sample")
    result["predicted"] = pd.to_numeric(result.predicted, errors="coerce")
    result["usable"] = result.status.eq("ok") & np.isfinite(result.predicted)
    result.loc[~result.usable & result.error.eq(""), "error"] = "Failed or nonfinite prior"
    return result


def fit_predict_benchmark_fold(train_frame, valid_frame, features, priors, *, task, output, seed=11):
    """Fit fixed recipes on train only; return predictions and training reports.

    features maps plm_surface/plm_analytic/legacy_v1/bilinear_pair_surface to
    NPZ paths or LoadedFeatures. priors maps stab/ba to unlabeled CSVs or frames.
    Missing features/priors produce explicit per-method/row failures; programming
    and schema errors propagate. No validation label is accessed or copied.
    """
    methods = TASK_METHODS[task]
    target = "pkd" if task.startswith("absolute_") else "ddg"
    train = train_frame.copy().reset_index(drop=True)
    valid_columns = [name for name in ("sample_id", "complex_id", "family_group", "interaction_group", "initial_interaction_group",
                                      "variant_group", "n_mutations", "complex_seen_in_training") if name in valid_frame]
    valid = valid_frame[valid_columns].copy().reset_index(drop=True)
    if train.sample_id.duplicated().any() or valid.sample_id.duplicated().any():
        raise ValueError("Fold sample IDs must be unique")
    if set(train.sample_id) & set(valid.sample_id):
        raise ValueError("Training and validation sample IDs overlap")
    if not np.isfinite(train[target].to_numpy(float)).all():
        raise ValueError("Training targets must be finite")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    priors = priors or {}
    started = time.monotonic()
    loaded, aligned = {}, {}
    predictions, reports = {}, {}

    def feature(key):
        if key not in features:
            return None
        if key not in loaded:
            block = LoadedFeatures.read(features[key])
            loaded[key] = block
            aligned[key] = (block.aligned(train.sample_id), block.aligned(valid.sample_id))
        return loaded[key], *aligned[key]

    def emit(name, values=None, *, usable=None, error="", report=None, model=None, **extra):
        result = valid.copy()
        result["method"], result["target"] = name, target
        result["predicted"] = np.nan if values is None else values
        usable = np.isfinite(result.predicted.to_numpy()) if usable is None else np.asarray(usable, bool)
        result["status"] = np.where(usable, "ok", "error")
        result["error"] = np.where(usable, "", error)
        for key, value in extra.items():
            result[key] = value
        metadata = {"task": task, "target": target, "method": name, "seed": seed,
                    "n_train_declared": len(train), "n_validation_declared": len(valid),
                    "n_predictions": int(usable.sum()), "n_train": 0, "train_sample_ids": [],
                    "validation_labels_read": False,
                    "fit_status": "not_fitted" if values is None else "ok", **(report or {})}
        destination = output / name
        destination.mkdir(exist_ok=True)
        if model is not None:
            joblib.dump({"model": model, "report": metadata}, destination / "model.joblib")
        result.to_csv(destination / "validation_predictions.csv", index=False)
        (destination / "report.json").write_text(json.dumps(metadata, indent=2) + "\n")
        predictions[name], reports[name] = result, metadata

    def learned(name, keys, fit, *, surface=True, prior=None):
        blocks = [feature(key) for key in keys]
        if any(block is None for block in blocks):
            emit(name, error="Required feature block not supplied", status="missing_features", report={"feature_blocks": keys})
            return
        train_ok = np.logical_and.reduce([np.isfinite(x).all(axis=1) for _, x, _ in blocks])
        valid_features_ok = np.logical_and.reduce([np.isfinite(xv).all(axis=1) for _, _, xv in blocks])
        valid_ok = valid_features_ok.copy()
        arrays = []
        for key, (block, x, xv) in zip(keys, blocks, strict=True):
            if not surface:
                if key == "legacy_v1":
                    groups = block.feature_groups if block.feature_groups is not None else np.asarray(
                        [feature_group(str(n)) for n in block.feature_names])
                    mask = np.asarray(groups) == "surface"
                    x, xv = x.copy(), xv.copy()
                    x[:, mask], xv[:, mask] = 0., 0.
                else:
                    x, xv = zero_surface_slots(x, block.feature_names), zero_surface_slots(xv, block.feature_names)
            arrays.append((block, x, xv))
        if prior is not None:
            train_ok &= prior[0].usable.to_numpy()
            valid_ok &= prior[1].usable.to_numpy()
        selected = train.loc[train_ok].reset_index(drop=True)
        metadata = {"feature_blocks": keys, "use_surface": surface, "n_train": len(selected),
                    "train_sample_ids": selected.sample_id.tolist(),
                    "excluded_train_sample_ids": train.loc[~train_ok, "sample_id"].tolist(),
                    "n_validation_features_usable": int(valid_features_ok.sum()),
                    "feature_names": {key: np.asarray(block[0].feature_names).astype(str).tolist()
                                      for key, block in zip(keys, blocks, strict=True)}}
        if "plm_analytic" in keys:
            metadata["curvature"] = "analytic_Gaussian_Goldman"
        if len(selected) < 2:
            emit(name, error="Fewer than two usable training rows", status="not_fitted", report=metadata)
            return
        inputs = [(block, x[train_ok], xv[valid_ok]) for block, x, xv in arrays]
        model, params, candidates, reason = fit(selected, inputs, train_ok)
        if model is None:
            emit(name, error=reason, status="not_fitted", report={**metadata, "reason": reason})
            return
        values = np.full(len(valid), np.nan)
        if valid_ok.any():
            if prior is not None:
                values[valid_ok] = model.predict(prior[1].predicted.to_numpy()[valid_ok], inputs[0][2])
            elif len(inputs) == 1:
                values[valid_ok] = model.predict(inputs[0][2])
            else:
                values[valid_ok] = model.predict(np.concatenate([item[2] for item in inputs], axis=1))
        statuses = np.where(valid_ok, "ok", np.where(valid_features_ok, "prior_failed", "missing_features"))
        emit(name, values, usable=valid_ok, error="Missing/nonfinite features or unsuccessful prior", status=statuses,
             report={**metadata, "parameters": params}, model=model)
        if candidates:
            pd.DataFrame(candidates).to_csv(output / name / "inner_candidates.csv", index=False)
        del model

    def frozen(family, params, *, odd=False):
        def fit(frame, inputs, _):
            x = inputs[0][1]
            if odd:
                model = fit_odd_regressor(x, frame[target].to_numpy(float), frame.complex_id,
                                          family=family, alpha=params["alpha"], gamma=params.get("gamma"))
            else:
                model = fit_frozen_regressor(x, frame[target], frame.complex_id, family=family,
                                             params=params, antisymmetric=False, seed=seed)
            return model, params, [], ""
        return fit

    def raw(name, prior):
        emit(name, np.where(prior.usable, prior.predicted, np.nan), usable=prior.usable,
             error=prior.error.to_numpy(), status=np.where(prior.usable, "ok", "prior_failed"),
             report={"fit_status": "frozen_prior", "n_train": 0}, model={"kind": "frozen_scalar_prior", "name": name})

    def ensemble(name, left, right):
        base, other = predictions[left], predictions[right]
        usable = base.status.eq("ok").to_numpy()
        both = usable & other.status.eq("ok").to_numpy()
        values = base.predicted.to_numpy().copy()
        values[both] = .5 * (values[both] + other.predicted.to_numpy()[both])
        components = np.where(both, "equal_ensemble", np.where(usable, "base_only", "unavailable"))
        emit(name, values, usable=usable, error=base.error.to_numpy(), components_used=components,
             status=np.where(usable, "ok", base.status),
             base_prediction=base.predicted.to_numpy(), complementary_prediction=other.predicted.to_numpy(),
             report={"combination": ".5 * base + .5 * complement; base_only on complement failure",
                     "source_methods": [left, right], "n_base_only": int((usable & ~both).sum()),
                     "n_train": reports[left].get("n_train", 0),
                     "n_complement_train": reports[right].get("n_train", 0),
                     "train_sample_ids": reports[left].get("train_sample_ids", []),
                     "complement_train_sample_ids": reports[right].get("train_sample_ids", [])},
             model={"component_models": [str(output / part / "model.joblib")
                     if (output / part / "model.joblib").exists() else None for part in (left, right)],
                    "weights": [.5, .5]})

    with threadpool_limits(limits=2):
        if target == "pkd":
            for surface in (True, False):
                label = "surface" if surface else "contact"

                def absolute_fit(frame, inputs, _, bilinear=False):
                    (p, x, _), (legacy_block, legacy, _) = inputs
                    groups = legacy_block.feature_groups if legacy_block.feature_groups is not None else np.asarray(
                        [feature_group(str(n)) for n in legacy_block.feature_names])
                    if not bilinear:
                        bundle = fit_legacy_bundle(x, p.feature_names, legacy, legacy_block.feature_names, groups, frame,
                                                   alpha=1., gamma_factor=.1, legacy_weight=.25, use_surface=True)
                        return bundle["model"], {"alpha": 1., "gamma_factor": .1, "weights": [.75, .25]}, [], ""
                    matrix, _, legacy_mask = concatenate_blocks(x, p.feature_names, legacy, legacy_block.feature_names,
                                                                 groups, use_surface=True)
                    pnames = np.asarray(p.feature_names).astype(str)
                    pair = np.concatenate([np.char.startswith(pnames, "pair.") | np.char.startswith(pnames, "surface.pair."),
                                           np.zeros(len(legacy_block.feature_names), dtype=bool)])
                    weights = complex_weights(frame.complex_id)
                    center = float(np.average(frame.pkd, weights=weights))
                    mixture, blocks = (.375, .125, .5), []
                    gram = np.zeros((len(frame), len(frame)))
                    for mask, fraction in zip((~legacy_mask & ~pair, legacy_mask, pair), mixture, strict=True):
                        block, component = fit_kernel_block(matrix, mask, weights, .1)
                        blocks.append(block)
                        gram += fraction * component
                        del component
                    estimator = KernelRidge(alpha=1., kernel="precomputed").fit(
                        gram, frame.pkd.to_numpy()-center, sample_weight=weights)
                    return WeightedKernelRegressor(estimator, tuple(blocks), mixture, center), {
                        "alpha": 1., "gamma_factor": .1, "weights": list(mixture)}, [], ""

                learned(f"topoppi_{label}", ["plm_surface", "legacy_v1"], absolute_fit, surface=surface)
                learned(f"bilinear_{label}", ["bilinear_pair_surface", "legacy_v1"],
                        lambda f, x, mask: absolute_fit(f, x, mask, bilinear=True), surface=surface)
            legacy = feature("legacy_v1")
            for name, alpha, factor in (("legacy_original_recipe", 10., 1.), ("legacy_alternative_recipe", 1., .1)):
                if legacy is None:
                    emit(name, error="Required legacy feature block not supplied")
                else:
                    learned(name, ["legacy_v1"], frozen("rbf", {
                        "alpha": alpha, "gamma": factor / len(legacy[0].feature_names), "gamma_factor": factor}))
        else:
            ba_valid = _prior_frame(priors.get("ba"), valid, "BA-Cycle-partner-MC8")
            raw("ba_raw", ba_valid)
            leaf = 2 if task == "mutation_variant" else 5
            for surface in (True, False):
                label = "surface" if surface else "contact"
                learned(f"legacy_{label}", ["legacy_v1"], frozen("trees", {"min_samples_leaf": leaf}), surface=surface)
            if task == "mutation_variant":
                stab_train = _prior_frame(priors.get("stab"), train, "StaB-ddG-stability-only")
                stab_valid = _prior_frame(priors.get("stab"), valid, "StaB-ddG-stability-only")
                raw("stab_raw", stab_valid)
                for surface in (True, False):
                    label = "surface" if surface else "contact"

                    def residual_fit(frame, inputs, mask):
                        model, params, candidates, reason = fit_benchmark_residual(
                            inputs[0][1], stab_train.predicted.to_numpy()[mask], frame, seed=seed)
                        if model is not None:
                            params = {**params, "origin_slope": model.slope, "inner_group": "variant_group", "inner_folds": 2}
                        return model, params, candidates, reason

                    learned(f"stab_{label}_residual", ["plm_surface"], residual_fit, surface=surface,
                            prior=(stab_train, stab_valid))
                    ensemble(f"known_{label}_ensemble", f"legacy_{label}", f"stab_{label}_residual")
            else:
                for surface in (True, False):
                    label = "surface" if surface else "contact"
                    learned(f"ridge_{label}", ["plm_analytic"], frozen("ridge", {"alpha": 10000.}, odd=True), surface=surface)
                    ensemble(f"generalization_{label}_ensemble", f"ridge_{label}", "ba_raw")
            emit("null_zero", np.zeros(len(valid)), report={"constant": 0., "n_train": 0})
        if len(train):
            mean = float(np.average(train[target], weights=complex_weights(train.complex_id)))
            emit("null_train_mean", np.full(len(valid), mean), report={"constant": mean, "n_train": len(train),
                 "weighting": "equal total weight per complex", "train_sample_ids": train.sample_id.tolist()})
        else:
            emit("null_train_mean", error="Empty training fold")
    result = pd.concat([predictions[name] for name in methods], ignore_index=True)
    mutation_dimension = max((len(block.feature_names) for block in loaded.values()), default=0) if target == "ddg" else 0
    report = {"task": task, "target": target, "methods": list(methods), "train_sample_ids": train.sample_id.tolist(),
              "validation_sample_ids": valid.sample_id.tolist(), "validation_labels_read": False,
              "preprocessing": "All fitted scalers, kernel centers, target means and prior calibration use training only",
              "memory_estimate": benchmark_memory_estimate(
                  len(train) if target == "pkd" else 0, len(train) if target == "ddg" else 0, mutation_dimension),
              "elapsed_seconds": time.monotonic()-started,
              "method_reports": reports}
    result.to_csv(output / "predictions.csv", index=False)
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return {"predictions": result, "report": report}
