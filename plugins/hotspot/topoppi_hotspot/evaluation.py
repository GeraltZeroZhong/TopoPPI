"""Nested family-held-out evaluation; no pretrained Affinity models are used."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import average_precision_score, brier_score_loss, mean_absolute_error, roc_auc_score
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def metrics(y, energy, p, predicted_energy):
    y, energy, p, predicted_energy = map(np.asarray, (y, energy, p, predicted_energy))
    both = len(np.unique(y)) == 2
    count = min(3, len(y))
    cutoff = np.sort(p)[-count]
    above, tied = p > cutoff, p == cutoff
    precision_at_3 = (y[above].sum() + (count - above.sum()) * y[tied].mean()) / count
    return {
        "n": len(y), "hotspots": int(y.sum()), "prevalence": float(np.mean(y)),
        "AP": float(average_precision_score(y, p)) if y.sum() else np.nan,
        "AUROC": float(roc_auc_score(y, p)) if both else np.nan,
        "Brier": float(brier_score_loss(y, p)),
        "MAE": float(mean_absolute_error(energy, predicted_energy)),
        "RMSE": float(np.sqrt(np.mean((energy - predicted_energy) ** 2))),
        "Spearman": float(spearmanr(energy, predicted_energy).statistic)
        if len(y) >= 3 and np.ptp(energy) > 0 and np.ptp(predicted_energy) > 0 else np.nan,
        "precision_at_3": float(precision_at_3),
    }


def _fit(X, y, value, endpoint):
    estimator = (LogisticRegression(C=value, solver="lbfgs", max_iter=3000)
                 if endpoint == "primary" else Ridge(alpha=value))
    return make_pipeline(StandardScaler(), estimator).fit(X, y)


def _tune(X, y, groups, endpoint):
    """Choose regularization entirely within outer training groups."""
    grid = [.01, .1, 1.] if endpoint == "primary" else [1., 10., 100.]
    candidates = []
    cv = list(GroupKFold(n_splits=min(3, len(np.unique(groups)))).split(X, y, groups))
    for value in grid:
        group_losses = []
        for train, valid in cv:
            if endpoint == "primary" and len(np.unique(y[train])) < 2:
                prediction = np.full(len(valid), np.mean(y[train]))
            else:
                model = _fit(X[train], y[train], value, endpoint)
                prediction = (model.predict_proba(X[valid])[:, 1] if endpoint == "primary"
                              else model.predict(X[valid]))
            loss = (y[valid] - prediction) ** 2 if endpoint == "primary" else abs(y[valid] - prediction)
            group_losses.extend(np.mean(loss[groups[valid] == group]) for group in np.unique(groups[valid]))
        candidates.append(float(np.mean(group_losses)))
    return grid[int(np.argmin(candidates))]


def evaluate(frame, feature_sets, *, output_dir=None):
    """Predict every observed residue exactly once with its entire family held out.

    Return predictions, pooled metrics, fold metrics, coefficients and fold metadata.
    Explanations are additive standardized-feature contributions, not causal energies.
    """
    frame = frame.reset_index(drop=True)
    groups = frame.family_group.to_numpy()
    if len(np.unique(groups)) < 4:
        raise ValueError("Nested group evaluation requires at least four independent groups.")
    all_features = sorted({col for values in feature_sets.values() for col in values})
    if not np.isfinite(frame[all_features].to_numpy(float)).all():
        raise ValueError("Nonfinite structural predictors; fix feature coverage before evaluation.")
    y, energy = frame.hotspot.to_numpy(int), frame.Er.to_numpy(float)
    predictions, coefficients, folds = [], [], []
    meta = [col for col in ["complex_id", "family_group", "residue_key", "mutation", "wt", "side",
                            "Er", "hotspot", "label_conflict", "n_measurements"] if col in frame]
    for fold, (train, test) in enumerate(LeaveOneGroupOut().split(frame, y, groups)):
        assert not set(groups[train]) & set(groups[test])
        held = str(groups[test][0])
        for name, columns in {"prevalence_mean": [], **feature_sets}.items():
            row = frame.loc[test, meta].copy()
            row["model"], row["fold"] = name, fold
            details = {"fold": fold, "held_out_group": held, "model": name,
                       "train_groups": sorted(set(groups[train])), "n_train": len(train), "n_test": len(test)}
            if not columns:
                row["P_hotspot"], row["predicted_Er"] = np.mean(y[train]), np.mean(energy[train])
            else:
                X = frame[columns].to_numpy(float)
                for endpoint, target, outcol in [("primary", y, "P_hotspot"),
                                                  ("secondary", energy, "predicted_Er")]:
                    if endpoint == "primary" and len(np.unique(target[train])) < 2:
                        row[outcol] = float(target[train].mean())
                        details[endpoint + "_constant_train_class"] = True
                        continue
                    value = _tune(X[train], target[train], groups[train], endpoint)
                    model = _fit(X[train], target[train], value, endpoint)
                    details[endpoint + "_regularization"] = value
                    row[outcol] = (model.predict_proba(X[test])[:, 1] if endpoint == "primary"
                                   else model.predict(X[test]))
                    coef = np.asarray(model[-1].coef_).ravel()
                    intercept = float(np.asarray(model[-1].intercept_).ravel()[0])
                    transformed = model[0].transform(X[test])
                    contribution = transformed * coef
                    row[endpoint + "_intercept"] = intercept
                    row[endpoint + "_feature_contributions"] = [json.dumps(
                        {columns[j]: float(values[j]) for j in np.argsort(-abs(values))[:6]})
                        for values in contribution]
                    row[endpoint + "_other_feature_contribution"] = [
                        float(values.sum() - values[np.argsort(-abs(values))[:6]].sum()) for values in contribution]
                    coefficients.extend({"fold": fold, "held_out_group": held, "model": name,
                                         "endpoint": endpoint, "feature": col,
                                         "standardized_coefficient": float(c), "intercept": intercept,
                                         "training_feature_mean": float(mean), "training_feature_scale": float(scale)}
                                        for col, c, mean, scale in zip(columns, coef, model[0].mean_,
                                                                      model[0].scale_, strict=True))
            predictions.append(row)
            folds.append(details)
        print(f"completed outer fold {fold + 1}/{len(np.unique(groups))}: {held}", flush=True)
    pred = pd.concat(predictions, ignore_index=True)
    pooled, per_group = [], []
    for name, subset in pred.groupby("model", sort=False):
        pooled.append({"model": name, **metrics(subset.hotspot, subset.Er, subset.P_hotspot, subset.predicted_Er)})
        for group, part in subset.groupby("family_group"):
            per_group.append({"model": name, "family_group": group, "complex_id": ";".join(part.complex_id.unique()),
                              **metrics(part.hotspot, part.Er, part.P_hotspot, part.predicted_Er)})
    pooled, per_group, coefficients = pd.DataFrame(pooled), pd.DataFrame(per_group), pd.DataFrame(coefficients)
    if output_dir is not None:
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        for name, data in [("predictions", pred), ("metrics", pooled), ("group_metrics", per_group),
                           ("coefficients", coefficients)]:
            data.to_csv(output / (name + ".csv"), index=False)
        (output / "folds.json").write_text(json.dumps(folds, indent=2) + "\n")
    return pred, pooled, per_group, coefficients, folds


def bootstrap_comparisons(predictions, comparisons, *, n_bootstrap=1000, seed=42):
    """Paired family resampling; conditions on already fitted OOF predictions."""
    by_model = {name: part.set_index(["family_group", "complex_id", "residue_key"]).sort_index()
                for name, part in predictions.groupby("model")}
    result = []
    rng = np.random.default_rng(seed)
    for treatment, control in comparisons:
        left, right = by_model[treatment], by_model[control]
        if not left.index.equals(right.index):
            raise ValueError("Paired model comparison has different evaluation residues.")
        groups = left.index.get_level_values(0).unique().to_numpy()
        indices = {g: np.flatnonzero(left.index.get_level_values(0) == g) for g in groups}
        observed_l = metrics(left.hotspot, left.Er, left.P_hotspot, left.predicted_Er)
        observed_r = metrics(right.hotspot, right.Er, right.P_hotspot, right.predicted_Er)
        draws = {m: [] for m in ["AP", "AUROC", "Brier", "MAE", "Spearman"]}
        for _ in range(n_bootstrap):
            chosen = np.concatenate([indices[g] for g in rng.choice(groups, len(groups), replace=True)])
            left_sample, right_sample = left.iloc[chosen], right.iloc[chosen]
            lm = metrics(left_sample.hotspot, left_sample.Er, left_sample.P_hotspot, left_sample.predicted_Er)
            rm = metrics(right_sample.hotspot, right_sample.Er, right_sample.P_hotspot, right_sample.predicted_Er)
            for metric in draws:
                draws[metric].append(lm[metric] - rm[metric])
        for metric, samples in draws.items():
            finite = np.asarray(samples)[np.isfinite(samples)]
            lo, hi = np.quantile(finite, [.025, .975]) if len(finite) else [np.nan, np.nan]
            result.append({"treatment": treatment, "control": control, "metric": metric,
                           "difference": observed_l[metric] - observed_r[metric],
                           "lower_95": lo, "upper_95": hi, "valid_bootstrap": len(finite),
                           "independent_groups": len(groups),
                           "favorable_direction": "negative" if metric in ["Brier", "MAE"] else "positive"})
    return pd.DataFrame(result)


def calibration_table(predictions):
    rows = []
    for name, subset in predictions.groupby("model"):
        bins = np.minimum((subset.P_hotspot.to_numpy() * 5).astype(int), 4)
        for b in range(5):
            part = subset.iloc[np.flatnonzero(bins == b)]
            rows.append({"model": name, "bin_low": b / 5, "bin_high": (b + 1) / 5, "n": len(part),
                         "mean_probability": part.P_hotspot.mean(), "observed_frequency": part.hotspot.mean()})
    return pd.DataFrame(rows)
