"""Regression summaries with complexes, rather than mutation rows, as groups."""

from __future__ import annotations

import numpy as np
from scipy.stats import pearsonr, spearmanr


def regression_metrics(y_true, y_pred, groups=None, *, n_bootstrap=0, seed=42):
    truth, predicted = np.asarray(y_true, dtype=float), np.asarray(y_pred, dtype=float)
    if truth.ndim != 1 or predicted.shape != truth.shape or len(truth) == 0:
        raise ValueError("Predictions and observations must be aligned nonempty vectors")
    if not np.isfinite(truth).all() or not np.isfinite(predicted).all():
        raise ValueError("Cannot evaluate nonfinite predictions or observations")
    residual = predicted - truth
    result = {"n": len(truth), "mae": float(np.mean(np.abs(residual))),
              "rmse": float(np.sqrt(np.mean(residual ** 2))), "bias": float(np.mean(residual))}
    variable = len(truth) >= 3 and np.ptp(truth) > 0 and np.ptp(predicted) > 0
    result["pearson"] = float(pearsonr(truth, predicted).statistic) if variable else None
    result["spearman"] = float(spearmanr(truth, predicted).statistic) if variable else None
    if groups is None:
        return result
    labels = np.asarray(groups).astype(str)
    if labels.shape != truth.shape:
        raise ValueError("Group labels must align with predictions")
    unique = np.unique(labels)
    indices = [np.flatnonzero(labels == group) for group in unique]
    result["n_groups"] = len(unique)
    result["macro_mae"] = float(np.mean([np.mean(np.abs(residual[idx])) for idx in indices]))
    result["macro_rmse"] = float(np.mean([np.sqrt(np.mean(residual[idx] ** 2)) for idx in indices]))
    within = [float(spearmanr(truth[idx], predicted[idx]).statistic) for idx in indices
              if len(idx) >= 3 and np.ptp(truth[idx]) > 0 and np.ptp(predicted[idx]) > 0]
    result["within_group_spearman"] = float(np.mean(within)) if within else None
    result["within_group_evaluable"] = len(within)
    if n_bootstrap and len(unique) >= 2:
        generator = np.random.default_rng(seed)
        values = []
        for _ in range(n_bootstrap):
            sampled = generator.integers(len(indices), size=len(indices))
            errors = residual[np.concatenate([indices[index] for index in sampled])]
            values.append([np.mean(np.abs(errors)), np.sqrt(np.mean(errors ** 2))])
        intervals = np.quantile(values, [0.025, 0.975], axis=0)
        result["group_bootstrap_95_ci"] = {"mae": intervals[:, 0].tolist(), "rmse": intervals[:, 1].tolist()}
    return result


def mutation_strata(y_true, y_pred, counts, groups):
    counts, groups = np.asarray(counts), np.asarray(groups)
    truth, predicted = np.asarray(y_true), np.asarray(y_pred)
    result = {}
    for name, selected in (("single", counts == 1), ("double", counts == 2),
                           ("triple", counts == 3), ("four_or_more", counts >= 4)):
        if np.any(selected):
            result[name] = regression_metrics(truth[selected], predicted[selected], groups[selected])
    return result

