"""Linear models for ranking measured residues within a complex interface.

The pairwise model returns ranking scores, not calibrated hotspot probabilities.
Only training residue rows fit the scaler; query labels construct training pairs.
"""
from __future__ import annotations

import warnings

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from sklearn.base import BaseEstimator
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.utils.validation import check_array, check_is_fitted, check_X_y


def _training_inputs(X, y, query_ids):
    X, y = check_X_y(X, y, dtype=float)
    queries = np.asarray(query_ids)
    if queries.ndim != 1 or len(queries) != len(y):
        raise ValueError("query_ids must contain one complex-interface ID per training residue.")
    if not np.isin(y, [0, 1]).all():
        raise ValueError("Hotspot training labels must be binary 0/1.")
    return X, y.astype(int), queries


def _training_pairs(X, y, query_ids):
    """Build both orientations, with equal total weight per informative query.

    A query with p positives and n negatives contributes 2*p*n oriented pairs.
    Each pair receives weight 1/(Q*2*p*n), where Q is the number of queries
    containing both classes. Thus all pair weights sum to one, independently
    of query size, label prevalence, and the number of training queries.
    """
    differences, targets, weights, rows = [], [], [], []
    for query in np.unique(query_ids):
        positive = np.flatnonzero((query_ids == query) & (y == 1))
        negative = np.flatnonzero((query_ids == query) & (y == 0))
        count = len(positive) * len(negative)
        row = {
            "query_id": query.item() if isinstance(query, np.generic) else query,
            "n_positive": len(positive), "n_negative": len(negative),
            "n_pairs": count, "n_oriented_pairs": 2 * count,
            "included": bool(count),
        }
        if count:
            delta = (X[positive, None, :] - X[negative][None, :, :]).reshape(count, X.shape[1])
            differences.extend([delta, -delta])
            targets.extend([np.ones(count), np.zeros(count)])
            weights.append(np.full(2 * count, 1. / (2 * count)))
        else:
            row["excluded_reason"] = "no_observed_positive" if not len(positive) else "no_observed_negative"
        rows.append(row)
    informative = sum(row["included"] for row in rows)
    if not informative:
        raise ValueError("Pairwise ranking requires a training interface with both hotspot and nonhotspot labels.")
    for row in rows:
        row["pair_weight_sum"] = 1. / informative if row["included"] else 0.
    summary = {
        "n_training_residues": len(y), "n_queries": len(rows),
        "n_informative_queries": informative, "n_excluded_queries": len(rows) - informative,
        "n_pairs": sum(row["n_pairs"] for row in rows),
        "n_oriented_pairs": sum(row["n_oriented_pairs"] for row in rows),
        "pair_weight_sum": 1., "queries": rows,
        "excluded_queries": [row for row in rows if not row["included"]],
    }
    return np.concatenate(differences), np.concatenate(targets), np.concatenate(weights) / informative, summary


class PairwiseHotspotRanker(BaseEstimator):
    """Fit a linear within-interface hotspot ranker with L2 regularization.

    The exact minimized objective is::

        mean_query(mean_positive_negative_pair(log(1 + exp(-w @ delta)))))
        + alpha / 2 * ||w||^2

    ``delta`` is the positive-minus-negative difference after fitting a
    StandardScaler on training residue rows. Both pair orientations enter the
    equivalent binary logistic loss, and no intercept is fitted. Interface
    constants therefore cancel in each pair. No pairs cross query boundaries.
    Interfaces with a single observed class contribute scaler rows but no
    pairwise labels; these exclusions are recorded in ``training_summary_``.

    Increasing alpha strengthens regularization. The average-loss definition
    keeps its meaning independent of the number of pairs and training queries.
    Callers must split by biological family before fitting, and supply complex
    interface IDs as query_ids. The class cannot infer either grouping.
    """

    def __init__(self, alpha=0.1, max_iter=1000, tol=1e-8):
        self.alpha = alpha
        self.max_iter = max_iter
        self.tol = tol

    def fit(self, X, y, query_ids):
        if not np.isfinite(self.alpha) or self.alpha <= 0:
            raise ValueError("alpha must be a finite positive L2 regularization strength.")
        X, y, queries = _training_inputs(X, y, query_ids)
        self.scaler_ = StandardScaler().fit(X)
        differences, targets, weights, summary = _training_pairs(self.scaler_.transform(X), y, queries)
        signs = 2 * targets - 1

        def objective(coef):
            margin = signs * (differences @ coef)
            loss = np.dot(weights, np.logaddexp(0, -margin)) + self.alpha / 2 * np.dot(coef, coef)
            gradient = differences.T @ (-weights * signs * expit(-margin)) + self.alpha * coef
            return float(loss), gradient

        result = minimize(objective, np.zeros(X.shape[1]), jac=True, method="L-BFGS-B",
                          options={"maxiter": self.max_iter, "gtol": self.tol, "ftol": self.tol})
        if not result.success:
            warnings.warn(f"Pairwise ranker did not converge: {result.message}", ConvergenceWarning, stacklevel=2)
        self.coef_ = result.x.reshape(1, -1)
        self.intercept_ = 0.
        self.n_features_in_ = X.shape[1]
        self.n_iter_ = np.asarray([result.nit])
        self.training_summary_ = {
            **summary, "alpha": float(self.alpha), "objective": float(result.fun),
            "converged": bool(result.success), "n_iterations": int(result.nit),
            "regularization": "mean_per_interface_pair_logloss_plus_alpha_over_2_L2",
        }
        return self

    def decision_function(self, X):
        """Return uncalibrated scores; compare scores within an interface."""
        check_is_fitted(self, ["coef_", "scaler_"])
        X = check_array(X, dtype=float)
        return self.scaler_.transform(X) @ self.coef_.ravel()


def fit_pointwise(X, y, query_ids, C=0.1, equal_interface_weight=True):
    """Fit a train-only scaled logistic comparator.

    When enabled, each of Q interfaces has total row weight n_train/Q:
    a row from interface i receives n_train/(Q*n_i), so weights sum to
    n_train and their overall mean is one. There is no extra class weighting.
    This pointwise comparator retains sklearn's usual inverse-regularization C.
    """
    X, y, queries = _training_inputs(X, y, query_ids)
    _, inverse, counts = np.unique(queries, return_inverse=True, return_counts=True)
    weights = len(y) / (len(counts) * counts[inverse]) if equal_interface_weight else np.ones(len(y))
    model = Pipeline([
        ("scaler", StandardScaler()),
        ("logisticregression", LogisticRegression(C=C, solver="lbfgs", max_iter=3000)),
    ])
    model.fit(X, y, logisticregression__sample_weight=weights)
    return model
