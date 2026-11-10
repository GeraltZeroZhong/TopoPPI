"""Exact pairwise linear ranking with train-fitted, balanced feature blocks."""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from sklearn.base import BaseEstimator
from sklearn.exceptions import ConvergenceWarning
from sklearn.preprocessing import StandardScaler
from sklearn.utils.validation import check_array, check_is_fitted

from .graph_ranking import _training_groups


class BlockPairwiseRanker(BaseEstimator):
    """Equal-query logistic ranker with a geometry and a full embedding block.

    Geometry columns are standardized on the supplied training rows, then
    divided by sqrt(number of geometry columns). The entire embedding block
    is divided by sqrt(mean_training(sum_over_columns(embedding**2))). Thus
    both blocks have a comparable mean squared row norm, without per-column
    whitening or PCA of the embeddings. A zero embedding block uses scale one.

    The caller must perform any label-free within-query embedding centering
    using its intended structural population before selecting training rows.
    This estimator never centers embeddings across training examples and never
    receives held-out rows during fit. Geometry may have zero columns.

    Full-batch L-BFGS minimizes the exact mean of each informative query's mean
    positive-negative log(1 + exp(-margin)), plus alpha/2 * ||coefficient||².
    There is no intercept, pair sampling, class reweighting, or early stopping
    selected by held-out labels. Scores are uncalibrated within-query ranks.
    """

    def __init__(self, alpha=0.1, max_iter=1000, tol=1e-8):
        self.alpha = alpha
        self.max_iter = max_iter
        self.tol = tol

    @staticmethod
    def _arrays(X_geometry, X_embedding):
        embedding = check_array(X_embedding, dtype=float)
        geometry = check_array(X_geometry, dtype=float, ensure_min_features=0)
        if len(geometry) != len(embedding):
            raise ValueError("Geometry and embedding blocks must have the same row count.")
        return geometry, embedding

    def _transform(self, geometry, embedding):
        if geometry.shape[1] != self.n_geometry_features_in_ or embedding.shape[1] != self.n_embedding_features_in_:
            raise ValueError("Feature block dimensions differ from the fitted model.")
        if self.n_geometry_features_in_:
            geometry = self.geometry_scaler_.transform(geometry) / np.sqrt(self.n_geometry_features_in_)
        return np.column_stack([geometry, embedding / self.embedding_scale_])

    def fit(self, X_geometry, X_embedding, y, query_ids):
        if not np.isfinite(self.alpha) or self.alpha <= 0:
            raise ValueError("alpha must be a finite positive L2 strength.")
        if (not isinstance(self.max_iter, (int, np.integer)) or self.max_iter < 1
                or not np.isfinite(self.tol) or self.tol <= 0):
            raise ValueError("max_iter must be a positive integer and tol must be finite and positive.")
        geometry, embedding = self._arrays(X_geometry, X_embedding)
        _, groups, summary = _training_groups(y, query_ids, np.arange(len(geometry)), len(geometry))
        self.n_geometry_features_in_, self.n_embedding_features_in_ = geometry.shape[1], embedding.shape[1]
        self.n_features_in_ = self.n_geometry_features_in_ + self.n_embedding_features_in_
        self.geometry_scaler_ = StandardScaler().fit(geometry) if self.n_geometry_features_in_ else None
        rms = float(np.sqrt(np.mean(np.square(embedding).sum(axis=1))))
        self.embedding_scale_ = rms if rms > 0 else 1.
        X = self._transform(geometry, embedding)
        differences, weights = [], []
        for positive, negative in groups:
            count = len(positive) * len(negative)
            differences.append((X[positive, None, :] - X[negative][None, :, :]).reshape(count, X.shape[1]))
            weights.append(np.full(count, 1. / (len(groups) * count)))
        differences, weights = np.concatenate(differences), np.concatenate(weights)

        def objective(coefficient):
            margin = differences @ coefficient
            loss = weights @ np.logaddexp(0., -margin) + self.alpha / 2 * (coefficient @ coefficient)
            gradient = -(differences.T @ (weights * expit(-margin))) + self.alpha * coefficient
            return float(loss), gradient

        result = minimize(objective, np.zeros(X.shape[1]), method="L-BFGS-B", jac=True,
                          options={"maxiter": self.max_iter, "ftol": self.tol, "gtol": self.tol})
        if not result.success:
            warnings.warn(f"Block ranker did not converge: {result.message}", ConvergenceWarning, stacklevel=2)
        self.coef_ = result.x.reshape(1, -1)
        self.intercept_ = 0.
        self.n_iter_ = np.asarray([result.nit])
        self.training_summary_ = {
            **summary, "alpha": self.alpha, "objective": float(result.fun),
            "converged": bool(result.success), "n_iterations": int(result.nit),
            "gradient_max_abs": float(np.max(np.abs(result.jac))),
            "query_weighting": "mean_of_informative_query_mean_pair_losses",
            "regularization": "exact_mean_pair_logloss_plus_alpha_over_2_coefficient_L2",
            "n_geometry_features": self.n_geometry_features_in_, "n_embedding_features": self.n_embedding_features_in_,
            "geometry_scaling": "training_StandardScaler_then_divide_sqrt_geometry_dimension",
            "embedding_scaling": "one_training_RMS_row_norm_scalar_without_PCA_or_per_dimension_whitening",
            "embedding_scale": self.embedding_scale_,
        }
        return self

    def decision_function(self, X_geometry, X_embedding):
        check_is_fitted(self, ["coef_", "embedding_scale_", "geometry_scaler_"])
        geometry, embedding = self._arrays(X_geometry, X_embedding)
        return self._transform(geometry, embedding) @ self.coef_.ravel()

    def save(self, directory):
        """Write plain JSON scaler/configuration and portable NumPy coefficients."""
        check_is_fitted(self, ["coef_", "embedding_scale_", "geometry_scaler_"])
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        scaler = None
        if self.geometry_scaler_ is not None:
            scaler = {key: getattr(self.geometry_scaler_, key + "_").tolist() for key in ("mean", "var", "scale")}
            scaler["n_samples_seen"] = int(self.geometry_scaler_.n_samples_seen_)
        metadata = {
            "estimator": type(self).__name__, "params": self.get_params(deep=False),
            "n_geometry_features": self.n_geometry_features_in_, "n_embedding_features": self.n_embedding_features_in_,
            "geometry_scaler": scaler, "embedding_scale": self.embedding_scale_,
            "training_summary": self.training_summary_,
        }
        (directory / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
        np.savez(directory / "weights.npz", coefficient=self.coef_.ravel())

    @classmethod
    def load(cls, directory):
        directory = Path(directory)
        metadata = json.loads((directory / "metadata.json").read_text())
        if metadata["estimator"] != cls.__name__:
            raise ValueError(f"Saved estimator is {metadata['estimator']}, not {cls.__name__}.")
        model = cls(**metadata["params"])
        model.n_geometry_features_in_ = metadata["n_geometry_features"]
        model.n_embedding_features_in_ = metadata["n_embedding_features"]
        model.n_features_in_ = model.n_geometry_features_in_ + model.n_embedding_features_in_
        model.embedding_scale_ = metadata["embedding_scale"]
        model.geometry_scaler_ = None
        if metadata["geometry_scaler"] is not None:
            model.geometry_scaler_ = StandardScaler()
            for key in ("mean", "var", "scale"):
                setattr(model.geometry_scaler_, key + "_", np.asarray(metadata["geometry_scaler"][key], dtype=float))
            model.geometry_scaler_.n_features_in_ = model.n_geometry_features_in_
            model.geometry_scaler_.n_samples_seen_ = metadata["geometry_scaler"]["n_samples_seen"]
        with np.load(directory / "weights.npz", allow_pickle=False) as weights:
            model.coef_ = weights["coefficient"].reshape(1, -1)
        model.intercept_ = 0.
        model.training_summary_ = metadata["training_summary"]
        model.n_iter_ = np.asarray([model.training_summary_["n_iterations"]])
        return model
