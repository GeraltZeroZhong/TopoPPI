"""Small CPU ranking heads for externally split, measured hotspot rows.

Both heads scale all training rows, but learn rankings only from interfaces with
both observed classes. Callers own family splits and hyperparameter selection.
Scores are uncalibrated and must be compared within an interface.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from sklearn.base import BaseEstimator
from sklearn.preprocessing import StandardScaler
from sklearn.utils.validation import check_array, check_is_fitted

from .ranking_models import _training_inputs


def _query_groups(y, queries):
    groups, rows = [], []
    for query in np.unique(queries):
        positive = np.flatnonzero((queries == query) & (y == 1))
        negative = np.flatnonzero((queries == query) & (y == 0))
        count = len(positive) * len(negative)
        row = {
            "query_id": query.item() if isinstance(query, np.generic) else query,
            "n_positive": len(positive), "n_negative": len(negative),
            "n_pairs": count, "included": bool(count),
        }
        if count:
            groups.append((positive, negative))
        else:
            row["excluded_reason"] = "no_observed_positive" if not len(positive) else "no_observed_negative"
        rows.append(row)
    if not groups:
        raise ValueError("Pairwise ranking requires a training interface with both hotspot and nonhotspot labels.")
    return groups, {
        "n_training_residues": len(y), "n_queries": len(rows),
        "n_informative_queries": len(groups), "n_excluded_queries": len(rows) - len(groups),
        "n_pairs": sum(row["n_pairs"] for row in rows),
        "queries": rows, "excluded_queries": [row for row in rows if not row["included"]],
    }


def _query_pair_loss(scores, groups):
    """Average all positive-negative softplus losses, then average queries."""
    import torch
    from torch.nn.functional import softplus

    return torch.stack([
        softplus(scores[negative][None, :] - scores[positive][:, None]).mean()
        for positive, negative in groups
    ]).mean()


class _PortableRanker:
    def _save_metadata(self, directory):
        check_is_fitted(self, ["model_", "scaler_"])
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        metadata = {
            "estimator": type(self).__name__, "params": self.get_params(deep=False),
            "scaler": {
                "mean": self.scaler_.mean_.tolist(), "var": self.scaler_.var_.tolist(),
                "scale": self.scaler_.scale_.tolist(),
                "n_samples_seen": int(self.scaler_.n_samples_seen_),
            },
            "training_summary": self.training_summary_,
        }
        (directory / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
        return directory

    @classmethod
    def _load_metadata(cls, directory):
        directory = Path(directory)
        metadata = json.loads((directory / "metadata.json").read_text(encoding="utf-8"))
        if metadata["estimator"] != cls.__name__:
            raise ValueError(f"Saved estimator is {metadata['estimator']}, not {cls.__name__}.")
        model = cls(**metadata["params"])
        model.scaler_ = StandardScaler()
        for name in ("mean", "var", "scale"):
            setattr(model.scaler_, name + "_", np.asarray(metadata["scaler"][name], dtype=float))
        model.n_features_in_ = len(model.scaler_.mean_)
        model.scaler_.n_features_in_ = model.n_features_in_
        model.scaler_.n_samples_seen_ = metadata["scaler"]["n_samples_seen"]
        model.training_summary_ = metadata["training_summary"]
        return model, directory


class XGBPairwiseRanker(_PortableRanker, BaseEstimator):
    """Original XGBRanker with stable query sorting and uniform group weights.

    ``rank:pairwise`` and binary ``rank:map`` retain XGBoost's own sampled-pair
    LambdaMART objectives and normalizations. Equal external query weights do
    not make these objectives identical to the MLP's exhaustive mean-pair loss.
    No row-level weights or early-stopping validation data are passed.
    """

    def __init__(self, objective="rank:pairwise", max_depth=2, n_estimators=100,
                 learning_rate=0.05, reg_lambda=1.0, min_child_weight=1.0,
                 lambdarank_num_pair_per_sample=8, random_state=0, n_jobs=1):
        self.objective = objective
        self.max_depth = max_depth
        self.n_estimators = n_estimators
        self.learning_rate = learning_rate
        self.reg_lambda = reg_lambda
        self.min_child_weight = min_child_weight
        self.lambdarank_num_pair_per_sample = lambdarank_num_pair_per_sample
        self.random_state = random_state
        self.n_jobs = n_jobs

    def _new_model(self):
        from xgboost import XGBRanker

        return XGBRanker(
            **self.get_params(deep=False), tree_method="hist", device="cpu",
            lambdarank_pair_method="mean", lambdarank_normalization=True,
            lambdarank_score_normalization=True, verbosity=0,
        )

    def fit(self, X, y, query_ids):
        if self.objective not in {"rank:pairwise", "rank:map"}:
            raise ValueError("objective must be rank:pairwise or rank:map.")
        X, y, queries = _training_inputs(X, y, query_ids)
        groups, summary = _query_groups(y, queries)
        self.scaler_ = StandardScaler().fit(X)
        # Stable ordering within each query; labels remain aligned with rows.
        ordered_groups = [np.sort(np.concatenate(group)) for group in groups]
        order = np.concatenate(ordered_groups)
        qid = np.repeat(np.arange(len(groups)), [len(group) for group in ordered_groups])
        self.model_ = self._new_model()
        self.model_.fit(self.scaler_.transform(X)[order], y[order], qid=qid,
                        sample_weight=np.ones(len(groups)), verbose=False)
        self.n_features_in_ = X.shape[1]
        self.training_summary_ = {
            **summary, "objective": self.objective,
            "query_weighting": "one_equal_sample_weight_per_informative_query",
            "pair_construction": "xgboost_mean_sampled_pairs",
            "n_pairs_description": "available_positive_negative_pairs_not_sampled_count",
            "n_estimators": self.n_estimators,
        }
        return self

    def decision_function(self, X):
        check_is_fitted(self, ["model_", "scaler_"])
        # Prediction requires no query sorting, and preserves the supplied order.
        return self.model_.predict(self.scaler_.transform(check_array(X, dtype=float)))

    def save(self, directory):
        """Save XGBoost's portable UBJSON model and plain-JSON scaler/params."""
        directory = self._save_metadata(directory)
        self.model_.save_model(directory / "model.ubj")

    @classmethod
    def load(cls, directory):
        model, directory = cls._load_metadata(directory)
        model.model_ = model._new_model()
        model.model_.load_model(directory / "model.ubj")
        return model


class TorchPairwiseMLP(_PortableRanker, BaseEstimator):
    """One hidden GELU layer trained with an exact equal-query pairwise loss.

    The objective is mean_query(mean_positive_negative_pair(softplus(-margin)))
    + alpha/2 * sum(parameter**2), including the hidden bias. The output has no
    bias because it cancels in pairwise scores. Full-batch Adam uses a fixed
    epoch count; no held-out rows select an epoch. Torch is imported lazily and
    every tensor stays on CPU. ``random_state`` controls initialization/dropout
    without changing the caller's CPU random-generator state after fitting.
    """

    def __init__(self, hidden_dim=32, dropout=0.2, alpha=0.01, epochs=200,
                 learning_rate=0.001, random_state=0):
        self.hidden_dim = hidden_dim
        self.dropout = dropout
        self.alpha = alpha
        self.epochs = epochs
        self.learning_rate = learning_rate
        self.random_state = random_state

    def _new_model(self):
        from torch import nn

        return nn.Sequential(
            nn.Linear(self.n_features_in_, self.hidden_dim), nn.GELU(),
            nn.Dropout(self.dropout), nn.Linear(self.hidden_dim, 1, bias=False),
        ).cpu()

    def fit(self, X, y, query_ids):
        import torch

        if not np.isfinite(self.alpha) or self.alpha < 0:
            raise ValueError("alpha must be a finite nonnegative L2 strength.")
        if not 0 <= self.dropout < 1:
            raise ValueError("dropout must lie in [0, 1).")
        if self.hidden_dim < 1 or self.epochs < 1 or not np.isfinite(self.learning_rate) or self.learning_rate <= 0:
            raise ValueError("hidden_dim, epochs and learning_rate must be positive.")
        X, y, queries = _training_inputs(X, y, query_ids)
        groups, summary = _query_groups(y, queries)
        self.scaler_ = StandardScaler().fit(X)
        self.n_features_in_ = X.shape[1]
        inputs = torch.as_tensor(self.scaler_.transform(X), dtype=torch.float32, device="cpu")
        pairs = [(torch.as_tensor(pos), torch.as_tensor(neg)) for pos, neg in groups]
        with torch.random.fork_rng(devices=[]):
            torch.random.default_generator.manual_seed(self.random_state)
            self.model_ = self._new_model()
            optimizer = torch.optim.Adam(self.model_.parameters(), lr=self.learning_rate)
            self.loss_curve_ = []
            self.model_.train()
            for _ in range(self.epochs):
                optimizer.zero_grad()
                scores = self.model_(inputs).flatten()
                pair_loss = _query_pair_loss(scores, pairs)
                penalty = self.alpha / 2 * sum(p.square().sum() for p in self.model_.parameters())
                loss = pair_loss + penalty
                loss.backward()
                optimizer.step()
                self.loss_curve_.append(float(loss.detach()))
            self.model_.eval()
            with torch.no_grad():
                final_loss = float(_query_pair_loss(self.model_(inputs).flatten(), pairs))
        self.training_summary_ = {
            **summary, "objective": "mean_per_interface_pair_logloss_plus_alpha_over_2_L2",
            "alpha": self.alpha, "epochs": self.epochs, "final_pair_logloss": final_loss,
            "query_weighting": "mean_of_informative_query_mean_pair_losses",
        }
        return self

    def decision_function(self, X):
        import torch

        check_is_fitted(self, ["model_", "scaler_"])
        values = self.scaler_.transform(check_array(X, dtype=float))
        self.model_.eval()
        with torch.no_grad():
            return self.model_(torch.as_tensor(values, dtype=torch.float32, device="cpu")).flatten().numpy()

    def save(self, directory):
        """Save CPU parameter arrays and plain-JSON architecture/scaler/params."""
        directory = self._save_metadata(directory)
        np.savez(directory / "weights.npz", **{
            key: value.detach().cpu().numpy() for key, value in self.model_.state_dict().items()
        })

    @classmethod
    def load(cls, directory):
        import torch

        model, directory = cls._load_metadata(directory)
        with torch.random.fork_rng(devices=[]):
            model.model_ = model._new_model()
        with np.load(directory / "weights.npz", allow_pickle=False) as weights:
            model.model_.load_state_dict({key: torch.from_numpy(weights[key].copy()) for key in weights.files})
        model.model_.eval()
        return model
