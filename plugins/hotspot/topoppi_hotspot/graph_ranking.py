"""One-layer relational mean-aggregation ranker for externally split graphs.

The caller supplies fixed row-normalized local/context operators and owns all
group splits. Features may include unlabelled nodes, but only ``train_indices``
contribute labels or fit the scaler. No edges or label-derived inputs are made
here. A block-diagonal graph stays disconnected throughout feature aggregation.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from sklearn.base import BaseEstimator
from sklearn.preprocessing import StandardScaler
from sklearn.utils.validation import check_array, check_is_fitted

from .nonlinear_ranking import _PortableRanker, _query_pair_loss


def _training_groups(y, query_ids, train_indices, n_nodes):
    """Read labels only at training indices and construct within-query pairs."""
    indices = np.asarray(train_indices)
    if indices.ndim != 1 or not np.issubdtype(indices.dtype, np.integer) or not len(indices):
        raise ValueError("train_indices must be a nonempty one-dimensional integer array.")
    if np.any(indices < 0) or np.any(indices >= n_nodes) or len(np.unique(indices)) != len(indices):
        raise ValueError("train_indices must be unique indices within X.")
    labels, queries = np.asarray(y), np.asarray(query_ids)
    if labels.ndim != 1 or len(labels) != n_nodes:
        raise ValueError("y must have one entry per node; non-training entries are ignored.")
    if queries.ndim != 1 or len(queries) != n_nodes:
        raise ValueError("query_ids must have one entry per node.")
    # Do not cast or validate held-out entries: they may be missing or arbitrary.
    train_labels = np.asarray(labels[indices], dtype=float)
    train_queries = queries[indices]
    if not np.isin(train_labels, [0, 1]).all():
        raise ValueError("Labels at train_indices must be binary 0/1.")
    groups, rows = [], []
    for query in np.unique(train_queries):
        positive = indices[(train_queries == query) & (train_labels == 1)]
        negative = indices[(train_queries == query) & (train_labels == 0)]
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
        raise ValueError("Pairwise ranking requires a training query with both binary classes.")
    return indices, groups, {
        "n_training_nodes": len(indices), "n_queries": len(rows),
        "n_informative_queries": len(groups), "n_excluded_queries": len(rows) - len(groups),
        "n_pairs": sum(row["n_pairs"] for row in rows),
        "queries": rows, "excluded_queries": [row for row in rows if not row["included"]],
    }


def _neighbor_mean(adjacency, values, name):
    """Apply a supplied sparse/dense mean operator without densifying a graph."""
    if sparse.issparse(adjacency):
        operator = adjacency.astype(np.float32, copy=False).tocsr()
        weights = operator.data
    else:
        operator = check_array(adjacency, dtype=np.float32)
        weights = operator
    if operator.shape != (len(values), len(values)):
        raise ValueError(f"{name} must be a square operator with one row per node.")
    if not np.isfinite(weights).all() or np.any(weights < 0):
        raise ValueError(f"{name} must contain finite nonnegative weights.")
    row_sums = np.asarray(operator.sum(axis=1)).ravel()
    if not np.all(np.isclose(row_sums, 1., atol=1e-5) | (row_sums == 0)):
        raise ValueError(f"{name} must be row-normalized, with zero rows for isolated nodes.")
    return np.asarray(operator @ values, dtype=np.float32)


class GraphRanker(_PortableRanker, BaseEstimator):
    """GraphSAGE-style relation means with a residual own-node projection.

    For standardized numeric features z, the relational encoder computes::

        h = GELU(W_self z + W_local (A_local z) + W_context (A_context z) + b)
        score = w_out Dropout(h)

    This adapts GraphSAGE's concatenated-mean aggregation (Hamilton et al.,
    NeurIPS 2017, Algorithm 1 and Section 3.3) with separate local and context
    relations, one hidden layer, and no embedding normalization:
    https://papers.nips.cc/paper_files/paper/2017/file/5dd9db5e033da9c6fb5ba83c7a7ebea9-Paper.pdf
    ``local`` omits the context branch; ``identity`` omits both graph branches
    and has the same function class as ``TorchPairwiseMLP`` at matching width.
    Each node's own features have an independent linear path into the hidden
    representation. There is no output bias because it cancels in ranking.

    The fixed-epoch full-batch CPU objective is the exact mean, over informative
    training queries, of each query's mean positive-negative softplus loss,
    plus alpha/2 times the sum of squared parameters (including hidden bias).
    Queries with only one training class fit the scaler but create no pairs.
    No validation labels or early stopping are used. All graph nodes may carry
    features, including unlabelled neighbors. Supplied operators must therefore
    already enforce the caller's permitted connectivity and split boundaries.

    X is a finite dense numeric array. A_local and A_context may be NumPy dense
    arrays or SciPy sparse matrices/arrays, with row sums one or zero. Inactive
    ablation operators are ignored and may be None. Prediction requires fresh
    operators aligned to X; graphs and labels are not stored in the model.
    """

    def __init__(self, hidden_dim=32, dropout=0.2, alpha=0.01, epochs=200,
                 learning_rate=0.001, random_state=0, mode="relational"):
        self.hidden_dim = hidden_dim
        self.dropout = dropout
        self.alpha = alpha
        self.epochs = epochs
        self.learning_rate = learning_rate
        self.random_state = random_state
        self.mode = mode

    def _new_model(self):
        from torch import nn

        class RelationalMeanEncoder(nn.Module):
            def __init__(self, n_features, width, dropout, mode):
                super().__init__()
                self.own = nn.Linear(n_features, width)
                self.local = nn.Linear(n_features, width, bias=False) if mode != "identity" else None
                self.context = nn.Linear(n_features, width, bias=False) if mode == "relational" else None
                self.activation = nn.GELU()
                self.dropout = nn.Dropout(dropout)
                self.output = nn.Linear(width, 1, bias=False)

            def forward(self, own, local=None, context=None):
                hidden = self.own(own)
                if self.local is not None:
                    hidden = hidden + self.local(local)
                if self.context is not None:
                    hidden = hidden + self.context(context)
                return self.output(self.dropout(self.activation(hidden))).flatten()

        return RelationalMeanEncoder(self.n_features_in_, self.hidden_dim, self.dropout, self.mode).cpu()

    def _inputs(self, X, A_local, A_context):
        import torch

        values = self.scaler_.transform(X).astype(np.float32)
        local = _neighbor_mean(A_local, values, "A_local") if self.mode != "identity" else None
        context = _neighbor_mean(A_context, values, "A_context") if self.mode == "relational" else None
        return tuple(torch.as_tensor(v, dtype=torch.float32, device="cpu") if v is not None else None
                     for v in (values, local, context))

    def fit(self, X, y, query_ids, train_indices, A_local, A_context):
        """Fit using only training labels; all supplied feature nodes may aggregate."""
        import torch

        if self.mode not in {"relational", "local", "identity"}:
            raise ValueError("mode must be relational, local or identity.")
        if not np.isfinite(self.alpha) or self.alpha < 0:
            raise ValueError("alpha must be a finite nonnegative L2 strength.")
        if not 0 <= self.dropout < 1:
            raise ValueError("dropout must lie in [0, 1).")
        if (not isinstance(self.hidden_dim, (int, np.integer)) or self.hidden_dim < 1
                or not isinstance(self.epochs, (int, np.integer)) or self.epochs < 1
                or not np.isfinite(self.learning_rate) or self.learning_rate <= 0):
            raise ValueError("hidden_dim and epochs must be positive integers; learning_rate must be positive.")
        X = check_array(X, dtype=float)
        indices, groups, summary = _training_groups(y, query_ids, train_indices, len(X))
        self.scaler_ = StandardScaler().fit(X[indices])
        self.n_features_in_ = X.shape[1]
        # Fixed inputs and a single graph layer make these means constant across
        # epochs; precomputation is algebraically identical to A @ X per forward.
        inputs = self._inputs(X, A_local, A_context)
        pairs = [(torch.as_tensor(pos), torch.as_tensor(neg)) for pos, neg in groups]
        with torch.random.fork_rng(devices=[]):
            torch.random.default_generator.manual_seed(self.random_state)
            self.model_ = self._new_model()
            optimizer = torch.optim.Adam(self.model_.parameters(), lr=self.learning_rate)
            self.loss_curve_ = []
            self.model_.train()
            for _ in range(self.epochs):
                optimizer.zero_grad()
                pair_loss = _query_pair_loss(self.model_(*inputs), pairs)
                penalty = self.alpha / 2 * sum(p.square().sum() for p in self.model_.parameters())
                loss = pair_loss + penalty
                loss.backward()
                optimizer.step()
                self.loss_curve_.append(float(loss.detach()))
            self.model_.eval()
            with torch.no_grad():
                final_loss = float(_query_pair_loss(self.model_(*inputs), pairs))
        self.training_summary_ = {
            **summary, "n_graph_nodes": len(X), "mode": self.mode,
            "objective": "mean_per_query_pair_logloss_plus_alpha_over_2_L2",
            "query_weighting": "mean_of_informative_query_mean_pair_losses",
            "alpha": self.alpha, "epochs": self.epochs, "final_pair_logloss": final_loss,
            "scaler_fit": "training_nodes_only", "device": "cpu",
        }
        return self

    def decision_function(self, X, A_local, A_context):
        """Return uncalibrated node scores in X order; compare within queries."""
        import torch

        check_is_fitted(self, ["model_", "scaler_"])
        inputs = self._inputs(check_array(X, dtype=float), A_local, A_context)
        self.model_.eval()
        with torch.no_grad():
            return self.model_(*inputs).numpy()

    def save(self, directory):
        """Save CPU parameter arrays and plain-JSON architecture/scaler metadata."""
        directory = self._save_metadata(directory)
        np.savez(directory / "weights.npz", **{
            key: value.detach().cpu().numpy() for key, value in self.model_.state_dict().items()
        })

    @classmethod
    def load(cls, directory):
        """Restore portable weights and train-fitted scaler on CPU."""
        import torch

        model, directory = cls._load_metadata(directory)
        with torch.random.fork_rng(devices=[]):
            model.model_ = model._new_model()
        with np.load(directory / "weights.npz", allow_pickle=False) as weights:
            model.model_.load_state_dict({key: torch.from_numpy(weights[key].copy()) for key in weights.files})
        model.model_.eval()
        return model
