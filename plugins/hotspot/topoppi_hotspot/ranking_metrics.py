"""Equal-interface ranking endpoints for experimentally measured residues.

AP uses sklearn's threshold/tie-block convention. Top-k precision/recall and
DCG average over uniformly random orders inside exact score ties; these are
expectations, not one arbitrary or seeded tie break. Binary relevance is used
for both NDCG endpoints. Unmeasured residues never become negative labels.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

RANKING_METRICS = (
    "AP", "precision_at_3", "precision_at_5", "recall_at_5",
    "NDCG_at_5", "NDCG", "enrichment_top_10pct",
)


def ranking_metrics(y, score):
    """Evaluate one measured interface, combining both partner sides.

    Effective k is min(k, n_measured). AP, recall, NDCG and enrichment are
    undefined without observed positives. Precision remains zero for an
    all-negative interface. Empty measured sets have no defined endpoint.
    ``random_expected_AP`` is the exact expectation under a uniform random
    permutation, which differs from AP of a constant-score predictor.
    """
    y, score = np.asarray(y, dtype=float), np.asarray(score, dtype=float)
    if y.ndim != 1 or score.shape != y.shape:
        raise ValueError("Labels and scores must be aligned one-dimensional arrays.")
    if not np.isfinite(y).all() or not np.isin(y, [0, 1]).all():
        raise ValueError("Ranking requires measured binary labels; omit unmeasured residues.")
    if not np.isfinite(score).all():
        raise ValueError("Ranking scores must be finite.")
    n, positives = len(y), int(y.sum())
    result = {
        "n": n, "hotspots": positives,
        "prevalence": positives / n if n else np.nan,
        "k3_effective": min(3, n), "k5_effective": min(5, n),
        "k10pct_effective": int(np.ceil(.1 * n)),
        "has_k3": n >= 3, "has_k5": n >= 5,
        "random_expected_AP": np.nan,
        **{metric: np.nan for metric in RANKING_METRICS},
    }
    if not n:
        return result

    order = np.argsort(-score, kind="stable")
    sorted_y, sorted_score = y[order], score[order]
    # Every position in a tied block has the same expected binary relevance.
    boundaries = np.r_[0, np.flatnonzero(np.diff(sorted_score)) + 1, n]
    expected_relevance = np.empty(n)
    for start, end in zip(boundaries[:-1], boundaries[1:], strict=True):
        expected_relevance[start:end] = sorted_y[start:end].mean()
    k3, k5, k10 = (result[key] for key in ["k3_effective", "k5_effective", "k10pct_effective"])
    result["precision_at_3"] = float(expected_relevance[:k3].mean())
    result["precision_at_5"] = float(expected_relevance[:k5].mean())
    if not positives:
        return result

    discount = 1 / np.log2(np.arange(n) + 2)
    result.update({
        "AP": float(average_precision_score(y, score)),
        "recall_at_5": float(expected_relevance[:k5].sum() / positives),
        "NDCG_at_5": float(np.dot(expected_relevance[:k5], discount[:k5])
                           / discount[:min(k5, positives)].sum()),
        "NDCG": float(np.dot(expected_relevance, discount) / discount[:positives].sum()),
        "enrichment_top_10pct": float(expected_relevance[:k10].mean() / (positives / n)),
    })
    harmonic = np.sum(1 / np.arange(1, n + 1))
    result["random_expected_AP"] = float(
        harmonic / n + (positives - 1) * (n - harmonic) / (n * (n - 1)) if n > 1 else 1.
    )
    return result


def summarize_rankings(predictions, query_columns=("complex_id",), score_col="rank_score"):
    """Return per-query metrics and equal-query means, separately by model.

    ``complex_id`` is the default query: both partner sides belong to the same
    ranking. Alternative queries, e.g. (complex_id, side), must be explicit.
    Missing hotspot labels are excluded and their count retained. Each summary
    endpoint reports its own finite-query denominator. Family membership and
    evaluated residue identities, when available, accompany the per-query rows
    so paired comparisons can reject mismatched evaluation cohorts.
    """
    query_columns = list(query_columns)
    rows = []
    for keys, part in predictions.groupby(["model", *query_columns], sort=False, dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        row = dict(zip(["model", *query_columns], keys, strict=True))
        measured = part.loc[part.hotspot.notna()]
        row["n_unmeasured"] = len(part) - len(measured)
        if "family_group" in part:
            families = part.family_group.unique()
            if len(families) != 1 or pd.isna(families[0]):
                raise ValueError("Each ranking query must belong to exactly one family_group.")
            row["family_group"] = families[0]
        if "residue_key" in part:
            if measured.residue_key.duplicated().any():
                raise ValueError("A measured residue occurs more than once in a ranking query.")
            ordered = measured.sort_values("residue_key")
            row["evaluated_residue_keys"] = tuple(ordered.residue_key)
            row["evaluated_labels"] = tuple(ordered.hotspot.astype(int))
        row.update(ranking_metrics(measured.hotspot, measured[score_col]))
        rows.append(row)
    per_query = pd.DataFrame(rows)
    summaries = []
    if per_query.empty:
        return per_query, pd.DataFrame()
    for model, part in per_query.groupby("model", sort=False):
        summary = {
            "model": model, "n_queries": int((part.n > 0).sum()),
            "n_empty_queries": int((part.n == 0).sum()),
            "n_positive_queries": int((part.hotspots > 0).sum()),
            "n_zero_positive_queries": int(((part.n > 0) & (part.hotspots == 0)).sum()),
            "n_residues": int(part.n.sum()), "n_unmeasured": int(part.n_unmeasured.sum()),
            "n_queries_at_least_3": int(part.has_k3.sum()),
            "n_queries_at_least_5": int(part.has_k5.sum()),
            "mean_k3_effective": part.loc[part.n > 0, "k3_effective"].mean(),
            "mean_k5_effective": part.loc[part.n > 0, "k5_effective"].mean(),
            "n_families": part.loc[part.n > 0, "family_group"].nunique()
            if "family_group" in part else np.nan,
        }
        for metric in [*RANKING_METRICS, "random_expected_AP"]:
            summary["macro_" + metric] = part[metric].mean()
            summary["n_" + metric] = int(part[metric].notna().sum())
        summaries.append(summary)
    return per_query, pd.DataFrame(summaries)


def bootstrap_ranking_comparisons(
    per_query, comparisons, *, query_columns=("complex_id",), n_bootstrap=2000, seed=42,
):
    """Paired family-cluster bootstrap of equal-interface metric differences.

    Each sampled family includes every query in that family; repeated draws
    retain multiplicity. Percentile intervals condition on the already fitted
    OOF predictions and do not represent uncertainty from refitting models.
    All requested endpoints are returned, including undefined comparisons.
    """
    if n_bootstrap < 1:
        raise ValueError("n_bootstrap must be positive.")
    by_model = {}
    for model, part in per_query.groupby("model", sort=False):
        if part.duplicated(list(query_columns)).any():
            raise ValueError("Duplicate evaluation queries in paired model comparison.")
        by_model[model] = part.set_index(list(query_columns)).sort_index()
    rows = []
    rng = np.random.default_rng(seed)
    for treatment, control in comparisons:
        left, right = by_model[treatment], by_model[control]
        if not left.index.equals(right.index):
            raise ValueError("Paired model comparison has different evaluation queries.")
        for column in ["family_group", "n", "hotspots", "evaluated_residue_keys", "evaluated_labels"]:
            if column in left and not left[column].equals(right[column]):
                raise ValueError(f"Paired model comparison has different evaluation {column}.")
        families = left.family_group.unique()
        indices = [np.flatnonzero(left.family_group.to_numpy() == family) for family in families]
        draws = rng.integers(len(families), size=(n_bootstrap, len(families)))
        for metric in RANKING_METRICS:
            difference = left[metric].to_numpy(float) - right[metric].to_numpy(float)
            finite = np.isfinite(difference)
            sums = np.asarray([np.sum(difference[idx][finite[idx]]) for idx in indices])
            counts = np.asarray([finite[idx].sum() for idx in indices])
            denominators = counts[draws].sum(axis=1)
            valid = denominators > 0
            samples = sums[draws].sum(axis=1)[valid] / denominators[valid]
            lower, upper = np.quantile(samples, [.025, .975]) if len(samples) else [np.nan, np.nan]
            rows.append({
                "treatment": treatment, "control": control, "metric": metric,
                "difference": float(difference[finite].mean()) if finite.any() else np.nan,
                "lower_95": lower, "upper_95": upper, "valid_bootstrap": int(valid.sum()),
                "n_bootstrap": n_bootstrap, "independent_groups": len(families),
                "n_queries": len(left), "n_metric_queries": int(finite.sum()),
                "n_metric_groups": int((counts > 0).sum()), "favorable_direction": "positive",
            })
    return pd.DataFrame(rows)
