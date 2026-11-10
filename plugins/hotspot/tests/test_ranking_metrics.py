from itertools import permutations

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import average_precision_score
from topoppi_hotspot.ranking_metrics import (
    RANKING_METRICS,
    bootstrap_ranking_comparisons,
    ranking_metrics,
    summarize_rankings,
)


def test_ties_use_expected_top_k_relevance_and_binary_dcg():
    y = np.array([1, 1, 0, 0, 0, 0])
    result = ranking_metrics(y, np.ones(6))
    discount = 1 / np.log2(np.arange(6) + 2)
    assert result["AP"] == pytest.approx(1 / 3)
    assert result["precision_at_3"] == pytest.approx(1 / 3)
    assert result["precision_at_5"] == pytest.approx(1 / 3)
    assert result["recall_at_5"] == pytest.approx(5 / 6)
    assert result["NDCG_at_5"] == pytest.approx(discount[:5].sum() / 3 / discount[:2].sum())
    assert result["NDCG"] == pytest.approx(discount.sum() / 3 / discount[:2].sum())
    assert result["enrichment_top_10pct"] == pytest.approx(1.)
    random_ap = np.mean([average_precision_score(order, -np.arange(6)) for order in set(permutations(y))])
    assert result["random_expected_AP"] == pytest.approx(random_ap)
    assert result["random_expected_AP"] > result["AP"]


def test_partial_tie_at_cutoff_is_not_broken_by_input_order():
    result = ranking_metrics([1, 1, 0, 0, 0, 1], [4, 3, 2, 2, 2, 1])
    permuted = ranking_metrics([1, 0, 1, 0, 0, 1], [4, 2, 3, 2, 2, 1])
    assert result == permuted
    assert result["precision_at_3"] == pytest.approx(2 / 3)


def test_no_hotspots_and_small_queries_have_explicit_coverage():
    result = ranking_metrics([0, 0], [1, 0])
    assert result["precision_at_3"] == result["precision_at_5"] == 0
    assert result["k3_effective"] == result["k5_effective"] == 2
    assert not result["has_k3"] and not result["has_k5"]
    for metric in ["AP", "recall_at_5", "NDCG_at_5", "NDCG", "enrichment_top_10pct"]:
        assert np.isnan(result[metric])
    result = ranking_metrics([1, 0], [1, 0])
    assert result["precision_at_5"] == .5
    assert result["recall_at_5"] == 1
    assert result["enrichment_top_10pct"] == 2


def _predictions():
    return pd.DataFrame({
        "model": ["one"] * 8, "complex_id": ["a"] * 2 + ["b"] * 6,
        "family_group": ["fa"] * 2 + ["fb"] * 6, "side": [0, 1] * 4,
        "residue_key": list("abcdefgh"), "hotspot": [1, 0, 0, 0, 0, 0, 0, 1],
        "rank_score": [1, 0, 6, 5, 4, 3, 2, 1],
    })


def test_macro_ranking_is_not_pooled_and_ignores_interface_score_shift():
    frame = _predictions()
    per_query, summary = summarize_rankings(frame)
    assert len(per_query) == 2  # Both partner sides comprise one query.
    assert summary.loc[0, "macro_AP"] == pytest.approx((1 + 1 / 6) / 2)
    pooled = average_precision_score(frame.hotspot, frame.rank_score)
    assert summary.loc[0, "macro_AP"] != pytest.approx(pooled)
    frame.loc[frame.complex_id == "a", "rank_score"] += 100
    shifted, shifted_summary = summarize_rankings(frame)
    pd.testing.assert_frame_equal(per_query, shifted)
    pd.testing.assert_frame_equal(summary, shifted_summary)
    assert pooled != pytest.approx(average_precision_score(frame.hotspot, frame.rank_score))


def test_unmeasured_residues_are_excluded_not_labeled_negative():
    frame = _predictions()
    frame.loc[frame.complex_id == "a", "hotspot"] = 0
    frame.loc[2, "hotspot"] = np.nan
    per_query, summary = summarize_rankings(frame)
    assert summary.loc[0, "n_queries"] == 2
    assert summary.loc[0, "n_positive_queries"] == 1
    assert summary.loc[0, "n_AP"] == 1
    assert summary.loc[0, "n_precision_at_5"] == 2
    assert summary.loc[0, "n_queries_at_least_5"] == 1
    assert summary.loc[0, "n_residues"] == 7
    assert per_query.loc[per_query.complex_id == "b", "n_unmeasured"].item() == 1


def test_paired_bootstrap_rejects_query_and_residue_mismatch():
    frame = _predictions()
    second = frame.assign(model="two")
    per_query, _ = summarize_rankings(pd.concat([frame, second.iloc[2:]]))
    with pytest.raises(ValueError, match="different evaluation queries"):
        bootstrap_ranking_comparisons(per_query, [("one", "two")], n_bootstrap=10)
    second.loc[0, "residue_key"] = "changed"
    per_query, _ = summarize_rankings(pd.concat([frame, second]))
    with pytest.raises(ValueError, match="different evaluation evaluated_residue_keys"):
        bootstrap_ranking_comparisons(per_query, [("one", "two")], n_bootstrap=10)


def test_paired_bootstrap_uses_whole_families_and_macro_query_weights():
    frame = _predictions()
    third = frame.loc[frame.complex_id == "a"].assign(complex_id="c")
    frame = pd.concat([frame, third], ignore_index=True)
    second = frame.assign(model="two", rank_score=-frame.rank_score)
    per_query, _ = summarize_rankings(pd.concat([frame, second]))
    result = bootstrap_ranking_comparisons(per_query, [("one", "two")], n_bootstrap=100)
    result_again = bootstrap_ranking_comparisons(per_query, [("one", "two")], n_bootstrap=100)
    pd.testing.assert_frame_equal(result, result_again)
    ap = result.loc[result.metric == "AP"].iloc[0]
    assert ap.n_queries == 3 and ap.independent_groups == 2
    # Family fa contributes two queries, each +0.5; fb contributes one, -5/6.
    assert ap.difference == pytest.approx((.5 + .5 - 5 / 6) / 3)
    assert ap.lower_95 == pytest.approx(-5 / 6)
    assert ap.upper_95 == pytest.approx(.5)
    assert set(result.metric) == set(RANKING_METRICS)
