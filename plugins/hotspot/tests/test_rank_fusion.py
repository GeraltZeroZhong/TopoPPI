import numpy as np
from topoppi_hotspot.rank_fusion import aggregate_ranks, partner_rank_percentiles


def test_fusion_uses_whole_candidate_query_and_ignores_context_only_nodes():
    queries = np.array(['a', 'a', 'a', 'a', 'b', 'b'])
    eligible = np.array([True, True, True, False, True, True])
    scores = np.array([1., 2., 2., 100., 5., 9.])
    percentiles = partner_rank_percentiles(scores, queries, eligible)
    np.testing.assert_allclose(percentiles[eligible], [1 / 6, 2 / 3, 2 / 3, .25, .75])
    assert np.isnan(percentiles[3])


def test_equal_rank_fusion_is_invariant_to_monotone_component_rescaling():
    queries = np.array(['a'] * 4)
    eligible = np.ones(4, bool)
    first, second = np.array([1., 4., 2., 3.]), np.array([4., 3., 1., 2.])
    baseline = aggregate_ranks([first, second], queries, eligible)
    transformed = aggregate_ranks([np.exp(first), second * 100 - 15], queries, eligible)
    np.testing.assert_allclose(baseline, transformed)
