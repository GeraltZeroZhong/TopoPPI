import numpy as np
import pytest
from sklearn.base import clone
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from topoppi_hotspot.ranking_models import PairwiseHotspotRanker, _training_pairs, fit_pointwise


def test_pairs_stay_within_interface_and_receive_equal_interface_weight():
    X = np.array([[3.], [1.], [110.], [108.], [100.], [99.], [1000.]])
    y = np.array([1, 0, 1, 1, 0, 0, 0])
    queries = np.array(["a", "a", "b", "b", "b", "b", "excluded"])
    delta, targets, weights, summary = _training_pairs(X, y, queries)
    np.testing.assert_array_equal(delta[:2].ravel(), [2, -2])
    np.testing.assert_array_equal(delta[2:6].ravel(), [10, 11, 8, 9])
    np.testing.assert_array_equal(delta[6:].ravel(), [-10, -11, -8, -9])
    np.testing.assert_array_equal(targets, [1, 0, 1, 1, 1, 1, 0, 0, 0, 0])
    assert weights[:2].sum() == pytest.approx(.5)
    assert weights[2:].sum() == pytest.approx(.5)
    assert weights.sum() == pytest.approx(1.)
    assert summary["n_pairs"] == 5
    assert summary["n_informative_queries"] == 2
    assert summary["excluded_queries"][0]["query_id"] == "excluded"
    assert summary["excluded_queries"][0]["excluded_reason"] == "no_observed_positive"


def test_ranker_orders_training_pairs_and_cancels_interface_constant_feature():
    X = np.array([[0., 5.], [2., 5.], [3., 80.], [5., 80.]])
    model = PairwiseHotspotRanker(alpha=.1).fit(X, [0, 1, 0, 1], ["a", "a", "b", "b"])
    scores = model.decision_function(X)
    assert scores[1] > scores[0]
    assert scores[3] > scores[2]
    assert model.coef_[0, 1] == pytest.approx(0., abs=1e-12)
    assert model.intercept_ == 0
    assert model.training_summary_["converged"]
    assert clone(model).alpha == .1
    assert not hasattr(model, "predict_proba")


def test_fitted_within_interface_ranking_is_invariant_to_query_constant_offsets():
    X = np.array([[0., 2.], [2., 1.], [3., 3.], [5., 1.]])
    model = PairwiseHotspotRanker().fit(X, [0, 1, 0, 1], ["a", "a", "b", "b"])
    shifted = X.copy()
    shifted[:2] += [10., 50.]
    shifted[2:] += [-20., -10.]
    original = model.decision_function(X).reshape(2, 2)
    updated = model.decision_function(shifted).reshape(2, 2)
    np.testing.assert_allclose(np.diff(original, axis=1), np.diff(updated, axis=1))


def test_scaler_fits_training_residues_and_heldout_values_do_not_change_fit():
    X = np.array([[0.], [1.], [3.], [5.], [1000.], [2000.]])
    y = np.array([0, 1, 0, 1, 0, 1])
    queries = np.array(["a", "a", "b", "b", "held", "held"])
    train = queries != "held"
    before = PairwiseHotspotRanker().fit(X[train], y[train], queries[train])
    np.testing.assert_allclose(before.scaler_.mean_, X[train].mean(axis=0))
    np.testing.assert_allclose(before.scaler_.var_, X[train].var(axis=0))
    y[~train] = 1 - y[~train]
    X[~train] *= -10
    after = PairwiseHotspotRanker().fit(X[train], y[train], queries[train])
    np.testing.assert_allclose(before.coef_, after.coef_)
    np.testing.assert_allclose(before.scaler_.mean_, after.scaler_.mean_)


def test_one_class_interfaces_excluded_and_no_mixed_query_raises():
    X = np.arange(6).reshape(-1, 1)
    model = PairwiseHotspotRanker().fit(X, [0, 1, 1, 1, 0, 0], ["a", "a", "b", "b", "c", "c"])
    summary = model.training_summary_
    assert summary["n_informative_queries"] == 1
    assert summary["n_excluded_queries"] == 2
    assert {row["excluded_reason"] for row in summary["excluded_queries"]} == {
        "no_observed_positive", "no_observed_negative",
    }
    with pytest.raises(ValueError, match="both hotspot and nonhotspot"):
        PairwiseHotspotRanker().fit(X, [0, 0, 1, 1, 0, 0], ["a", "a", "b", "b", "c", "c"])


def test_pointwise_equal_interface_weights_match_explicit_fit():
    X = np.array([[0.], [1.], [2.], [2.5], [3.], [4.]])
    y = np.array([0, 1, 0, 0, 1, 1])
    queries = np.array(["a", "a", "b", "b", "b", "b"])
    model = fit_pointwise(X, y, queries, C=.3)
    transformed = StandardScaler().fit_transform(X)
    reference = LogisticRegression(C=.3, solver="lbfgs", max_iter=3000).fit(
        transformed, y, sample_weight=[1.5, 1.5, .75, .75, .75, .75])
    np.testing.assert_allclose(model[-1].coef_, reference.coef_)
    np.testing.assert_allclose(model[-1].intercept_, reference.intercept_)
    unweighted = fit_pointwise(X, y, queries, C=.3, equal_interface_weight=False)
    ordinary = LogisticRegression(C=.3, solver="lbfgs", max_iter=3000).fit(transformed, y)
    np.testing.assert_allclose(unweighted[-1].coef_, ordinary.coef_)
