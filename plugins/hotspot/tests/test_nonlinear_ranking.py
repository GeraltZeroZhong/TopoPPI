import numpy as np
import pytest
from sklearn.base import clone
from topoppi_hotspot.nonlinear_ranking import (
    TorchPairwiseMLP,
    XGBPairwiseRanker,
    _query_groups,
    _query_pair_loss,
)


def _ranker(kind):
    if kind == "xgb":
        pytest.importorskip("xgboost")
        return XGBPairwiseRanker(n_estimators=25, random_state=7)
    pytest.importorskip("torch")
    return TorchPairwiseMLP(epochs=30, random_state=7)


@pytest.mark.parametrize("kind", ["xgb", "torch"])
def test_training_scaler_includes_single_class_queries_but_never_heldout(kind):
    model = _ranker(kind)
    X = np.array([[0.], [2.], [3.], [5.], [20.], [22.]])
    y = [0, 1, 0, 1, 1, 1]
    queries = ["a", "a", "b", "b", "single", "single"]
    model.fit(X, y, queries)
    np.testing.assert_allclose(model.scaler_.mean_, X.mean(axis=0))
    np.testing.assert_allclose(model.scaler_.var_, X.var(axis=0))
    assert model.training_summary_["n_informative_queries"] == 2
    assert model.training_summary_["excluded_queries"][0]["query_id"] == "single"
    before = model.decision_function(X)
    model.decision_function([[1e6], [-1e6]])
    np.testing.assert_array_equal(model.decision_function(X), before)
    np.testing.assert_allclose(model.scaler_.mean_, X.mean(axis=0))
    assert clone(model).get_params() == model.get_params()


@pytest.mark.parametrize("kind", ["xgb", "torch"])
def test_separate_single_class_queries_cannot_create_cross_query_pairs(kind):
    # Globally both classes exist, but there are no legitimate within-query pairs.
    with pytest.raises(ValueError, match="both hotspot and nonhotspot"):
        _ranker(kind).fit([[0.], [1.], [5.], [6.]], [0, 0, 1, 1], ["a", "a", "b", "b"])


def test_torch_loss_ignores_query_offsets_and_gives_equal_query_weight():
    torch = pytest.importorskip("torch")
    # One easy pair in a and six hard pairs in b. Neither cross-query offsets
    # nor the number of pairs should alter the relative weighting of a and b.
    y = np.array([1, 0, 1, 1, 0, 0, 0])
    groups, _ = _query_groups(y, np.array(["a", "a", "b", "b", "b", "b", "b"]))
    scores = torch.tensor([2., 0., -1., -1., 1., 1., 1.], requires_grad=True)
    loss = _query_pair_loss(scores, groups)
    expected = (np.logaddexp(0, -2) + np.logaddexp(0, 2)) / 2
    assert float(loss.detach()) == pytest.approx(expected)
    shifted = scores.detach().clone()
    shifted[:2] += 1000
    shifted[2:] -= 1000
    assert float(_query_pair_loss(shifted, groups)) == pytest.approx(expected)
    loss.backward()
    assert float(scores.grad[:2].sum()) == pytest.approx(0., abs=1e-7)
    assert float(scores.grad[2:].sum()) == pytest.approx(0., abs=1e-7)


@pytest.mark.parametrize("objective", ["rank:pairwise", "rank:map"])
def test_xgb_receives_sorted_aligned_queries_and_one_weight_per_query(monkeypatch, objective):
    xgb = pytest.importorskip("xgboost")
    original_fit = xgb.XGBRanker.fit
    seen = {}

    def capture_fit(self, X, y, **kwargs):
        seen.update(X=X.copy(), y=np.array(y), qid=np.array(kwargs["qid"]),
                    weights=np.array(kwargs["sample_weight"]))
        return original_fit(self, X, y, **kwargs)

    monkeypatch.setattr(xgb.XGBRanker, "fit", capture_fit)
    X = np.array([[5.], [0.], [3.], [2.], [4.], [20.]])
    y = np.array([1, 0, 0, 1, 0, 1])
    q = ["z", "a", "z", "a", "z", "single"]
    # MAP gradients on two/three-row queries can be too small for the default
    # minimum leaf Hessian. Let this integration toy form its intended splits.
    model = XGBPairwiseRanker(objective=objective, n_estimators=40, min_child_weight=.01).fit(X, y, q)
    order = [1, 3, 0, 2, 4]
    np.testing.assert_allclose(seen["X"], model.scaler_.transform(X)[order])
    np.testing.assert_array_equal(seen["y"], y[order])
    np.testing.assert_array_equal(seen["qid"], [0, 0, 1, 1, 1])
    np.testing.assert_array_equal(seen["weights"], [1, 1])
    scores = model.decision_function(X)
    assert scores[3] > scores[1]
    assert scores[0] > max(scores[2], scores[4])
    permutation = [4, 2, 0, 5, 1, 3]
    np.testing.assert_array_equal(model.decision_function(X[permutation]), scores[permutation])


def test_mlp_learns_a_nonlinear_ordering_with_fixed_seed():
    torch = pytest.importorskip("torch")
    # Each interface contains XOR labels, which no linear readout can separate.
    corners = np.array([[-2., -2.], [-2., 2.], [2., -2.], [2., 2.]])
    X = np.tile(corners, (3, 1))
    y = np.tile([1, 0, 0, 1], 3)
    queries = np.repeat(["a", "b", "c"], 4)
    before_rng = torch.random.get_rng_state().clone()
    model = TorchPairwiseMLP(alpha=.001, random_state=4).fit(X, y, queries)
    assert torch.equal(before_rng, torch.random.get_rng_state())
    unseen = corners * .9
    scores = model.decision_function(unseen)
    assert min(scores[0], scores[3]) > max(scores[1], scores[2])
    assert model.training_summary_["final_pair_logloss"] < .3
    same = clone(model).fit(X, y, queries)
    np.testing.assert_array_equal(model.decision_function(unseen), same.decision_function(unseen))


@pytest.mark.parametrize("kind", ["xgb", "torch"])
def test_portable_model_roundtrip_restores_scaler_and_predictions(tmp_path, kind):
    X = np.array([[-2., 1.], [2., 1.], [-1., 4.], [3., 4.]])
    model = _ranker(kind).fit(X, [0, 1, 0, 1], ["a", "a", "b", "b"])
    heldout = np.array([[4., 3.], [-3., 8.], [0., 2.]])
    directory = tmp_path / kind
    model.save(directory)
    restored = type(model).load(directory)
    assert restored.get_params() == model.get_params()
    assert restored.training_summary_ == model.training_summary_
    np.testing.assert_array_equal(restored.scaler_.mean_, model.scaler_.mean_)
    np.testing.assert_array_equal(restored.scaler_.scale_, model.scaler_.scale_)
    np.testing.assert_array_equal(restored.decision_function(heldout), model.decision_function(heldout))
