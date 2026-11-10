import numpy as np
import pytest
from scipy import sparse
from sklearn.base import clone
from topoppi_hotspot.graph_ranking import GraphRanker, _training_groups
from topoppi_hotspot.nonlinear_ranking import TorchPairwiseMLP, _query_pair_loss

torch = pytest.importorskip("torch")


@pytest.fixture
def graph_data():
    # Two disconnected graphs and an isolated node, with artificial features.
    X = np.array([[-2., 1.], [2., 1.], [-1., 4.], [3., 4.],
                  [-4., -2.], [4., -2.], [1000., 2000.]])
    local = sparse.block_diag([
        np.array([[0., 1., 0., 0.], [1., 0., 0., 0.],
                  [0., 0., 0., 1.], [0., 0., 1., 0.]]),
        np.array([[0., 1.], [1., 0.]]), np.zeros((1, 1)),
    ], format="csr")
    context = sparse.block_diag([
        np.array([[0., 0., 1., 0.], [0., 0., 0., 1.],
                  [1., 0., 0., 0.], [0., 1., 0., 0.]]),
        np.array([[0., 1.], [1., 0.]]), np.zeros((1, 1)),
    ], format="csr")
    y = np.array([0., 1., 0., 1., np.nan, np.nan, 1.])
    queries = np.array(["a", "a", "b", "b", "heldout", "heldout", "single"])
    return X, y, queries, np.array([0, 1, 2, 3, 6]), local, context


def _ranker(**kwargs):
    return GraphRanker(epochs=12, random_state=7, **kwargs)


@pytest.mark.parametrize("dense", [False, True])
def test_node_and_graph_permutation_equivariance(graph_data, dense):
    X, y, queries, train, local, context = graph_data
    if dense:
        local, context = local.toarray(), context.toarray()
    model = _ranker().fit(X, y, queries, train, local, context)
    scores = model.decision_function(X, local, context)
    # Swap the graphs and reorder nodes within the first graph.
    order = np.array([5, 4, 6, 2, 0, 3, 1])
    permuted = model.decision_function(X[order], local[order][:, order], context[order][:, order])
    np.testing.assert_allclose(permuted, scores[order], atol=2e-7, rtol=2e-6)


def test_heldout_labels_are_ignored_and_scaler_uses_only_training_nodes(graph_data):
    X, y, queries, train, local, context = graph_data
    before_rng = torch.random.get_rng_state().clone()
    model = _ranker().fit(X, y, queries, train, local, context)
    assert torch.equal(before_rng, torch.random.get_rng_state())
    changed = y.astype(object)
    changed[4:6] = ["not even a label", None]
    other = clone(model).fit(X, changed, queries, train, local, context)
    np.testing.assert_array_equal(model.decision_function(X, local, context),
                                  other.decision_function(X, local, context))
    np.testing.assert_allclose(model.scaler_.mean_, X[train].mean(axis=0))
    np.testing.assert_allclose(model.scaler_.var_, X[train].var(axis=0))
    assert model.scaler_.n_samples_seen_ == len(train)
    assert model.training_summary_["n_informative_queries"] == 2
    assert model.training_summary_["excluded_queries"][0]["query_id"] == "single"
    assert len(model.loss_curve_) == model.epochs


def test_disconnected_graph_features_cannot_mix_into_training_or_predictions(graph_data):
    X, y, queries, train, local, context = graph_data
    model = _ranker().fit(X, y, queries, train, local, context)
    changed = X.copy()
    changed[4:6] = [[1e6, -1e6], [-1e6, 1e6]]
    other = clone(model).fit(changed, y, queries, train, local, context)
    original = model.decision_function(X, local, context)
    np.testing.assert_array_equal(other.decision_function(X, local, context), original)
    np.testing.assert_array_equal(model.decision_function(changed, local, context)[train], original[train])
    # Predicting a component alone gives exactly the same graph-local function.
    isolated = model.decision_function(X[:4], local[:4, :4], context[:4, :4])
    np.testing.assert_allclose(isolated, original[:4], atol=2e-7, rtol=2e-6)


def test_exact_pair_loss_has_equal_query_weight_and_no_heldout_gradient():
    # Query a has one pair; query b has six. Last two nodes are held out even
    # though their query IDs coincide with a training query.
    y = [1, 0, 1, 1, 0, 0, 0, 1, 0]
    queries = np.array(["a", "a", "b", "b", "b", "b", "b", "a", "a"])
    _, groups, summary = _training_groups(y, queries, np.arange(7), len(y))
    scores = torch.tensor([2., 0., -1., -1., 1., 1., 1., -100., 100.], requires_grad=True)
    loss = _query_pair_loss(scores, groups)
    expected = (np.logaddexp(0., -2.) + np.logaddexp(0., 2.)) / 2
    assert float(loss.detach()) == pytest.approx(expected)
    assert summary["n_pairs"] == 7
    loss.backward()
    np.testing.assert_array_equal(scores.grad[-2:].numpy(), 0.)
    assert float(scores.grad[:2].sum()) == pytest.approx(0., abs=1e-7)
    assert float(scores.grad[2:7].sum()) == pytest.approx(0., abs=1e-7)


def test_separate_single_class_queries_cannot_create_cross_query_pairs():
    with pytest.raises(ValueError, match="both binary classes"):
        _ranker(mode="identity").fit([[0.], [1.], [5.], [6.]], [0, 0, 1, 1],
                                     ["a", "a", "b", "b"], [0, 1, 2, 3], None, None)


def test_identity_matches_existing_mlp_at_same_seed_and_width():
    # Matching forward values and training prove that this ablation does not
    # secretly add another layer, residual output term, or graph contribution.
    X = np.array([[-2., 1.], [2., 1.], [-1., 4.], [3., 4.]])
    y, queries = [0, 1, 0, 1], ["a", "a", "b", "b"]
    graph = _ranker(mode="identity").fit(X, y, queries, np.arange(4), None, None)
    mlp = TorchPairwiseMLP(epochs=12, random_state=7).fit(X, y, queries)
    np.testing.assert_array_equal(graph.decision_function(X, None, None), mlp.decision_function(X))
    np.testing.assert_array_equal(graph.loss_curve_, mlp.loss_curve_)


def test_local_ablation_ignores_context_but_uses_unlabelled_neighbor_features():
    X = np.array([[0., 0.], [0., 0.], [-2., 1.], [2., 1.]])
    local = sparse.csr_matrix(([1., 1.], ([0, 1], [2, 3])), shape=(4, 4))
    y, queries = [0, 1, np.nan, np.nan], ["a", "a", "unused", "unused"]
    model = GraphRanker(mode="local", epochs=100, random_state=7).fit(X, y, queries, [0, 1], local, None)
    scores = model.decision_function(X, local, None)
    assert scores[1] > scores[0] + 1
    np.testing.assert_array_equal(model.decision_function(X, local, np.ones((4, 4))), scores)
    swapped = X.copy()
    swapped[[2, 3]] = swapped[[3, 2]]
    np.testing.assert_allclose(model.decision_function(swapped, local, None)[:2], scores[[1, 0]])


@pytest.mark.parametrize("mode", ["relational", "local", "identity"])
def test_portable_roundtrip_preserves_predictions_and_scaler(tmp_path, graph_data, mode):
    X, y, queries, train, local, context = graph_data
    model = _ranker(mode=mode).fit(X, y, queries, train, local, context)
    model.save(tmp_path / mode)
    restored = GraphRanker.load(tmp_path / mode)
    assert restored.get_params() == model.get_params()
    assert restored.training_summary_ == model.training_summary_
    np.testing.assert_array_equal(restored.scaler_.mean_, model.scaler_.mean_)
    np.testing.assert_array_equal(restored.scaler_.scale_, model.scaler_.scale_)
    np.testing.assert_array_equal(restored.decision_function(X, local, context),
                                  model.decision_function(X, local.toarray(), context.toarray()))


def test_unnormalized_operator_is_rejected(graph_data):
    X, y, queries, train, local, context = graph_data
    with pytest.raises(ValueError, match="row-normalized"):
        _ranker().fit(X, y, queries, train, local * 2, context)
