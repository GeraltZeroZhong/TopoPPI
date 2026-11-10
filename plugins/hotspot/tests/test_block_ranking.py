import numpy as np
import pytest
from sklearn.base import clone
from topoppi_hotspot.block_ranking import BlockPairwiseRanker


def _data():
    geometry = np.array([[-2., 1.], [2., 1.], [-1., 4.], [3., 4.], [20., -2.]])
    embedding = np.array([[-3., 1., 0.], [3., -1., 0.], [-2., 1., 0.], [2., -1., 0.], [0., 0., 0.]])
    return geometry, embedding, np.array([0, 1, 0, 1, 1]), np.array(["a", "a", "b", "b", "single"])


def test_train_only_scaling_has_declared_block_norms_and_ignores_prediction_data():
    geometry, embedding, y, queries = _data()
    model = BlockPairwiseRanker().fit(geometry, embedding, y, queries)
    np.testing.assert_allclose(model.geometry_scaler_.mean_, geometry.mean(axis=0))
    assert model.embedding_scale_ == pytest.approx(np.sqrt(np.mean(np.sum(embedding ** 2, axis=1))))
    transformed = model._transform(geometry, embedding)
    assert np.mean(np.sum(transformed[:, :2] ** 2, axis=1)) == pytest.approx(1.)
    assert np.mean(np.sum(transformed[:, 2:] ** 2, axis=1)) == pytest.approx(1.)
    assert model.training_summary_["n_training_nodes"] == len(y)
    assert model.training_summary_["n_informative_queries"] == 2
    before = model.decision_function(geometry, embedding)
    model.decision_function(geometry * 10000, embedding * 10000)
    np.testing.assert_array_equal(model.decision_function(geometry, embedding), before)
    assert clone(model).get_params() == model.get_params()


def test_embedding_units_and_node_order_do_not_change_the_fitted_function():
    geometry, embedding, y, queries = _data()
    model = BlockPairwiseRanker().fit(geometry, embedding, y, queries)
    rescaled = clone(model).fit(geometry, embedding * 1000, y, queries)
    scores = model.decision_function(geometry, embedding)
    np.testing.assert_allclose(rescaled.decision_function(geometry, embedding * 1000), scores, atol=1e-10)
    order = np.array([4, 2, 1, 3, 0])
    np.testing.assert_allclose(model.decision_function(geometry[order], embedding[order]), scores[order])
    reordered = clone(model).fit(geometry[order], embedding[order], y[order], queries[order])
    np.testing.assert_allclose(reordered.decision_function(geometry, embedding), scores, atol=1e-8)


def test_exact_equal_query_objective_and_optimum_gradient():
    # One positive-negative pair in query a, six in query b, and an excluded
    # single-class query. Reconstruct the independently weighted objective.
    embedding = np.array([[-2., 1.], [2., -1.], [-1., 0.], [-2., -1.], [1., 2.], [2., 1.], [3., 0.], [20., -4.]])
    y = np.array([0, 1, 0, 0, 1, 1, 1, 1])
    queries = np.array(["a", "a", "b", "b", "b", "b", "b", "single"])
    model = BlockPairwiseRanker(alpha=.2, tol=1e-12).fit(np.empty((8, 0)), embedding, y, queries)
    standardized = embedding / model.embedding_scale_

    def objective(coefficient):
        scores = standardized @ coefficient
        first = np.logaddexp(0., scores[0] - scores[1])
        second = np.logaddexp(0., scores[2:4, None] - scores[None, 4:7]).mean()
        return (first + second) / 2 + .2 / 2 * np.dot(coefficient, coefficient)

    coefficient = model.coef_.ravel()
    assert model.training_summary_["objective"] == pytest.approx(objective(coefficient), abs=1e-12)
    steps = np.eye(len(coefficient)) * 1e-5
    gradient = [(objective(coefficient + step) - objective(coefficient - step)) / 2e-5 for step in steps]
    np.testing.assert_allclose(gradient, 0., atol=1e-7)


@pytest.mark.parametrize("geometry_dimension", [0, 2])
def test_portable_roundtrip_with_and_without_geometry(tmp_path, geometry_dimension):
    geometry, embedding, y, queries = _data()
    geometry = geometry[:, :geometry_dimension]
    model = BlockPairwiseRanker().fit(geometry, embedding, y, queries)
    model.save(tmp_path)
    restored = BlockPairwiseRanker.load(tmp_path)
    np.testing.assert_array_equal(restored.decision_function(geometry, embedding),
                                  model.decision_function(geometry, embedding))
    assert restored.get_params() == model.get_params()
    assert restored.training_summary_ == model.training_summary_


def test_zero_embedding_block_is_finite_and_single_class_queries_do_not_cross_pair():
    geometry, embedding, y, queries = _data()
    model = BlockPairwiseRanker().fit(geometry, np.zeros_like(embedding), y, queries)
    assert model.embedding_scale_ == 1.
    assert np.isfinite(model.decision_function(geometry, np.zeros_like(embedding))).all()
    with pytest.raises(ValueError, match="both binary classes"):
        BlockPairwiseRanker().fit(geometry, embedding, y, ["negative", "positive", "negative", "positive", "positive"])
