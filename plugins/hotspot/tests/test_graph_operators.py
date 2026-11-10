"""Independent synthetic checks for reusable numerical primitives."""
import numpy as np
import pandas as pd
import pytest
from topoppi_hotspot import graph_operators


def _graph():
    # Full registry includes a non-interface row, shuffled relative to frame.
    return {
        "residue_keys": np.array(["outside", "right", "left", "alanine"]),
        "footprint_edges": np.array([[0, 2], [2, 3]]), "footprint_weights": np.array([99., 4.]),
        "spatial_edges": np.array([[2, 3]]),
        "cross_edges": np.array([[1, 2], [1, 3]]), "cross_weights": np.array([1., 3.]),
    }

def test_graph_mapping_retains_unlabelled_nodes_and_complex_boundaries():
    frame = pd.DataFrame({
        "complex_id": ["b", "a", "b", "a", "a", "b"],
        "residue_key": ["right", "left", "alanine", "right", "alanine", "left"],
        "side": [1, 0, 0, 1, 0, 0], "wt": ["X", "X", "A", "X", "A", "X"],
        "family_group": ["fb", "fa", "fb", "fa", "fa", "fb"],
    }, index=[20, 10, 70, 90, 50, 30])
    operators, metadata = graph_operators.build_graph_operators(frame, {"a": _graph(), "b": _graph()})
    assert metadata["n_nodes"] == 6
    assert metadata["n_native_alanine_nodes"] == 2
    expected = np.zeros((6, 6))
    expected[0, [5, 2]] = [.25, .75]
    expected[3, [1, 4]] = [.25, .75]
    expected[[5, 2], 0] = 1.
    expected[[1, 4], 3] = 1.
    np.testing.assert_allclose(operators["context"].toarray(), expected)
    expected_local = np.zeros((6, 6))
    expected_local[1, 4] = expected_local[4, 1] = 1.
    expected_local[2, 5] = expected_local[5, 2] = 1.
    for name in ("native", "spatial"):
        np.testing.assert_allclose(operators[name].toarray(), expected_local)
    # Row permutation must change only matrix ordering, including interleaved
    # graphs and frame indexes that are not consecutive positional indices.
    permutation = [5, 4, 3, 2, 1, 0]
    permuted, _ = graph_operators.build_graph_operators(frame.iloc[permutation], {"a": _graph(), "b": _graph()})
    for name in operators:
        np.testing.assert_array_equal(permuted[name].toarray(), operators[name].toarray()[permutation][:, permutation])

def test_mismatched_graph_key_does_not_silently_drop_interface_nodes():
    frame = pd.DataFrame({"complex_id": ["a"], "residue_key": ["missing"], "side": [0]})
    with pytest.raises(ValueError, match="Every interface row"):
        graph_operators.build_graph_operators(frame, {"a": _graph()})
