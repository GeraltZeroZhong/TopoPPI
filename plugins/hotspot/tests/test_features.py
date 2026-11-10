"""Structural checks for independent, label-free Hotspot feature extraction."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from topoppi_hotspot.features import (
    FEATURE_SETS,
    extract_features,
    neighbor_means,
    shuffled_neighbor_means,
)

from topoppi.interface import prepare_interface

FIXTURE = Path(__file__).resolve().parent / "fixtures/tiny_complex.pdb"


@pytest.fixture(scope="module")
def result():
    return extract_features(FIXTURE, "A", "B")


def test_residue_identity_and_area_conservation(result):
    frame, graph = result
    assert frame.residue_key.tolist() == ["A:GLY:1", "A:ALA:2", "B:GLY:1", "B:ALA:2"]
    assert frame.wt.tolist() == ["G", "A", "G", "A"]
    assert frame.insertion_code.eq("").all()
    assert np.isfinite(frame[FEATURE_SETS["topology"]].to_numpy()).all()
    for side in (0, 1):
        rows = frame[frame.side == side]
        assert rows["surface.area"].sum() == pytest.approx(graph["metadata"]["surface_area_by_side"][side])
        assert rows["surface.footprint_area"].sum() == pytest.approx(graph["metadata"]["footprint_area_by_side"][side])
    assert graph["footprint_weights"].sum() > 0
    native = graph["footprint_edges"]
    assert np.all(graph["side"][native[:, 0]] == graph["side"][native[:, 1]])
    cross = graph["cross_edges"]
    assert np.all(graph["side"][cross[:, 0]] != graph["side"][cross[:, 1]])
    assert (native[:, 0] != native[:, 1]).all()


def test_neighbors_exclude_focal_residue_and_duplicates():
    values = np.array([[1.], [3.], [9.], [100.]])
    means, count = neighbor_means(values, [[0, 1], [1, 0], [1, 2], [1, 1]])
    np.testing.assert_array_equal(count, [1, 2, 1, 0])
    np.testing.assert_allclose(means[:, 0], [3, 5, 3, 0])


def test_shuffle_preserves_side_pool_and_never_uses_labels():
    values = np.array([[1.], [3.], [100.], [300.]])
    edges = [[0, 1], [2, 3]]
    means, counts = shuffled_neighbor_means(values, edges, [0, 0, 1, 1], [1, 1, 1, 1], seed=42)
    assert sorted(means[:2, 0]) == [1., 3.]
    assert sorted(means[2:, 0]) == [100., 300.]
    np.testing.assert_array_equal(counts, [1, 1, 1, 1])
    again, _ = shuffled_neighbor_means(values, edges, [0, 0, 1, 1], [1, 1, 1, 1], seed=42)
    np.testing.assert_array_equal(means, again)


def test_swapping_prepared_partners_changes_neither_residues_nor_features(result):
    frame, graph = result
    reversed_preparation = prepare_interface(FIXTURE, "B", "A")
    swapped, swapped_graph = extract_features(FIXTURE, "B", "A", prepared=reversed_preparation)
    pd.testing.assert_frame_equal(frame, swapped)
    for name in ("positions", "side", "cross_edges", "footprint_edges", "footprint_weights"):
        np.testing.assert_allclose(graph[name], swapped_graph[name])
