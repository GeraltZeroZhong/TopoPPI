"""Graph regressions that would corrupt geometry, mutation mapping, or ablations."""
import json
from pathlib import Path

import numpy as np
import pytest

from topoppi_affinity.graph_features import (
    FEATURE_NAMES,
    extract_complex_graph,
    load_graph,
    mutant_node_features,
    save_graph,
)

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture(scope="module")
def graph():
    return extract_complex_graph(FIXTURES / "tiny_complex.pdb", "A", "B", max_surface_points=128)


def test_surface_mass_topological_edges_and_complete_graph(graph):
    assert set(graph["chain"]) == {"A", "B"}
    assert graph["x"].shape == (len(graph["residue_key"]), len(FEATURE_NAMES))
    assert np.isfinite(graph["x"]).all()
    assert graph["point_pos"].shape[0] <= 256
    area = graph["x"][:, FEATURE_NAMES.index("surface.area")].sum()
    footprint = graph["x"][:, FEATURE_NAMES.index("surface.footprint_area")].sum()
    assert np.isclose(graph["point_weight"].sum(), area, rtol=1e-5)
    assert np.isclose(graph["point_interface_weight"].sum(), footprint, rtol=1e-5)
    edges = graph["edge_index"]
    surface = graph["edge_attr"][:, 6] > 0
    assert surface.any()
    assert np.all(graph["side"][edges[0, surface]] == graph["side"][edges[1, surface]])
    assert {tuple(edge) for edge in edges.T} == {tuple(edge[::-1]) for edge in edges.T}
    assert not any("ddg" in str(key).lower() or "pkd" in str(key).lower() for key in graph)


def test_mutations_change_chemistry_only_joint_sites_and_noop(graph):
    changed, identities, mask = mutant_node_features(graph, "GA1W,AA2D")
    assert mask.sum() == 2
    np.testing.assert_array_equal(changed[:, 26:], graph["x"][:, 26:])
    assert np.count_nonzero(identities != graph["aa_index"]) == 2
    unchanged, _, noop = mutant_node_features(graph, "GA1G")
    np.testing.assert_array_equal(unchanged, graph["x"])
    assert not noop.any()
    with pytest.raises(ValueError, match="WT identity mismatch"):
        mutant_node_features(graph, "WA1F")


def test_partner_swap_and_portable_npz(graph, tmp_path):
    swapped = extract_complex_graph(FIXTURES / "tiny_complex.pdb", "B", "A", max_surface_points=128)
    for key in ("x", "pos", "edge_index", "edge_attr", "point_pos"):
        np.testing.assert_allclose(swapped[key], graph[key])
    path = tmp_path / "graph.npz"
    save_graph(graph, path)
    loaded = load_graph(path)
    np.testing.assert_array_equal(loaded["residue_key"], graph["residue_key"])
    assert json.loads(str(loaded["metadata_json"]))["geometry_protocol"] == "fixed_WT_identity_substitution"
