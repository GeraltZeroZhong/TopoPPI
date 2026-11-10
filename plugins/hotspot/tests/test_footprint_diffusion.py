import numpy as np
import pandas as pd
import pytest
from topoppi_hotspot.features import NEIGHBOR_FEATURES
from topoppi_hotspot.footprint_diffusion import (
    TIMES,
    add_graph_diffusion,
    build_diffusion_operators,
)


def _fixture():
    rows = []
    # A has unequal masses, an ordinary isolate, a zero-footprint/full-area
    # fallback, a wholly zero-area row, and an excluded noninterface row.
    definitions = [
        ("a", "r0", 0, 1., 1., 0., 3., "E"),
        ("a", "r1", 0, 3., 3., 4., 3., "A"),
        ("a", "r2", 0, 2., 2., 8., 3., "K"),
        ("a", "r3", 1, 2., 2., 100., 3., "E"),
        ("a", "r4", 1, 0., 4., 104., 3., "A"),
        ("a", "r5", 0, 0., 0., 1000., 3., "K"),
        ("a", "r6", 0, 1., 1., 10000., 8., "E"),
        ("b", "r0", 0, 1., 1., 50., 3., "E"),
        ("b", "r1", 0, 1., 1., 60., 3., "K"),
    ]
    for complex_id, key, side, footprint, area, value, distance, wt in definitions:
        row = {channel: value + number for number, channel in enumerate(NEIGHBOR_FEATURES)}
        row.update(complex_id=complex_id, residue_key=key, side=side, wt=wt,
                   **{"surface.footprint_area": footprint, "surface.area": area,
                      "contact.partner_distance": distance})
        rows.append(row)
    frame = pd.DataFrame(rows)
    # Edges crossing partners or leaving the interface must not transport data.
    keys_a = np.array(["r4", "r1", "r6", "r0", "r5", "r3", "r2"])
    pairs = [("r0", "r1"), ("r3", "r4"), ("r1", "r3"),
             ("r0", "r5"), ("r1", "r6")]
    indices = {key: i for i, key in enumerate(keys_a)}
    edges_a = np.array([(indices[a], indices[b]) for a, b in pairs])
    graphs = {
        "a": {"residue_keys": keys_a, "footprint_edges": edges_a,
              "footprint_weights": np.array([2., 5., 100., 10., 100.]),
              "spatial_edges": edges_a},
        "b": {"residue_keys": np.array(["r1", "r0"]),
              "footprint_edges": np.array([[0, 1]]), "footprint_weights": np.array([3.]),
              "spatial_edges": np.array([[0, 1]])},
    }
    return frame, graphs


def test_unequal_area_heat_matches_analytic_two_node_solution():
    frame, graphs = _fixture()
    operators = build_diffusion_operators(frame, graphs)
    values = frame.iloc[operators.source_rows]["chemical.hydropathy"].to_numpy()
    first = operators.queries[0]
    # W=2, M=(1,3), rate=mean(2,2/3)=4/3. The nonzero
    # normalized eigenvalue is 2 and the area-weighted equilibrium is 3.
    assert first.rate == pytest.approx(4 / 3)
    for time in TIMES:
        output = operators.apply(values, time)
        decay = np.exp(-2 * time)
        np.testing.assert_allclose(output[:2], [3 - 3 * decay, 3 + decay], atol=1e-12)
    np.testing.assert_array_equal(operators.apply(values, 0), values)


@pytest.mark.parametrize("kind", ["native", "spatial", "native_unweighted"])
def test_heat_preserves_constants_area_mass_and_isolated_rows(kind):
    frame, graphs = _fixture()
    operators = build_diffusion_operators(frame, graphs, kind)
    values = np.random.default_rng(3).normal(size=(len(operators.row_keys), 5))
    for time in TIMES:
        result = operators.apply(values, time)
        np.testing.assert_allclose(operators.apply(np.ones(len(values)), time), 1., atol=1e-12)
        for query in operators.queries:
            rows = query.row_indices
            np.testing.assert_allclose(query.mass @ result[rows], query.mass @ values[rows], atol=1e-11)
            np.testing.assert_array_equal(result[rows][query.isolated], values[rows][query.isolated])
    query_a0, query_a1 = operators.queries[:2]
    assert query_a0.metadata["n_zero_area_identity"] == 1
    assert query_a0.metadata["n_incident_edges_removed_for_zero_area"] == 1
    assert query_a1.metadata["n_full_surface_mass_fallback"] == 1
    np.testing.assert_array_equal(query_a1.mass, [2., 4.])


def test_no_partner_complex_or_noninterface_mixing():
    frame, graphs = _fixture()
    operators = build_diffusion_operators(frame, graphs)
    values = frame.iloc[operators.source_rows]["chemical.hydropathy"].to_numpy()
    baseline = operators.apply(values, 1)
    altered = values.copy()
    other = (operators.row_keys.side == 1) | (operators.row_keys.complex_id == "b")
    altered[other] += 1e6
    np.testing.assert_array_equal(operators.apply(altered, 1)[~other], baseline[~other])
    assert not (operators.row_keys.residue_key == "r6").any()
    assert operators.metadata["n_interface_rows"] == 8
    # Native alanine is part of the context, not dropped as an unmutated site.
    assert (operators.row_keys.residue_key == "r1").sum() == 2


def test_row_mapping_is_positional_and_independent_of_input_and_graph_order():
    frame, graphs = _fixture()
    shuffled = frame.sample(frac=1, random_state=9)
    shuffled.index = np.zeros(len(shuffled), dtype=int)  # duplicate dataframe labels are immaterial
    operators = build_diffusion_operators(shuffled, graphs)
    expected = frame.loc[frame["contact.partner_distance"] <= 6].sort_values(["complex_id", "residue_key"])
    pd.testing.assert_frame_equal(operators.row_keys, expected[["complex_id", "residue_key", "side"]].reset_index(drop=True))
    mapped = shuffled.iloc[operators.source_rows]
    np.testing.assert_array_equal(mapped.residue_key, operators.row_keys.residue_key)
    np.testing.assert_array_equal(mapped.complex_id, operators.row_keys.complex_id)
    reference = build_diffusion_operators(frame, graphs)
    x = mapped[NEIGHBOR_FEATURES].to_numpy()
    np.testing.assert_allclose(operators.apply(x, 1), reference.apply(expected[NEIGHBOR_FEATURES], 1), atol=1e-12)


def test_uniform_boundary_scaling_cancels_under_declared_rate_normalization():
    frame, graphs = _fixture()
    changed = {name: dict(graph) for name, graph in graphs.items()}
    for graph in changed.values():
        graph["footprint_weights"] = graph["footprint_weights"] * 7
    before = build_diffusion_operators(frame, graphs)
    after = build_diffusion_operators(frame, changed)
    x = frame.iloc[before.source_rows][NEIGHBOR_FEATURES].to_numpy()
    np.testing.assert_allclose(before.apply(x, 1), after.apply(x, 1), atol=1e-11)
    assert after.queries[0].rate == pytest.approx(7 * before.queries[0].rate)


def test_feature_blocks_preserve_input_and_name_heat_and_modulus_separately():
    frame, graphs = _fixture()
    original = frame.copy(deep=True)
    result, blocks, metadata = add_graph_diffusion(frame, graphs)
    pd.testing.assert_frame_equal(frame, original)
    expected = original.loc[original["contact.partner_distance"] <= 6].sort_values(["complex_id", "residue_key"]).reset_index(drop=True)
    pd.testing.assert_frame_equal(result[original.columns], expected)
    assert {name: len(columns) for name, columns in blocks.items()} == {
        "native_heat": 48, "native_scattering": 32, "spatial_heat": 48, "spatial_scattering": 32,
    }
    column = "scattering.native.t0p25_t1.chemical.hydropathy"
    # Analytic first-node modulus contrast, not a signed linear residual.
    assert result.loc[0, column] == pytest.approx(3 * (np.exp(-.5) - np.exp(-2)))
    assert not metadata["labels_read"]
    assert not metadata["operators"]["native"]["mesh_laplace_beltrami"]
    assert metadata["operators"]["spatial"]["edge_weight"] == "unit"
    changed = frame.copy()
    changed["hotspot"] = np.arange(len(frame)) % 2
    again, _, _ = add_graph_diffusion(changed, graphs)
    all_new = sum(blocks.values(), [])
    np.testing.assert_array_equal(again[all_new], result[all_new])


def test_zero_area_only_query_is_an_identity_without_invented_mass():
    frame, graphs = _fixture()
    frame = frame.loc[(frame.complex_id == "a") & (frame.residue_key == "r5")]
    operators = build_diffusion_operators(frame, graphs)
    np.testing.assert_array_equal(operators.apply([[5., 9.]], 4), [[5., 9.]])
    np.testing.assert_array_equal(operators.queries[0].mass, [0.])
    assert operators.queries[0].rate == 0


def test_absent_residue_mapping_fails_instead_of_using_wrong_graph_row():
    frame, graphs = _fixture()
    frame.loc[0, "residue_key"] = "missing"
    with pytest.raises(ValueError, match="map to a graph residue key"):
        build_diffusion_operators(frame, graphs)
