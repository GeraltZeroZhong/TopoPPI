import numpy as np

from topoppi_affinity.surface_transport import interface_measure, transport_features


def test_area_marginals_detect_mass_mismatch_missed_by_nearest_neighbor():
    # Both clouds have the same support: Chamfer is zero, but 70% of the mass
    # must move ten Angstroms under the balanced global matching constraint.
    graph = {
        "point_pos": np.array([[0., 0., 0.], [10., 0., 0.], [0., 0., 0.], [10., 0., 0.]]),
        "point_normal": np.array([[0., 0., 1.], [0., 0., 1.], [0., 0., -1.], [0., 0., -1.]]),
        "side": np.array([0, 0, 1, 1]), "point_node_index": np.arange(4),
        "point_interface_weight": np.array([.9, .1, .2, .8]),
    }
    features, metadata = transport_features(graph)
    assert features["geometry.area_chamfer"] == 0.
    np.testing.assert_allclose(features["geometry.exact_ot_cost"], 70.)
    np.testing.assert_allclose(features["geometry.sinkhorn_linear_cost"], 70., atol=1e-6)
    assert features["normal.under_geometry_plan"] == 0.
    assert metadata["geometry_solver"]["marginal_max_abs_error"] < 1e-6


def test_subsampling_preserves_total_native_footprint_area():
    graph = {"point_pos": np.arange(60).reshape(20, 3), "point_normal": np.tile([0., 0., 1.], (20, 1)),
             "side": np.zeros(20, dtype=int), "point_node_index": np.arange(20),
             "point_interface_weight": np.arange(1, 21, dtype=float)}
    measure = interface_measure(graph, 0, 8)
    assert measure["n_points"] <= 8
    assert measure["area"] == 210.
    np.testing.assert_allclose(measure["mass"].sum(), 1.)
