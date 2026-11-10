"""Observed-site nulls must not turn unmeasured interface sites into negatives."""
import numpy as np
import pandas as pd
import pytest
from topoppi_hotspot.organization import _bh, analyze_organization


def fixture():
    n = 10
    graph = {
        "residue_keys": np.asarray([f"A:{i}" for i in range(n)]),
        "side": np.zeros(n, dtype=int),
        "positions": np.column_stack(([0., 1., 2., 30., 45., 60., 75., 90., .5, 100.], np.zeros((n, 2)))),
        "ca_resolved": np.ones(n, dtype=bool),
        "interface_mask": np.ones(n, dtype=bool),
        "footprint_area": np.arange(1., n + 1),
        "contact_degree": np.asarray([8., 7., 6., 1., 1., 1., 1., 1., 100., 100.]),
        "contact_betweenness": np.asarray([.8, .7, .6, .1, .1, .1, .1, .1, 1., 1.]),
        "footprint_edges": np.asarray([[0, 1], [1, 2], [2, 8], [8, 3]]),
        "footprint_weights": np.ones(4),
    }
    labels = pd.DataFrame({"complex_id": ["example"] * 8, "residue_key": graph["residue_keys"][:8],
                           "Er": [3., 3., 3., 0., 0., 0., 0., 0.],
                           "hotspot": [True, True, True, False, False, False, False, False]})
    return labels, {"example": graph}


def test_observed_metrics_coverage_and_fixed_seed_null():
    labels, graphs = fixture()
    rows, summary = analyze_organization(labels, graphs, n_permutations=99)
    repeat, _ = analyze_organization(labels.sample(frac=1., random_state=7), graphs, n_permutations=99)
    pd.testing.assert_frame_equal(rows, repeat)
    values = rows.set_index("metric")["value"]
    assert values["hotspot_nearest_ca_distance_angstrom"] == pytest.approx(1.)
    assert values["hotspot_ca10_edge_count"] == 3
    assert values["hotspot_footprint_island_count"] == 1
    assert values["hotspot_footprint_largest_island_fraction"] == 1
    assert values["hotspot_contact_degree_difference"] == pytest.approx(6.)
    assert values["hotspot_contact_betweenness_difference"] == pytest.approx(.6)
    assert values["hotspot_footprint_area_share"] == pytest.approx(6. / 36.)
    assert values["positive_er_hhi"] == pytest.approx(1. / 3.)
    assert values["positive_er_footprint_area_concentration"] == pytest.approx(36. / 9. * (1. + .5 + 1. / 3.))
    assert summary["coverage"][0]["n_interface"] == 10
    assert summary["coverage"][0]["n_measured"] == 8
    assert summary["coverage"][0]["measured_footprint_area_fraction"] == pytest.approx(36. / 55.)
    hhi = rows.set_index("metric").loc["positive_er_hhi"]
    assert np.isnan(hhi.p_value)
    assert hhi.reason == "relabeling_invariant_descriptive_metric"
    tested = rows[rows.test_status == "tested"]
    assert len(tested) > 0
    assert (tested.p_value >= .01).all()
    assert (tested.q_value_bh >= tested.p_value).all()


def test_unmeasured_nodes_never_bridge_hotspot_islands_or_become_negatives():
    labels, graphs = fixture()
    labels.loc[:, "hotspot"] = [True, False, False, True, False, False, False, False]
    labels.loc[:, "Er"] = [3., 0., 0., 3., 0., 0., 0., 0.]
    graphs["example"]["footprint_edges"] = np.asarray([[0, 8], [8, 3]])
    graphs["example"]["footprint_weights"] = np.ones(2)
    rows, summary = analyze_organization(labels, graphs, n_permutations=19)
    values = rows.set_index("metric")["value"]
    assert values["hotspot_footprint_island_count"] == 2
    assert values["hotspot_footprint_largest_island_fraction"] == .5
    assert set(rows.n_eligible) == {8}
    assert summary["coverage"][0]["n_hotspot"] == 2


@pytest.mark.parametrize("hotspot,reason", [
    ([False] * 8, "no_measured_hotspots"),
    ([True] * 8, "all_measured_sites_are_hotspots"),
    ([True] + [False] * 7, "fewer_than_two_measured_hotspots"),
    ([True] * 7 + [False], "fewer_than_two_measured_nonhotspots"),
])
def test_degenerate_label_sets_have_no_pseudotests(hotspot, reason):
    labels, graphs = fixture()
    labels["hotspot"] = hotspot
    rows, summary = analyze_organization(labels, graphs, n_permutations=9)
    assert rows.p_value.isna().all()
    assert rows.q_value_bh.isna().all()
    assert summary["n_tests"] == 0
    assert set(rows.loc[rows.alternative != "descriptive", "reason"]) == {reason}


def test_interface_ca_and_footprint_eligibility_are_explicit():
    labels, graphs = fixture()
    graph = graphs["example"]
    graph["interface_mask"][7] = False
    graph["ca_resolved"][0] = False
    graph["footprint_area"][1] = 0.
    rows, summary = analyze_organization(labels, graphs, n_permutations=19)
    by_metric = rows.set_index("metric")
    assert set(rows.n_measured) == {7}
    assert by_metric.loc["hotspot_nearest_ca_distance_angstrom", "n_eligible"] == 6
    assert by_metric.loc["hotspot_footprint_island_count", "n_eligible"] == 6
    assert by_metric.loc["hotspot_contact_degree_difference", "n_eligible"] == 7
    coverage = summary["coverage"][0]
    assert coverage["n_measured_noninterface_excluded"] == 1
    assert coverage["n_measured_missing_ca"] == 1
    assert coverage["n_measured_zero_footprint"] == 1


def test_two_partner_sides_are_not_pooled():
    labels, graphs = fixture()
    graphs["example"]["side"][4:] = 1
    rows, summary = analyze_organization(labels, graphs, n_permutations=19)
    assert set(rows.side) == {"0", "1"}
    assert len(summary["coverage"]) == 2
    assert rows.p_value.isna().all()
    side_one = rows[rows.side == "1"]
    assert set(side_one.n_measured) == {4}
    assert set(side_one.n_hotspot) == {0}


def test_missing_structure_and_mapping_are_reported_and_duplicates_rejected():
    labels, graphs = fixture()
    extra = pd.DataFrame([dict(complex_id="absent", residue_key="Z:99", Er=3., hotspot=True),
                          dict(complex_id="example", residue_key="Z:99", Er=3., hotspot=True)])
    _, summary = analyze_organization(pd.concat([labels, extra]), graphs, n_permutations=9)
    assert summary["missing_structure_complexes"] == ["absent"]
    assert summary["unmapped_measured_residues"] == [{"complex_id": "example", "residue_key": "Z:99"}]
    with pytest.raises(ValueError, match="duplicate residue labels"):
        analyze_organization(pd.concat([labels, labels.iloc[[0]]]), graphs)


def test_bh_uses_only_testable_comparisons_and_is_monotonic():
    values = _bh([.01, .04, np.nan, .03])
    np.testing.assert_allclose(values, [.03, .04, np.nan, .04], equal_nan=True)
