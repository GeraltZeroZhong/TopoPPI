"""Check structural normalization populations and ranking feature provenance."""
import numpy as np
import pandas as pd
import pytest
from topoppi_hotspot.features import FEATURE_SETS
from topoppi_hotspot.relative_features import SELF_DESCRIPTORS, add_relative_features


@pytest.fixture
def all_features():
    # I1 contains two measured and two unmeasured interface residues, plus one
    # noninterface residue. I2 has one singleton partner and one constant pair.
    frame = pd.DataFrame({
        "complex_id": ["I1"] * 5 + ["I2"] * 3,
        "residue_key": ["A:1", "A:2", "B:1", "B:2", "B:3", "A:1", "B:1", "B:2"],
        "side": [0, 0, 1, 1, 1, 0, 1, 1],
        "measured": [True, False, True, False, False, True, False, True],
    })
    for source in sorted({c for columns in FEATURE_SETS.values() for c in columns}):
        frame[source] = [2., 4., 10., 16., 1000., 3., 9., 9.]
    frame["contact.partner_distance"] = [2., 3., 4., 6., 7., 4., 4., 4.]
    frame["chemical.aromatic"] = 1.
    return frame


def test_all_structural_interface_residues_define_population(all_features):
    frame, _, meta = add_relative_features(all_features)
    first = frame[frame.complex_id == "I1"]
    assert len(first) == 4
    assert meta["n_interface_residues"] == 7
    assert first["context.interface.mean.surface.area"].tolist() == [8.] * 4
    np.testing.assert_allclose(first["relative.interface.percentile.surface.area"],
                               [.125, .375, .625, .875])
    measured = first[first.measured]
    # Measured-only normalization would produce mean 6 and percentiles .25,.75.
    np.testing.assert_allclose(measured["relative.interface.percentile.surface.area"], [.125, .625])
    assert meta["interface_populations"][0]["n_by_partner"] == {"0": 2, "1": 2}


def test_partner_normalization_is_separate_and_uses_population_std(all_features):
    frame, _, _ = add_relative_features(all_features)
    np.testing.assert_allclose(frame.iloc[:4]["relative.partner.zscore.surface.area"], [-1, 1, -1, 1])
    np.testing.assert_allclose(frame.iloc[:4]["relative.partner.percentile.surface.area"], [.25, .75, .25, .75])
    changed = all_features.copy()
    changed.loc[(changed.complex_id == "I1") & (changed.side == 1), "surface.area"] += 100.
    other, _, _ = add_relative_features(changed)
    columns = [c for c in frame if c.startswith("relative.partner.")]
    pd.testing.assert_frame_equal(frame.loc[frame.side == 0, columns], other.loc[other.side == 0, columns])


def test_constant_and_singleton_strata_have_zero_zscore_and_mid_percentile(all_features):
    frame, _, _ = add_relative_features(all_features)
    assert (frame["relative.interface.zscore.chemical.aromatic"] == 0).all()
    assert (frame["relative.interface.percentile.chemical.aromatic"] == .5).all()
    second = frame[frame.complex_id == "I2"]
    assert (second["relative.partner.zscore.surface.area"] == 0).all()
    assert (second["relative.partner.percentile.surface.area"] == .5).all()
    # Interface-level tied ranks are averaged, independently of partner strata.
    np.testing.assert_allclose(second["relative.interface.percentile.surface.area"], [1 / 6, 2 / 3, 2 / 3])


def test_labels_and_input_order_do_not_change_relative_features(all_features):
    frame, feature_sets, metadata = add_relative_features(all_features)
    altered = all_features.sample(frac=1., random_state=19).copy()
    altered["hotspot"] = np.arange(len(altered)) % 2
    altered["Er"] = np.arange(len(altered)) * 100.
    altered["measured"] = ~altered.measured
    other, other_sets, other_meta = add_relative_features(altered)
    columns = [c for c in frame if c.startswith(("context.", "relative."))]
    pd.testing.assert_frame_equal(frame[columns].sort_index(), other[columns].sort_index())
    assert feature_sets == other_sets
    assert metadata == other_meta


def test_noninterface_and_other_interface_values_do_not_change_normalization(all_features):
    frame, _, _ = add_relative_features(all_features)
    altered = all_features.copy()
    altered.loc[altered["contact.partner_distance"] > 6, "surface.area"] = -1e8
    altered.loc[altered.complex_id == "I2", "surface.area"] = 1234.
    other, _, _ = add_relative_features(altered)
    pd.testing.assert_frame_equal(frame[frame.complex_id == "I1"], other[other.complex_id == "I1"])


def test_decomposition_and_feature_sets_preserve_ablation_boundaries(all_features):
    frame, sets, metadata = add_relative_features(all_features)
    original = all_features[all_features["contact.partner_distance"] <= 6]
    pd.testing.assert_frame_equal(frame[original.columns], original)
    for source, mapping in metadata["feature_columns"].items():
        np.testing.assert_allclose(frame[source], frame[mapping["context"]] + frame[mapping["deviation"]])
    for kind in ("self", "spatial", "topology", "shuffled"):
        columns = sets[f"relative_{kind}"]
        assert columns[:len(FEATURE_SETS[kind])] == FEATURE_SETS[kind]
        assert not any(column.startswith("context.") for column in columns)
        assert all(column in frame for column in columns)
        assert len(columns) == len(set(columns))
    assert len(sets["relative_self"]) == len(FEATURE_SETS["self"]) + 5 * len(SELF_DESCRIPTORS)
    assert len(sets["relative_topology"]) == len(FEATURE_SETS["topology"]) + 5 * (len(SELF_DESCRIPTORS) + 5)
    assert not any("topology" in column or "spatial" in column or "shuffled" in column
                   for column in sets["relative_self"])
    assert not any("shuffled" in column or "spatial" in column for column in sets["relative_topology"])
