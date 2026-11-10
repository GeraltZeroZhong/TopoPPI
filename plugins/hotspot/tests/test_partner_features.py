"""Independent synthetic checks for reusable numerical primitives."""
import numpy as np
import pandas as pd
from topoppi_hotspot.features import FEATURE_SETS, NEIGHBOR_FEATURES
from topoppi_hotspot.partner_features import feature_specifications, prepare_partner_features


def _fixture():
    rows, labels, graphs = [], [], {}
    columns = sorted(set(FEATURE_SETS["spatial"] + FEATURE_SETS["topology"]))
    for family in range(3):
        complex_id = f"complex{family}"
        for i in range(6):
            values = {name: float((i + family + number) % 7) for number, name in enumerate(columns)}
            values["contact.partner_distance"] = 3.
            values["surface.area"] = float(i + 10 * (i // 3))
            rows.append({"complex_id": complex_id, "residue_key": f"R{i}", "side": i // 3,
                         "wt": "A" if i % 3 == 2 else "E", **values})
            if i % 3 != 2:
                labels.append({"complex_id": complex_id, "residue_key": f"R{i}",
                               "family_group": f"family{family}", "hotspot": i % 3})
        graphs[complex_id] = {"residue_keys": np.asarray([f"R{i}" for i in range(6)]),
                              "spatial_edges": np.asarray([[0, 1], [0, 2], [3, 4], [3, 5]]),
                              "native_neighbor_edges": np.asarray([[0, 1], [0, 2], [3, 4], [3, 5]])}
    return pd.DataFrame(rows), pd.DataFrame(labels), graphs

def test_relative_population_includes_unmeasured_alanine_and_excludes_other_partner():
    frame, _, _ = _fixture()
    prepared = prepare_partner_features(frame)
    chosen = (prepared.complex_id == "complex0") & (prepared.side == 0)
    column = "relative.partner.percentile.surface.area"
    np.testing.assert_allclose(prepared.loc[chosen, column], [1 / 6, .5, 5 / 6])
    frame.loc[frame.side == 1, "surface.area"] *= 1000
    changed = prepare_partner_features(frame)
    np.testing.assert_allclose(prepared.loc[chosen, feature_specifications()["relative_self"]],
                               changed.loc[chosen, feature_specifications()["relative_self"]])

def test_neighbor_shuffle_uses_same_partner_source_and_is_paired_across_graphs():
    frame, _, graphs = _fixture()
    frame[NEIGHBOR_FEATURES] = np.repeat((frame.side * 100 + frame.index % 3).to_numpy()[:, None], len(NEIGHBOR_FEATURES), axis=1)
    prepared = prepare_partner_features(frame, graphs=graphs, shuffle_seed=2)
    assert prepared.loc[prepared.side == 0, "spatial.mean.chemical.volume"].max() <= 2
    assert prepared.loc[prepared.side == 1, "spatial.mean.chemical.volume"].min() >= 100
    np.testing.assert_allclose(prepared["spatial.mean.chemical.volume"], prepared["topology.mean.chemical.volume"])
    again = prepare_partner_features(frame.sample(frac=1, random_state=9), graphs=graphs, shuffle_seed=2)
    pd.testing.assert_frame_equal(prepared, again)
