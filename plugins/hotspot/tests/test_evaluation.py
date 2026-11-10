import numpy as np
import pandas as pd
import pytest
from topoppi_hotspot.evaluation import bootstrap_comparisons, evaluate, metrics


def _cohort():
    rows = []
    for group in range(4):
        for residue in range(8):
            energy = residue * .5 + (group - 1.5) * .2
            rows.append({"family_group": f"g{group}", "complex_id": f"c{group}",
                         "residue_key": f"A:TYR:{residue}", "Er": energy,
                         "hotspot": int(energy >= 2), "x": residue + group / 10})
    return pd.DataFrame(rows)


def test_held_out_labels_cannot_change_own_predictions():
    original = _cohort()
    before, _, _, _, folds = evaluate(original, {"base": ["x"]})
    changed = original.copy()
    mask = changed.family_group == "g0"
    changed.loc[mask, "Er"] += 80
    changed.loc[mask, "hotspot"] = 1 - changed.loc[mask, "hotspot"]
    after, _, _, _, _ = evaluate(changed, {"base": ["x"]})
    cols = ["P_hotspot", "predicted_Er"]
    np.testing.assert_allclose(before.loc[before.family_group == "g0", cols],
                               after.loc[after.family_group == "g0", cols])
    assert before.P_hotspot.between(0, 1).all()
    assert all(f["held_out_group"] not in f["train_groups"] for f in folds)
    assert not before.duplicated(["model", "complex_id", "residue_key"]).any()


def test_paired_bootstrap_rejects_incomparable_residues():
    frame = _cohort().assign(P_hotspot=.4, predicted_Er=1., model="one")
    second = frame.iloc[1:].assign(model="two")
    with pytest.raises(ValueError, match="different evaluation residues"):
        bootstrap_comparisons(pd.concat([frame, second]), [("one", "two")], n_bootstrap=10)


def test_constant_probabilities_have_no_artificial_top_three_order():
    result = metrics([1, 1, 1, 0, 0, 0], [3, 3, 3, 0, 0, 0], [.5] * 6, [1.5] * 6)
    assert result["precision_at_3"] == .5
