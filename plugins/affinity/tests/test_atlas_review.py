"""Scientific correspondence failures and portable review state."""

from types import SimpleNamespace
from unittest.mock import Mock

import pandas as pd
import pytest


@pytest.mark.parametrize("prediction", [float("nan"), float("inf"), ""])
def test_failed_prediction_review_uses_readable_labels_and_preserves_data(prediction):
    from topoppi_affinity.atlas_review.gui import ReviewWindow

    predictions = pd.DataFrame([{"variant_id": "failed", "predicted": prediction, "sites": "[]",
                                 "mutations": "A:G1A", "status": "error", "target": "ddg"}])
    original = predictions.copy(deep=True)
    session = SimpleNamespace(document=SimpleNamespace(predictions=predictions, report={}),
                              variant_id="failed", selected=set(), select_variant=lambda _: None)
    output = Mock()
    window = SimpleNamespace(session=session, variant=SimpleNamespace(get=lambda: "failed"),
                             choose=lambda _: None, output=output)
    ReviewWindow.variant_selected(window, None)
    assert "Prediction: Unavailable kcal/mol" in output.set.call_args.args[0]
    assert "Failed" in output.set.call_args.args[0]
    pd.testing.assert_frame_equal(predictions, original)












def test_contact_comparison_does_not_turn_unmapped_sites_into_losses():
    from topoppi_affinity.atlas_review.comparison import contact_delta

    mapping = [dict(reference_key="A:GLY:1", mobile_key="X:ALA:8"), dict(reference_key="B:SER:2", mobile_key="Y:SER:3")]
    left = [
        dict(source="A:GLY:1", target="B:SER:2", min_distance=5.9),
        dict(source="A:GLY:1", target="B:GLY:9", min_distance=4.0),
    ]
    right = [dict(source="X:ALA:8", target="Y:SER:3", min_distance=6.1)]
    assert len(contact_delta(mapping, left, right)) == 1
    assert contact_delta(mapping, left, right)[0]["status"] == "reference_only"
    assert contact_delta(mapping, left, right, 6.5)[0]["status"] == "shared"



