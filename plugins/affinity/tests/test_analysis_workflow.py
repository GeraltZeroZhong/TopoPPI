"""Visible workflow failures: loss of predictions, identity, joint effects or review."""
import json
import os
from pathlib import Path

import joblib
import numpy as np
import pytest
from sklearn.dummy import DummyRegressor

from topoppi_affinity.models import AffinityRegressor
from topoppi_affinity.workflow import analyze, load_analysis

FIXTURE = Path(__file__).parent / "fixtures/tiny_complex.pdb"


@pytest.fixture
def model(tmp_path):
    regressor = AffinityRegressor(DummyRegressor(strategy="constant", constant=1.5).fit([[1]], [1.5]),
        ("contact.residue_pairs",), "ddg", "test constant", {"scope": "test only"})
    path = tmp_path / "model.joblib"
    joblib.dump(regressor, path)
    return path


@pytest.mark.parametrize("failure", ["missing", "timeout", "solver_error"])
def test_layout_failures_preserve_predictions_and_reopen_without_models(tmp_path, model, monkeypatch, failure):
    if failure != "missing":
        def fail(*args, **kwargs):
            config = kwargs["optcuts"]
            assert config.residue_fragmentation_weight == 20.0
            assert config.optcuts_mode == config.optcuts_headless_mode
            raise TimeoutError("test timeout") if failure == "timeout" else RuntimeError("solver failed")
        monkeypatch.setattr("topoppi.pipeline.layout_interface", fail)
    output = tmp_path / "analysis"
    report = analyze(FIXTURE, "A", "B", model, output,
        mutations=["GA1W", "GA1F", "GA1W,GB1A", "GA99W"], wt_kd=1e-8,
        layout=True, optcuts_bin=str(tmp_path / "absent"), log=lambda _: None)
    assert report["prediction_status"] == "partial"
    assert report["layout_status"] == "error"
    model.unlink()
    doc = load_analysis(output)
    predictions = doc.predictions
    np.testing.assert_allclose(predictions.predicted.iloc[:3].astype(float), 1.5)
    assert predictions.iloc[3].status == "error"
    assert "absent" in predictions.iloc[3].error
    assert float(predictions.iloc[0].anchored_mutant_kd_M) > 1e-8
    joint, left, hidden = doc.selection("variant_003", 0)
    assert len(left) == 1 and len(hidden) == 1
    assert len(json.loads(joint.sites)) == 2
    assert len(doc.selection("variant_003", 1)[1]) == 1
    assert len(predictions[predictions.mutations.str.startswith("GA1")]) == 3
    patches, _, style = doc.view(1, "surface", "variant_003", "contact_degree")
    assert patches and style["highlight_residues"] == ["B:GLY:1"]
    assert "predicted_ddg" not in doc.residues
    assert not any("uv" in p.metadata for p in patches)


def test_prediction_failure_retains_interface_and_report(tmp_path, model, monkeypatch):
    def fail(*args, **kwargs):
        raise RuntimeError("inference failure")
    monkeypatch.setattr("topoppi_affinity.inference.predict_features", fail)
    report = analyze(FIXTURE, "A", "B", model, tmp_path / "analysis", mutations=["GA1W"], log=lambda _: None)
    assert report["interface_status"] == "ok"
    assert report["prediction_status"] == "error"
    doc = load_analysis(tmp_path / "analysis")
    assert doc.view()[0]
    assert doc.predictions.predicted.eq("").all()
    with pytest.raises(ValueError, match="already contains"):
        analyze(FIXTURE, "A", "B", model, tmp_path / "analysis", mutations=["GA1W"])


@pytest.mark.skipif(not os.environ.get("DISPLAY"), reason="Tk display required")
def test_gui_reopens_joint_variants_changes_side_and_keeps_unmapped_rows(tmp_path, model):
    import tkinter as tk

    from topoppi_affinity.gui import AffinityPanel

    output = tmp_path / "analysis"
    analyze(FIXTURE, "A", "B", model, output, mutations=["GA1W,GB1A", "GA99W"], log=lambda _: None)
    root = tk.Tk()
    root.withdraw()
    panel = AffinityPanel(root, analysis=output)
    try:
        root.update()
        assert panel.figure is not None
        panel.table.selection_set("variant_001")
        panel.draw()
        assert "Joint variant total" in panel.details.get()
        panel.side.set("1")
        panel.draw()
        assert "Substitution sites on this partner: 1" in panel.details.get()
        panel.table.selection_set("variant_002")
        panel.draw()
        assert "absent from the structure" in panel.details.get()
        assert len(panel.table.get_children()) == 2
        panel.mode.set("atlas")
        panel.draw()
        assert "No atlas available" in panel.details.get()
        assert panel.figure is None  # Do not leave a stale opposite-side/3D image under the new view label.
    finally:
        panel.close()
        root.destroy()
