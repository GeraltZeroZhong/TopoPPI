"""Cancelled work remains distinguishable from failed or completed results."""
from __future__ import annotations

import threading
from pathlib import Path
from types import SimpleNamespace

import joblib
import pandas as pd
import pytest
from sklearn.dummy import DummyRegressor

from topoppi_affinity import cli, workflow
from topoppi_affinity.models import AffinityRegressor

FIXTURE = Path(__file__).parent / "fixtures/tiny_complex.pdb"


@pytest.fixture
def pilot_model(tmp_path):
    model = AffinityRegressor(DummyRegressor(strategy="constant", constant=1.5).fit([[1]], [1.5]),
        ("contact.residue_pairs",), "ddg", "test constant", {"scope": "test only"})
    path = tmp_path / "pilot.joblib"
    joblib.dump(model, path)
    return path


@pytest.mark.parametrize("stage", ["interface", "prediction"])
def test_cancelled_rows_keep_validation_failures_distinct(tmp_path, pilot_model, monkeypatch, stage):
    cancel = threading.Event()

    def stop(*args, **kwargs):
        cancel.set()
        raise InterruptedError("Analysis cancelled")

    target = ("topoppi_affinity.workflow.prepare_interface" if stage == "interface"
              else "topoppi_affinity.inference.build_features")
    monkeypatch.setattr(target, stop)
    output = tmp_path / "result"
    report = workflow.analyze(FIXTURE, "A", "B", pilot_model, output,
        mutations=["GA1W", "GA99W"], cancel_event=cancel, log=lambda _: None)
    assert report["prediction_status"] == "cancelled"
    predictions = workflow.load_analysis(output).predictions
    assert predictions.status.tolist() == ["cancelled", "error"]
    assert "cancelled" in predictions.error.iloc[0]
    assert "absent" in predictions.error.iloc[1]


def test_cancelled_second_atlas_keeps_first_atlas_and_predictions(tmp_path, pilot_model, monkeypatch):
    cancel = threading.Event()

    def layout(prepared, *, side, **kwargs):
        if side == 1:
            cancel.set()
            raise InterruptedError("Analysis cancelled")
        return prepared.partners[side].patches, {}

    monkeypatch.setattr("topoppi.pipeline.layout_interface", layout)
    monkeypatch.setattr("topoppi.visualization.atlas_io.save_atlas",
                        lambda path, *args, **kwargs: Path(path).write_bytes(b"completed atlas"))
    output = tmp_path / "result"
    report = workflow.analyze(FIXTURE, "A", "B", pilot_model, output, mutations=["GA1W"],
        layout=True, cancel_event=cancel, log=lambda _: None)
    assert report["prediction_status"] == "ok"
    assert report["layout_status"] == "partial"
    assert report["layout"]["0"]["status"] == "ok"
    assert report["layout"]["1"]["status"] == "cancelled"
    assert (output / "atlas_0.npz").read_bytes() == b"completed atlas"
    assert workflow.load_analysis(output).predictions.status.tolist() == ["ok"]


@pytest.mark.parametrize("command", ["analyze", "layout", "example"])
def test_cli_cancellation_returns_130_with_available_partial_results(tmp_path, monkeypatch, capsys, command):
    report = {"output_directory": str(tmp_path), "prediction_status": "ok", "layout_status": "partial"}
    document = SimpleNamespace(path=tmp_path / "analysis.json", report=report,
        predictions=pd.DataFrame([{"variant_id": "001", "mutations": "", "predicted": 7.,
                                   "target": "pkd", "status": "ok", "error": ""}]))

    def stop(*args, cancel_event, **kwargs):
        cancel_event.set()
        return report

    monkeypatch.setattr(workflow, "analyze", stop)
    monkeypatch.setattr(workflow, "layout_existing", stop)
    monkeypatch.setattr(workflow, "load_analysis", lambda path: document)
    arguments = {
        "analyze": ["analyze", "--structure", str(FIXTURE), "--chains-a", "A", "--chains-b", "B"],
        "layout": ["layout", str(tmp_path)],
        "example": ["example", "--output", str(tmp_path / "example"), "--run"],
    }
    assert cli.main(arguments[command]) == 130
    assert "Prediction: ok | Layout: partial" in capsys.readouterr().out
    assert report["layout_status"] == "partial"


def test_download_progress_throttles_actual_downloading_events(monkeypatch):
    times = iter([0., .2, .9, 1.1, 1.2])
    messages = []
    monkeypatch.setattr(cli.time, "monotonic", lambda: next(times))
    monkeypatch.setattr(cli, "_log", messages.append)
    progress = cli._DownloadProgress()
    for number in range(4):
        progress({"event": "downloading", "resource": "ESMC", "message": f"chunk {number}"})
    progress({"event": "downloaded", "resource": "ESMC", "message": "complete"})
    assert messages == ["chunk 0", "chunk 3", "complete"]


@pytest.mark.parametrize("failure", ["solver", "save", "cancelled"])
def test_regeneration_failure_reports_attempt_and_keeps_viewable_old_atlases(
        tmp_path, pilot_model, monkeypatch, capsys, failure):
    output = tmp_path / "result"
    report = workflow.analyze(FIXTURE, "A", "B", pilot_model, output, mutations=["GA1W"], log=lambda _: None)
    report["layout_status"] = "ok"
    report["layout"] = {str(side): {"status": "ok", "atlas": f"atlas_{side}.npz"} for side in (0, 1)}
    report["layout_attempt_warning"] = "Previous layout warning"
    report["warnings"] = ["Model applicability note", "Previous layout warning"]
    for side in (0, 1):
        (output / f"atlas_{side}.npz").write_bytes(f"original atlas {side}".encode())
    workflow._write_report(report, output)

    def layout(prepared, *, side, **kwargs):
        if failure == "solver":
            raise RuntimeError("New layout parameters failed")
        if failure == "cancelled":
            raise InterruptedError("Atlas regeneration cancelled")
        return prepared.partners[side].patches, {}

    def save(path, *args, **kwargs):
        Path(path).write_bytes(b"incomplete new archive")
        raise OSError("No space left for new archive")

    monkeypatch.setattr("topoppi.pipeline.layout_interface", layout)
    monkeypatch.setattr("topoppi.visualization.atlas_io.save_atlas", save)
    assert cli.main(["layout", str(output)]) == (130 if failure == "cancelled" else 2)
    text = capsys.readouterr().out
    attempt = "cancelled" if failure == "cancelled" else "error"
    assert f"Latest atlas attempt: {attempt}" in text
    assert "Earlier atlases for A, B remain viewable" in text
    assert "Previous layout warning" not in text
    saved = workflow.load_analysis(output)
    assert saved.report["layout_status"] == "ok"
    assert saved.report["layout_attempt_status"] == attempt
    assert "Model applicability note" in saved.report["warnings"]
    assert saved.predictions.status.tolist() == ["ok"]
    for side in (0, 1):
        assert saved.report["layout"][str(side)]["status"] == "ok"
        assert saved.report["layout"][str(side)]["last_attempt_error"]
        assert (output / f"atlas_{side}.npz").read_bytes() == f"original atlas {side}".encode()
    assert not list(output.glob(".atlas_*.pending.npz"))

    if failure == "solver":
        monkeypatch.setattr("topoppi.pipeline.layout_interface",
                            lambda prepared, *, side, **kwargs: (prepared.partners[side].patches, {}))
        monkeypatch.setattr("topoppi.visualization.atlas_io.save_atlas",
                            lambda path, *args, **kwargs: Path(path).write_bytes(b"complete new archive"))
        assert cli.main(["layout", str(output)]) == 0
        recovered = workflow.load_analysis(output)
        assert recovered.report["layout_attempt_status"] == "ok"
        assert recovered.report["warnings"] == ["Model applicability note"]
        assert "layout_attempt_warning" not in recovered.report
        assert all("last_attempt_error" not in entry for entry in recovered.report["layout"].values())
