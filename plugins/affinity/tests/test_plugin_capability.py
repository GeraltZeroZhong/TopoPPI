"""The reference capability consumes existing Core geometry without rereading it."""
import shutil
from pathlib import Path

import joblib
import numpy as np
from sklearn.dummy import DummyRegressor
from topoppi import compute_geometry, get_interface, load_structure, plugins
from topoppi.execution import ExecutionContext

from topoppi_affinity.features import extract_complex_features
from topoppi_affinity.models import AffinityRegressor
from topoppi_affinity.workflow import analyze, load_analysis

FIXTURE = Path(__file__).parent / "fixtures/tiny_complex.pdb"


def test_object_capability_matches_path_predictions_after_source_removal(tmp_path, monkeypatch):
    source = tmp_path / "toy.pdb"
    shutil.copyfile(FIXTURE, source)
    geometry = compute_geometry(get_interface(load_structure(source), "A", "B"))
    bundle = extract_complex_features(source, "A", "B", prepared=geometry)
    regressor = AffinityRegressor(DummyRegressor(strategy="constant", constant=1.5).fit([[1]], [1.5]),
        ("contact.residue_pairs",), "ddg", "synthetic constant", {"n_train": 1})
    model = tmp_path / "toy.joblib"
    joblib.dump(regressor, model)
    expected = analyze(source, "A", "B", model, tmp_path / "path-result", mutations=["GA1W", "GA99W"],
                       log=lambda _: None)
    source.unlink()

    def unexpected(*args, **kwargs):
        raise AssertionError("The capability reopened the structure or recomputed geometry")

    monkeypatch.setattr("topoppi_affinity.workflow.PDBLoader", unexpected)
    monkeypatch.setattr("topoppi_affinity.workflow.prepare_interface", unexpected)
    monkeypatch.setattr("topoppi_affinity.features.prepare_interface", unexpected)
    observed = extract_complex_features(source, "A", "B", prepared=geometry)
    assert observed == bundle
    result = plugins.run("affinity.predict", geometry, model=model, mutations=["GA1W", "GA99W"],
                         context=ExecutionContext(output_dir=tmp_path / "object-result"))
    assert result["prediction_status"] == expected["prediction_status"] == "partial"
    before, after = [load_analysis(tmp_path / name).predictions for name in ("path-result", "object-result")]
    np.testing.assert_array_equal(before.predicted, after.predicted)
    assert before.status.tolist() == after.status.tolist()
    assert result["analysis_id"] == geometry.analysis_id


def test_graph_consumes_existing_task_descriptors(tmp_path, monkeypatch):
    from topoppi_affinity.graph_features import extract_complex_graph

    geometry = compute_geometry(get_interface(load_structure(FIXTURE), "A", "B"))
    descriptors = extract_complex_features(FIXTURE, "A", "B", prepared=geometry)
    reference = extract_complex_graph(FIXTURE, "A", "B", prepared=geometry)
    monkeypatch.setattr("topoppi_affinity.graph_features.extract_complex_features",
                        lambda *a, **kw: (_ for _ in ()).throw(AssertionError("Repeated descriptor calculation")))
    reused = extract_complex_graph(FIXTURE, "A", "B", prepared=geometry, descriptor_bundle=descriptors)
    for name in reference:
        np.testing.assert_array_equal(reference[name], reused[name])


def test_explicit_model_and_output_override_namespaced_defaults(monkeypatch, tmp_path):
    from types import SimpleNamespace

    from topoppi_affinity.capabilities import predict

    geometry = SimpleNamespace(source={"path": "synthetic.pdb"},
                               partners=[SimpleNamespace(chains=("A",)), SimpleNamespace(chains=("B",))])
    calls = []
    monkeypatch.setattr("topoppi_affinity.workflow.analyze", lambda *args, **kwargs: calls.append(kwargs) or kwargs)
    context = ExecutionContext(plugin_config={"affinity": {"model": "default.joblib", "output": "default-output", "device": "cpu"}},
                               output_dir=tmp_path / "context-output")
    result = predict(geometry, model="chosen.joblib", output=tmp_path / "chosen-output", context=context)
    assert result["model"] == "chosen.joblib"
    assert result["output"] == tmp_path / "chosen-output"
    assert result["device"] == "cpu" and len(calls) == 1
