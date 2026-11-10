"""Desktop contracts: task setup, preserved sample IDs and usable partial results."""
import os
import threading
import time
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from topoppi_affinity.gui import AffinityPanel, variant_lines

FIXTURE = Path(__file__).parent / "fixtures/tiny_complex.pdb"


def write_model(path, *, target="pkd", kind="fitted"):
    import joblib

    bundle = {"model": None, "feature_names": ["x"], "report": {"target": target}}
    if kind == "known":
        bundle.update(base={}, residual={})
    elif kind == "unseen":
        bundle.update(base={}, prior_method="test")
    joblib.dump(bundle, path)
    return path


def test_variant_editor_preserves_csv_ids_and_plain_row_joint_semantics():
    frame = variant_lines("sample_alpha\tA:G1W,B:G1A\nA:G1F\n")
    assert frame.variant_id.tolist() == ["sample_alpha", "variant_002"]
    assert frame.mutations.tolist() == ["A:G1W,B:G1A", "A:G1F"]
    assert variant_lines("\n") is None
    with pytest.raises(ValueError, match="unique"):
        variant_lines("same\tA:G1W\nsame\tA:G1F")


def test_bundled_default_is_ready_and_custom_selection_preserves_inputs(panel, tmp_path, monkeypatch):
    root, view = panel
    local = tmp_path / "provided.joblib"
    write_model(local)
    view.vars["structure"].set(str(FIXTURE))
    view.variants.insert("1.0", "sample1\tA:G1W")
    root.update()
    assert view.run_button["text"] == "Predict"
    assert "ready" in view.model_status.get().lower()
    assert not view.prepare_button.winfo_ismapped()
    assert str(view.open_button["state"]) == "normal"
    monkeypatch.setattr("topoppi_affinity.gui.filedialog.askopenfilename", lambda **_: str(local))
    monkeypatch.setattr("topoppi_affinity.resources.download_models", lambda **_: pytest.fail("Unconfigured source downloaded"))
    view.choose_local_model()
    assert view.model_source.get() == "Local model file"
    assert view.task.get() == "Binding affinity (pKd)"
    assert view.vars["model"].get() == str(local)
    assert view.vars["structure"].get() == str(FIXTURE)
    assert "sample1" in view.variants.get("1.0", "end")
    assert view.run_button["text"] == "Predict" and view.worker is None


@pytest.mark.parametrize("value", [float("nan"), float("inf"), ""])
def test_missing_prediction_values_have_readable_table_labels(panel, tmp_path, monkeypatch, value):
    _, view = panel
    document = _failed_document(tmp_path)
    document.predictions.loc[0, "predicted"] = value
    monkeypatch.setattr("topoppi_affinity.workflow.load_analysis", lambda _: document)
    view.load(tmp_path)
    assert view.table.item("sample1", "values")[2] == "Unavailable"
    assert "Residue absent" in view.details.get()


def test_headless_gui_reports_command_line_alternative(monkeypatch):
    import tkinter as tk

    from topoppi_affinity.gui import main

    def no_display():
        raise tk.TclError("No display")
    monkeypatch.setattr("topoppi_affinity.gui.tk.Tk", no_display)
    with pytest.raises(RuntimeError, match="analyze / inspect / render"):
        main()


@pytest.fixture
def panel(tmp_path, monkeypatch):
    if not os.environ.get("DISPLAY"):
        pytest.skip("Tk display required")
    import tkinter as tk

    monkeypatch.setenv("TOPOPPI_AFFINITY_CACHE", str(tmp_path / "resources"))
    root = tk.Tk()
    root.withdraw()
    monkeypatch.setattr(AffinityPanel, "_missing_dependencies", staticmethod(lambda task: []))
    host = SimpleNamespace(entry_optcuts_bin=SimpleNamespace(get=lambda: "/chosen/OptCuts_bin"))
    view = AffinityPanel(root, host=host)
    try:
        root.update()
        yield root, view
    finally:
        view.cancel.set()
        if view.worker is not None:
            view.worker.join(timeout=5)
        view.close()
        root.destroy()


def _pump(root, predicate, timeout=8):
    until = time.monotonic() + timeout
    while not predicate():
        root.update()
        if time.monotonic() > until:
            pytest.fail("GUI did not process the worker event")
        time.sleep(.01)
    root.update()


def test_discarded_plot_tk_resources_are_collected_before_worker_starts(panel):
    import gc
    import weakref

    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure

    root, view = panel
    owner = threading.get_ident()
    events = []

    def released(name):
        events.append((name, threading.get_ident()))

    def install_plot():
        figure = Figure(figsize=(2, 2))
        figure.subplots().plot([0, 1], [0, 1])
        view.figure = figure
        view.canvas = FigureCanvasTkAgg(figure, master=view.plot)
        toolbar = NavigationToolbar2Tk(view.canvas, view.plot, pack_toolbar=False)
        view.canvas.draw()
        weakref.finalize(view.canvas._tkphoto, released, "image")
        weakref.finalize(toolbar.message, released, "variable")
        return weakref.ref(view.canvas._tkphoto), weakref.ref(toolbar.message)

    image_ref, variable_ref = install_plot()
    # Promote the live view before discarding it so subsequent small GUI
    # allocations do not decide which thread collects its reference cycles.
    gc.collect()
    view._clear_plot()

    def action():
        events.append(("worker", threading.get_ident()))
        gc.collect()

    def finish():
        if view.worker is not None and not view.worker.is_alive():
            root.quit()
        else:
            root.after(10, finish)

    deadline = root.after(5000, root.quit)
    root.after(0, lambda: view._launch(action))
    root.after(10, finish)
    # Run an actual Tk main loop: the regression concerns resource ownership,
    # independently of test helpers that repeatedly call root.update().
    root.mainloop()
    root.after_cancel(deadline)
    assert view.worker is not None and not view.worker.is_alive()
    assert image_ref() is None and variable_ref() is None
    assert {name for name, _ in events} == {"image", "variable", "worker"}
    assert all(thread == owner for name, thread in events if name != "worker")
    worker_index = next(i for i, item in enumerate(events) if item[0] == "worker")
    assert worker_index == 2


@pytest.fixture
def configured_catalog(panel, monkeypatch):
    """Exercise explicit resource-preparation flows with a configured catalog."""
    from topoppi_affinity.resources import model_catalog

    catalog = model_catalog()
    for task, model in catalog["models"].items():
        model["url"] = f"https://example.invalid/{task}.zip"
    monkeypatch.setattr("topoppi_affinity.resources.model_catalog", lambda *args, **kwargs: catalog)
    panel[1]._task_changed()


def test_chain_picker_and_structured_mutation_support_numeric_chain(panel, tmp_path):
    _, view = panel
    content = "\n".join(line[:21] + ("1" if line[21:22] == "A" else line[21:22]) + line[22:]
                        if line.startswith(("ATOM", "HETATM")) else line for line in FIXTURE.read_text().splitlines())
    structure = tmp_path / "numeric.pdb"
    structure.write_text(content + "\n")
    view.vars["structure"].set(str(structure))
    view.read_structure()
    assert view.vars["optcuts_bin"].get() == "/chosen/OptCuts_bin"
    view.chain_choices.selection_set(0)
    view._assign_chains("a")
    assert view.vars["chains_a"].get() == "1"
    assert view.mutation_vars["wt"].get() == "G"
    view.mutation_vars["mt"].set("W")
    view.add_mutation()
    frame = variant_lines(view.variants.get("1.0", "end"))
    assert frame.mutations.tolist() == ["1:G1W"]
    from topoppi_affinity.features import parse_mutations
    assert parse_mutations(frame.mutations.iloc[0])[0]["chain"] == "1"
    view.mutation_vars["position"].set("999")
    assert view.mutation_vars["wt"].get() == ""
    view.vars["structure"].set(str(FIXTURE))
    assert view.structure_records is None
    assert view.chain_choices.size() == 0
    assert view.mutation_vars["wt"].get() == ""


def test_csv_import_keeps_ids_and_export_keeps_units_and_failures(panel, tmp_path, monkeypatch):
    _, view = panel
    source = tmp_path / "variants.csv"
    pd.DataFrame({"variant_id": ["experiment_1", "joint_2"], "mutations": ["A:G1W", "A:G1F,B:G1A"]}).to_csv(source, index=False)
    monkeypatch.setattr("topoppi_affinity.gui.filedialog.askopenfilename", lambda **_: str(source))
    view.import_variants()
    frame = variant_lines(view.variants.get("1.0", "end"))
    assert frame.variant_id.tolist() == ["experiment_1", "joint_2"]
    frame["predicted_ddg_kcal_mol"] = [1.2, ""]
    frame["status"], frame["error"] = ["ok", "error"], ["", "Missing residue"]
    view.document = SimpleNamespace(predictions=frame, path=tmp_path / "analysis.json")
    target = tmp_path / "export.csv"
    monkeypatch.setattr("topoppi_affinity.gui.filedialog.asksaveasfilename", lambda **_: str(target))
    view.export_predictions()
    restored = pd.read_csv(target, keep_default_na=False)
    assert restored.columns.tolist() == frame.columns.tolist()
    assert restored.error.tolist() == ["", "Missing residue"]


@pytest.mark.parametrize("partner_a, expected_chain", [("B", "B"), ("A,B", "A")])
def test_manual_mutation_builder_starts_in_selected_partners(panel, partner_a, expected_chain):
    _, view = panel
    view.vars["structure"].set(str(FIXTURE))
    view.vars["chains_a"].set(partner_a)
    view.vars["chains_b"].set("C")
    view.task.set("Mutation effects (ddG)")
    records = {"chains": [{"id": chain, "residue_count": 1} for chain in "ABC"],
               "residues": [{"chain": chain, "position": 1, "insertion": "", "wt": "G"} for chain in "ABC"]}
    view.read_structure(records=records)
    assert tuple(view.mutation_chain.cget("values")) == ("A", "B", "C")
    assert view.mutation_vars["chain"].get() == expected_chain
    view.mutation_vars["mt"].set("W")
    view.add_mutation()
    variants = variant_lines(view.variants.get("1.0", "end-1c"))
    assert variants.mutations.tolist() == [f"{expected_chain}:G1W"]
    view._check_structure_inputs({"structure": str(FIXTURE), "chains_a": partner_a,
                                  "chains_b": "C", "variants": variants, "task": "unseen"})


def _failed_document(tmp_path, *, known=False):
    predictions = pd.DataFrame([{"variant_id": "sample1", "mutations": "A:G1W", "predicted": "",
                                 "status": "base_only" if known else "error", "error": "Residue absent", "sites": "[]"}])
    report = {"prediction_status": "ok" if known else "error", "layout_status": "unavailable",
              "prediction_error": "Model could not be read", "layout": {}, "input": {}}
    if known:
        report.update(model={"target": "ddg", "kind": "known", "scope": "Known-system ensemble"},
                      scope_assignment={"user_basis": "Matched development system"},
                      prior={"error": "Prior unavailable"})
    def missing_view(*args):
        raise ValueError("No prepared interface")
    return SimpleNamespace(path=tmp_path / "analysis.json", report=report, predictions=predictions, interface=None,
                           selection=lambda *_: (predictions.iloc[0], [], []), view=missing_view)


def test_unknown_target_and_view_error_do_not_hide_prediction_failure(panel, tmp_path, monkeypatch):
    _, view = panel
    document = _failed_document(tmp_path)
    monkeypatch.setattr("topoppi_affinity.workflow.load_analysis", lambda _: document)
    view.load(tmp_path)
    assert "target not yet determined" in view.summary.get()
    assert "ddG" not in view.summary.get()
    assert "Model could not be read" in view.errors.get()
    assert "Residue absent" in view.details.get()
    assert "No prepared interface" in view.details.get()
    assert view.figure is None
    assert "Reuse inputs" in view.status.get()
    assert "review its interface" not in view.status.get()


def test_base_only_and_scope_are_visible_above_results(panel, tmp_path, monkeypatch):
    _, view = panel
    document = _failed_document(tmp_path, known=True)
    monkeypatch.setattr("topoppi_affinity.workflow.load_analysis", lambda _: document)
    view.load(tmp_path)
    assert "BASE MODEL ONLY" in view.summary.get()
    assert "Matched development system" in view.summary.get()
    assert "Prior unavailable" in view.errors.get()
    assert view.table.item("sample1", "values")[-1] == "Base model only"


@pytest.mark.parametrize("attempt", ["error", "cancelled"])
def test_failed_layout_attempt_is_visible_while_retained_atlas_stays_viewable(panel, tmp_path, monkeypatch, attempt):
    from matplotlib.figure import Figure

    _, view = panel
    document = _failed_document(tmp_path)
    warning = f"Latest atlas generation {attempt}; earlier atlas remains viewable."
    document.report.update(layout_status="ok", layout_attempt_status=attempt,
                           warnings=[warning, "Training applicability remains unverified"])
    calls = []
    def retained_view(side, mode, variant, descriptor):
        calls.append(mode)
        return [object()], SimpleNamespace(plot_patches=lambda *args, **kwargs: Figure(figsize=(2, 2))), {}
    document.view = retained_view
    monkeypatch.setattr("topoppi_affinity.workflow.load_analysis", lambda _: document)
    view.mode.set("atlas")
    view.load(tmp_path)
    assert "Atlas: Ready" in view.summary.get()
    expected = "Generation failed" if attempt == "error" else "Cancelled"
    assert f"Latest atlas attempt: {expected}" in view.summary.get()
    assert warning in view.errors.get()
    assert "Training applicability remains unverified" in view.errors.get()
    assert "Model could not be read" in view.errors.get()
    assert calls[-1] == "atlas" and view.figure is not None
    document.report.update(layout_attempt_status="ok", warnings=["Training applicability remains unverified"])
    view.load(tmp_path)
    assert warning not in view.errors.get()
    assert "Latest atlas attempt: Ready" in view.summary.get()


def test_new_run_clears_old_results_and_shows_prediction_before_layout_finishes(panel, tmp_path, monkeypatch):
    root, view = panel
    document = _failed_document(tmp_path)
    monkeypatch.setattr("topoppi_affinity.workflow.load_analysis", lambda _: document)
    view.load(tmp_path)
    release = threading.Event()
    submitted = {}
    def analyze(**kwargs):
        submitted.update(kwargs)
        kwargs["event_callback"]({"stage": "prediction_ready", "status": "ok", "analysis_json": str(document.path),
                                  "message": "Prediction ready; layout still running"})
        release.wait(5)
        return {"output_directory": str(tmp_path)}
    monkeypatch.setattr("topoppi_affinity.workflow.analyze", analyze)
    monkeypatch.setattr("topoppi_affinity.resources.download_models", lambda **_: pytest.fail("Custom model triggered an unrelated resource download"))
    view.vars["structure"].set(str(FIXTURE))
    view.model_source.set("Local model file")
    write_model(tmp_path / "custom.joblib")
    view.vars["model"].set(str(tmp_path / "custom.joblib"))
    view.vars["resource_cache"].set(str(tmp_path / "custom-cache"))
    try:
        view.start()
        assert view.document is None
        assert not view.table.get_children()
        _pump(root, lambda: view.document is document)
        assert view.worker.is_alive()
        assert view.table.get_children() == ("sample1",)
        assert submitted["device"] == "auto"
        assert submitted["retry"] is True
        assert submitted["cache_dir"] == str(tmp_path / "custom-cache")
        view.cancel_run()
        assert "Cancellation requested" in view.status.get()
        assert view.cancel.is_set()
    finally:
        release.set()
        _pump(root, lambda: str(view.run_button["state"]) == "normal")


def test_prepare_uses_selected_task_and_shared_resource_progress(panel, monkeypatch, configured_catalog):
    root, view = panel
    calls = []
    def download_models(**kwargs):
        calls.append(kwargs)
        kwargs["progress"]({"event": "progress", "resource": "ESMC", "current": 1, "total": 2,
                            "message": "Downloading from official ESMC source"})
        return {}
    monkeypatch.setattr("topoppi_affinity.resources.download_models", download_models)
    view.task.set("Mutation effects (ddG)")
    view.mutation_context.set("unseen")
    view.prepare_models()
    _pump(root, lambda: str(view.run_button["state"]) == "normal")
    assert calls[0]["tasks"] == ("unseen",)
    assert calls[0]["include_esmc"] is True
    assert "resources are ready" in view.status.get()


def test_layout_retry_uses_saved_result_and_preserves_model_free_review(panel, tmp_path, monkeypatch):
    root, view = panel
    document = _failed_document(tmp_path)
    monkeypatch.setattr("topoppi_affinity.workflow.load_analysis", lambda _: document)
    view.load(tmp_path)
    calls = []
    def layout_existing(path, **kwargs):
        calls.append((path, kwargs))
        return document.report
    monkeypatch.setattr("topoppi_affinity.workflow.layout_existing", layout_existing)
    view.retry_layout()
    _pump(root, lambda: str(view.run_button["state"]) == "normal")
    assert calls[0][0] == document.path
    assert calls[0][1]["optcuts_bin"] == "/chosen/OptCuts_bin"
    assert view.document.predictions.to_json() == document.predictions.to_json()


def test_invalid_numeric_settings_do_not_start_worker_or_clear_saved_result(panel, tmp_path, monkeypatch, configured_catalog):
    _, view = panel
    sentinel = SimpleNamespace(path=tmp_path / "analysis.json")
    view.document = sentinel
    view.vars["structure"].set(str(FIXTURE))
    view.vars["temperature"].set("nan")
    errors = []
    monkeypatch.setattr("topoppi_affinity.gui.messagebox.showerror", lambda title, text, **_: errors.append(text))
    view.start()
    assert "finite and positive" in errors[0]
    assert view.worker is None
    assert view.document is sentinel
    view.vars["temperature"].set("298.15")
    view.vars["layout_timeout"].set("invalid old layout setting")
    assert view._values()["layout_timeout"] == 600
    view.layout.set(True)
    with pytest.raises(ValueError):
        view._values()


def test_small_window_keeps_primary_actions_visible(panel):
    root, view = panel
    view.window.geometry("900x620")
    root.update()
    right = view.window.winfo_rootx() + view.window.winfo_width()
    bottom = view.window.winfo_rooty() + view.window.winfo_height()
    for widget in (view.run_button, view.cancel_button, view.layout_button, view.partner_selector):
        assert widget.winfo_rootx() + widget.winfo_width() <= right
        assert widget.winfo_rooty() + widget.winfo_height() <= bottom


def test_repair_opens_shared_installer_and_requires_restart_before_compute(panel, monkeypatch):
    _, view = panel
    calls = []
    monkeypatch.setattr("topoppi_affinity.installer.open_installer", lambda host, **kwargs: calls.append((host, kwargs)))
    view.repair_dependencies()
    assert calls == [(view, {"download": False})]
    view._affinity_installer = SimpleNamespace(worker=None, restart_required=True, window=SimpleNamespace(winfo_exists=lambda: False))
    assert view._can_compute() is False
    assert "Restart TopoPPI" in view.errors.get()


def test_dependency_restart_restores_edited_task_variants_and_thermodynamics(panel):
    _, view = panel
    view.task.set("Mutation effects (ddG)")
    view.mutation_context.set("known")
    view.vars["structure"].set(str(FIXTURE))
    view.vars["temperature"].set("310")
    view.vars["wt_kd"].set("1e-8")
    view.vars["scope_note"].set("A matched development system")
    view.variants.insert("1.0", "sample_with_ID\tA:G1W,B:G1A")
    view.layout.set(True)
    setup = view.export_setup()
    view.vars["temperature"].set("298")
    view.vars["wt_kd"].set("")
    view.variants.delete("1.0", "end")
    view.restore_setup(setup)
    assert view.export_setup() == setup


def test_load_installed_example_prefills_successful_task_and_reuses_original_files(panel, tmp_path, monkeypatch):
    _, view = panel
    destination = tmp_path / "example"
    monkeypatch.setattr("topoppi_affinity.gui.filedialog.askdirectory", lambda **_: str(destination))
    def no_download(**_):
        pytest.fail("Loading the installed example must not download model resources")
    monkeypatch.setattr("topoppi_affinity.resources.download_models", no_download)
    view.task.set("Mutation effects (ddG)")
    view.mutation_context.set("known")
    view.vars["wt_kd"].set("1e-8")
    view.vars["prior_csv"].set("missing-prior.csv")
    view.variants.insert("1.0", "old_variant\tA:G1W")
    view.mode.set("atlas")
    view.layout.set(True)
    view.load_example()
    assert all((destination / name).is_file() for name in ("1bvk.pdb", "variants.csv", "README.md"))
    values = view._values()
    assert values["structure"] == str(destination / "1bvk.pdb")
    assert (values["chains_a"], values["chains_b"], values["task"]) == ("A,B", "C", "absolute")
    assert values["variants"] is None and values["wt_kd"] is None and values["prior_csv"] is None
    assert values["layout"] is False and values["retry"] is True
    assert view.mode.get() == "surface" and view.document is None
    # The crystallographic file contains another assembly; the example selects A,B/C.
    assert {"A", "B", "C"}.issubset({row["id"] for row in view.structure_records["chains"]})
    original = (destination / "1bvk.pdb").read_bytes()
    retained = destination / "results" / "user-result.txt"
    retained.parent.mkdir()
    retained.write_text("Keep this result")
    view.load_example()
    assert (destination / "1bvk.pdb").read_bytes() == original
    assert retained.read_text() == "Keep this result"


def test_load_example_does_not_overwrite_edited_structure_or_change_current_form(panel, tmp_path, monkeypatch):
    _, view = panel
    from topoppi_affinity.resources import export_example

    destination = tmp_path / "example"
    export_example(destination)
    changed = (destination / "1bvk.pdb").read_text() + "REMARK User edited this structure\n"
    (destination / "1bvk.pdb").write_text(changed)
    monkeypatch.setattr("topoppi_affinity.gui.filedialog.askdirectory", lambda **_: str(destination))
    view.vars["structure"].set(str(FIXTURE))
    view.variants.insert("1.0", "keep_ID\tA:G1W")
    before = view.export_setup()
    view.load_example()
    assert view.export_setup() == before
    assert (destination / "1bvk.pdb").read_text() == changed
    assert "choose a new directory" in view.errors.get()
    assert "current inputs" in view.status.get()


def test_task_switch_shows_relevant_inputs_and_ignores_hidden_mutation_settings(panel):
    root, view = panel
    view.vars["structure"].set(str(FIXTURE))
    view.task.set("Mutation effects (ddG)")
    view.mutation_context.set("known")
    root.update()
    assert view.mutation_box.winfo_ismapped()
    assert not view.custom_box.winfo_ismapped()
    view.variants.insert("1.0", "chosen_ID\tA:G1W")
    view.vars["wt_kd"].set("1e-8")
    view.vars["prior_csv"].set("a mutation-specific prior.csv")
    view.task.set("Binding affinity (pKd)")
    root.update()
    assert not view.mutation_box.winfo_ismapped()
    values = view._values()
    assert values["variants"] is None and values["wt_kd"] is None and values["prior_csv"] is None
    view.task.set("Mutation effects (ddG)")
    view.mutation_context.set("known")
    assert variant_lines(view.variants.get("1.0", "end")).variant_id.tolist() == ["chosen_ID"]
    assert view.vars["wt_kd"].get() == "1e-8"


def test_partner_assignment_moves_chains_between_groups_and_preflight_checks_variants(panel):
    _, view = panel
    view.vars["structure"].set(str(FIXTURE))
    view.read_structure()
    view.chain_choices.selection_set(0, 1)
    view._assign_chains("a")
    assert view.vars["chains_a"].get() == "A,B"
    assert view.vars["chains_b"].get() == ""
    view.vars["chains_b"].set("B")
    with pytest.raises(ValueError, match="overlap"):
        view._values()
    view.vars["chains_a"].set("A")
    view.vars["chains_b"].set("X")
    with pytest.raises(ValueError, match="Missing: X"):
        view._values()
    view.vars["chains_b"].set("B")
    view.task.set("Mutation effects (ddG)")
    view.mutation_context.set("unseen")
    with pytest.raises(ValueError, match="at least one substitution"):
        view._values()
    view.variants.insert("1.0", "bad\tA:G999W")
    with pytest.raises(ValueError, match="bad:.*absent"):
        view._values()
    view.variants.insert("end", "\ngood\tA:G1W")
    # Keep row-specific failures when another variant can be predicted.
    assert view._values()["variants"].variant_id.tolist() == ["bad", "good"]


def test_predict_prepares_selected_resources_then_runs_without_second_click(panel, tmp_path, monkeypatch, configured_catalog):
    root, view = panel
    document = _failed_document(tmp_path, known=True)
    monkeypatch.setattr("topoppi_affinity.workflow.load_analysis", lambda _: document)
    view.vars["structure"].set(str(FIXTURE))
    view.task.set("Mutation effects (ddG)")
    view.mutation_context.set("unseen")
    view.variants.insert("1.0", "experiment\tA:G1W")
    calls = []
    def download(**kwargs):
        calls.append(("prepare", kwargs))
        kwargs["progress"]({"resource": "ESMC", "current": 5, "total": 10})
    def analyze(**kwargs):
        calls.append(("predict", kwargs))
        return {"output_directory": str(tmp_path)}
    monkeypatch.setattr("topoppi_affinity.resources.download_models", download)
    monkeypatch.setattr("topoppi_affinity.workflow.analyze", analyze)
    assert view.run_button["text"] == "Prepare and predict"
    view.start()
    _pump(root, lambda: not view._working)
    assert [name for name, _ in calls] == ["prepare", "predict"]
    assert calls[0][1]["tasks"] == ("unseen",)
    assert calls[0][1]["include_esmc"] is True
    assert calls[1][1]["variants"].variant_id.tolist() == ["experiment"]
    assert view.document is document


def test_preparation_cancellation_stops_before_prediction_and_restores_inputs(panel, monkeypatch, configured_catalog):
    root, view = panel
    view.vars["structure"].set(str(FIXTURE))
    view.task.set("Mutation effects (ddG)")
    view.mutation_context.set("unseen")
    view.variants.insert("1.0", "experiment\tA:G1W")
    entered, release = threading.Event(), threading.Event()
    def download(**kwargs):
        entered.set()
        release.wait(5)
    monkeypatch.setattr("topoppi_affinity.resources.download_models", download)
    monkeypatch.setattr("topoppi_affinity.workflow.analyze", lambda **_: pytest.fail("Prediction started after cancellation"))
    try:
        view.start()
        _pump(root, entered.is_set)
        assert str(view.task_selector["state"]) == "disabled"
        view.cancel_run()
        release.set()
        _pump(root, lambda: not view._working)
        assert str(view.task_selector["state"]) == "readonly"
        assert view.vars["structure"].get() == str(FIXTURE)
        assert "Cancelled" in view.status.get()
        assert str(view.progress["mode"]) == "determinate"
        assert float(view.progress["value"]) == 0
    finally:
        release.set()


def test_prior_override_prepares_only_model_and_embeddings(panel, tmp_path, monkeypatch, configured_catalog):
    root, view = panel
    view.task.set("Mutation effects (ddG)")
    view.mutation_context.set("known")
    (tmp_path / "prior.csv").touch()
    view.vars["prior_csv"].set(str(tmp_path / "prior.csv"))
    calls = []
    monkeypatch.setattr("topoppi_affinity.resources.download_models", lambda **kwargs: calls.append(kwargs))
    view.prepare_models()
    _pump(root, lambda: not view._working)
    assert calls[0]["tasks"] == ("known",)
    assert calls[0]["include_priors"] is False


def test_missing_dependencies_open_selected_task_setup_before_download(panel, monkeypatch, configured_catalog):
    _, view = panel
    monkeypatch.setattr(view, "_missing_dependencies", lambda task: ["ESMC"] if task else [])
    calls = []
    monkeypatch.setattr("topoppi_affinity.installer.open_installer", lambda host, **kwargs: calls.append((host, kwargs)))
    monkeypatch.setattr("topoppi_affinity.resources.download_models", lambda **_: pytest.fail("Downloaded before dependency setup"))
    view.task.set("Mutation effects (ddG)")
    view.mutation_context.set("unseen")
    view.variants.insert("1.0", "saved_ID\tA:G1W")
    before = view.export_setup()
    assert view.run_button["text"] == "Set up prediction..."
    view.start()
    assert calls == [(view, {"download": True, "tasks": ("unseen",)})]
    assert view.export_setup() == before
    assert view.worker is None
    view.model_source.set("Local model file")
    assert view._prediction_dependencies_ready()


def test_task_only_restart_preset_and_model_free_empty_review(panel):
    root, view = panel
    view.restore_setup({"task": "known"})
    root.update()
    assert view.task.get() == "Mutation effects (ddG)"
    assert view.mutation_context.get() == "known"
    assert view.mutation_box.winfo_ismapped()
    assert view.vars["temperature"].get() == "298.15"
    assert str(view.layout_button["state"]) == "disabled"
    assert str(view.linked_button["state"]) == "disabled"
    assert str(view.open_button["state"]) == "normal"


@pytest.mark.parametrize("key", ["prior_csv", "prior_checkpoint", "prior_repository"])
def test_invalid_prior_path_stops_before_resource_download(panel, tmp_path, monkeypatch, key, configured_catalog):
    _, view = panel
    view.vars["structure"].set(str(FIXTURE))
    view.task.set("Mutation effects (ddG)")
    view.mutation_context.set("known")
    view.variants.insert("1.0", "A:G1W")
    view.vars[key].set(str(tmp_path / "missing"))
    errors = []
    monkeypatch.setattr("topoppi_affinity.gui.messagebox.showerror", lambda title, text, **_: errors.append(text))
    monkeypatch.setattr("topoppi_affinity.resources.download_models", lambda **_: pytest.fail("Downloaded with an invalid prior path"))
    view.start()
    assert errors and "does not exist" in errors[0]
    assert view.worker is None
    view.prepare_models()
    assert "does not exist" in view.errors.get()
    assert view.worker is None


def test_host_output_directory_is_used_for_affinity_subdirectory(tmp_path):
    if not os.environ.get("DISPLAY"):
        pytest.skip("Tk display required")
    import tkinter as tk

    root = tk.Tk()
    root.withdraw()
    host = SimpleNamespace(var_output_dir=SimpleNamespace(get=lambda: str(tmp_path)))
    view = AffinityPanel(root, host=host)
    try:
        assert Path(view.vars["output"].get()) == tmp_path / "affinity"
    finally:
        view.close()
        root.destroy()


def test_ready_resources_predict_directly_without_downloading(panel, tmp_path, monkeypatch):
    from topoppi_affinity.resources import resource_status

    root, view = panel
    ready = resource_status()
    ready["esmc"]["installed"] = True
    for model in ready["models"].values():
        model.update(ready=True, installed=True, prior_ready=True, missing_resources=[])
    monkeypatch.setattr("topoppi_affinity.resources.resource_status", lambda **_: ready)
    monkeypatch.setattr("topoppi_affinity.resources.download_models", lambda **_: pytest.fail("Ready resources were downloaded again"))
    document = _failed_document(tmp_path)
    monkeypatch.setattr("topoppi_affinity.workflow.load_analysis", lambda _: document)
    calls = []
    def analyze(**kwargs):
        calls.append(kwargs)
        return {"output_directory": str(tmp_path)}
    monkeypatch.setattr("topoppi_affinity.workflow.analyze", analyze)
    view.vars["structure"].set(str(FIXTURE))
    view._task_changed()
    assert view.run_button["text"] == "Predict"
    view.start()
    _pump(root, lambda: not view._working)
    assert len(calls) == 1 and calls[0]["task"] == "absolute"


def test_clearing_plot_cancels_queued_redraw(panel, tmp_path, monkeypatch):
    from matplotlib.figure import Figure
    root, view = panel
    document = _failed_document(tmp_path)
    document.view = lambda *_: ([object()], SimpleNamespace(
        plot_patches=lambda *args, **kwargs: Figure(figsize=(2, 2))), {})
    monkeypatch.setattr("topoppi_affinity.workflow.load_analysis", lambda _: document)
    view.load(tmp_path)
    pending = view.canvas._idle_draw_id
    assert pending in root.tk.call('after', 'info')
    view._clear_plot()
    assert pending not in root.tk.call('after', 'info')
    assert view.figure is None and view.canvas is None


def _main_window_fields(view, *, structure, chains_a="A", chains_b="B", output=""):
    fields = {"var_input_path": structure, "var_chain_a": chains_a, "var_chain_b": chains_b, "var_output_dir": output}
    for name, value in fields.items():
        setattr(view.host, name, SimpleNamespace(get=lambda value=value: str(value)))
    view.host._affinity_panel = view
    return view.host


def test_reopening_affinity_exposes_main_input_change_without_overwriting_setup(panel, tmp_path):
    from topoppi_affinity.gui import open_panel

    _, view = panel
    view.vars["structure"].set(str(FIXTURE))
    view.task.set("Mutation effects (ddG)")
    view.mutation_context.set("known")
    view.variants.insert("1.0", "edited_ID\tA:G1W")
    view.vars["temperature"].set("310")
    newer = tmp_path / "new-complex.pdb"
    newer.write_text(FIXTURE.read_text())
    host = _main_window_fields(view, structure=newer)
    before = view.export_setup()
    assert open_panel(host) is view
    assert view.export_setup() == before
    assert "new-complex.pdb" in view.main_input_context.get()
    assert "next prediction" in view.main_input_context.get()
    assert str(view.use_main_input_button["state"]) == "normal"


def test_explicit_main_input_handoff_preserves_variants_model_settings_and_displayed_result(panel, tmp_path, monkeypatch):
    _, view = panel
    document = _failed_document(tmp_path, known=True)
    document.report["input"] = {"path": str(tmp_path / "displayed-complex.pdb")}
    monkeypatch.setattr("topoppi_affinity.workflow.load_analysis", lambda _: document)
    view.load(document.path)
    view.task.set("Mutation effects (ddG)")
    view.mutation_context.set("known")
    view.variants.insert("1.0", "edited_ID\tA:G1W")
    view.vars["temperature"].set("310")
    view.vars["scope_note"].set("Matched system")
    summary, rows = view.summary.get(), view.table.get_children()
    content = "\n".join(line[:21] + ("1" if line[21:22] == "A" else line[21:22]) + line[22:]
                        if line.startswith(("ATOM", "HETATM")) else line for line in FIXTURE.read_text().splitlines())
    newer = tmp_path / "numeric-complex.pdb"
    newer.write_text(content + "\n")
    _main_window_fields(view, structure=newer, chains_a="1", output=tmp_path / "figures")
    view.use_main_window_input()
    assert view.vars["structure"].get() == str(newer)
    assert (view.vars["chains_a"].get(), view.vars["chains_b"].get()) == ("1", "B")
    assert Path(view.vars["output"].get()) == tmp_path / "figures" / "affinity"
    assert view.task.get() == "Mutation effects (ddG)"
    assert view.mutation_context.get() == "known"
    assert view.vars["temperature"].get() == "310" and view.vars["scope_note"].get() == "Matched system"
    assert view.variants.get("1.0", "end-1c") == "edited_ID\tA:G1W"
    assert view.document is document and view.table.get_children() == rows and view.summary.get() == summary
    assert "displayed-complex.pdb" in summary
    assert "Variant rows are retained" in view.status.get() and "displayed result is unchanged" in view.status.get()
    assert view.mutation_vars["chain"].get() == "1"


def test_main_input_handoff_is_gated_during_work_and_preserves_setup_for_unreadable_structure(panel, tmp_path):
    _, view = panel
    view.vars["structure"].set(str(FIXTURE))
    _main_window_fields(view, structure=tmp_path / "missing.pdb")
    before = view.export_setup()
    view._busy(True)
    try:
        view._refresh_main_input_context()
        assert str(view.use_main_input_button["state"]) == "disabled"
        view.use_main_window_input()
        assert view.export_setup() == before
    finally:
        view._busy(False)
    assert str(view.use_main_input_button["state"]) == "normal"
    view.use_main_window_input()
    assert view.export_setup() == before
    assert "Check the structure selected in the main window" in view.errors.get()


def test_prediction_action_is_primary_and_cancel_tracks_active_work(panel):
    root, view = panel
    assert view.run_button.winfo_width() > view.prepare_button.winfo_width()
    assert not view.cancel_button.winfo_ismapped()
    view._busy(True)
    root.update()
    assert view.cancel_button.winfo_ismapped()
    assert view.progress.winfo_ismapped()
    view._busy(False)
    root.update()
    assert not view.cancel_button.winfo_ismapped()
    assert not view.progress.winfo_ismapped()


def test_scientific_task_and_model_source_keep_optional_advanced_scope(panel):
    _, view = panel
    assert tuple(view.task_selector['values']) == ('Binding affinity (pKd)', 'Mutation effects (ddG)')
    view.restore_setup({'prediction_task': 'mutation'})
    assert view.mutation_context.get() == ''
    view.vars['structure'].set(str(FIXTURE))
    view.variants.insert('1.0', 'draft\tA:G1W')
    assert view._values()['task'] == 'mutation'
    view.mutation_context.set('known')
    view.restore_setup({'prediction_task': 'absolute'})
    assert view._values()['variants'] is None
    view.restore_setup({'prediction_task': 'mutation'})
    assert view.mutation_context.get() == 'known'
    assert view._values()['task'] == 'known'
    assert view.variants.get('1.0', 'end-1c') == 'draft\tA:G1W'


def test_mutation_model_selector_can_restore_bundled_default(panel):
    root, view = panel
    view.restore_setup({'prediction_task': 'mutation'})
    view.tabs.select(1)
    root.update()
    assert view.scope_selector.winfo_ismapped()
    view.scope_selector.set('Unseen-system neural ensemble')
    view.scope_selector.event_generate('<<ComboboxSelected>>')
    root.update()
    assert view.mutation_context.get() == 'unseen'
    view.scope_selector.set('Bundled physical model')
    view.scope_selector.event_generate('<<ComboboxSelected>>')
    root.update()
    assert view.mutation_context.get() == ''
    assert view._workflow_task() == 'mutation'
    assert view.run_button['text'] == 'Predict'


def test_bundled_model_preparation_keeps_network_resources_optional(panel, monkeypatch):
    _, view = panel
    monkeypatch.setattr('topoppi_affinity.resources.download_models',
                        lambda **_: pytest.fail('Bundled defaults requested optional neural weights'))
    for task in ('absolute', 'mutation'):
        view.restore_setup({'prediction_task': task})
        view.prepare_models()
        assert view.worker is None


def test_local_model_selection_infers_target_scope_and_detects_later_mismatch(panel, tmp_path, monkeypatch):
    _, view = panel
    model = write_model(tmp_path / 'known.joblib', target='ddg', kind='known')
    view.vars['structure'].set(str(FIXTURE))
    view.variants.insert('1.0', 'draft\tA:G1W')
    monkeypatch.setattr('topoppi_affinity.gui.filedialog.askopenfilename', lambda **_: str(model))
    view.choose_local_model()
    assert view.task.get() == 'Mutation effects (ddG)'
    assert view.model_source.get() == 'Local model file'
    assert view.mutation_context.get() == 'known'
    assert 'ddG (kcal/mol)' in view.local_model_details.get()
    assert view._values()['task'] is None
    assert view._values()['mutation_context'] == 'known'
    view.mutation_context.set('unseen')
    with pytest.raises(ValueError, match='for known systems'):
        view._values()
    view.restore_setup({'prediction_task': 'absolute'})
    assert 'does not match' in view.model_status.get()
    with pytest.raises(ValueError, match='task needs pkd'):
        view._values()
    assert view.variants.get('1.0', 'end-1c') == 'draft\tA:G1W'


def test_unscoped_local_model_does_not_establish_training_membership(panel, tmp_path):
    _, view = panel
    view.vars['model'].set(str(write_model(tmp_path / 'custom.joblib', target='ddg')))
    view.inspect_local_model()
    assert view.task.get() == 'Mutation effects (ddG)'
    assert view.mutation_context.get() == ''
    assert 'not declared' in view.local_model_details.get()


def _handoff(*rows, structure=FIXTURE, **extra):
    return {'prediction_task': 'mutation', 'structure': str(structure), 'chains_a': 'A', 'chains_b': 'B',
            'append_variants': [{'variant_id': name, 'mutations': mutations} for name, mutations in rows], **extra}


def test_residue_handoff_appends_single_and_joint_variants_without_changing_model_setup(panel):
    _, view = panel
    view.vars['structure'].set(str(FIXTURE))
    view.variants.insert('1.0', 'selected\tA:G1W')
    view.mutation_context.set('known')
    view.vars['temperature'].set('310')
    view.restore_setup(_handoff(('selected', 'A:G1F'), ('joint', 'A:G1W,B:G1A')))
    rows = variant_lines(view.variants.get('1.0', 'end'))
    assert rows.variant_id.tolist() == ['selected', 'selected_2', 'joint']
    assert rows.mutations.tolist() == ['A:G1W', 'A:G1F', 'A:G1W,B:G1A']
    assert view.mutation_context.get() == 'known'
    assert view.vars['temperature'].get() == '310'
    assert 'Added 2 variants' in view.status.get()


@pytest.mark.parametrize('conflict', ['structure', 'partners', 'worker', 'invalid_residue'])
def test_residue_handoff_rejects_conflicts_atomically(panel, tmp_path, conflict):
    _, view = panel
    view.vars['structure'].set(str(FIXTURE))
    view.variants.insert('1.0', 'my_draft\tA:G1W')
    incoming = _handoff(('selected', 'A:G1F'))
    if conflict == 'structure':
        newer = tmp_path / 'different.pdb'
        newer.write_text(FIXTURE.read_text())
        incoming['structure'] = str(newer)
    elif conflict == 'partners':
        incoming.update(chains_a='B', chains_b='A')
    elif conflict == 'worker':
        view._busy(True)
    else:
        incoming['append_variants'][0]['mutations'] = 'A:G999F'
    before = view.export_setup()
    try:
        with pytest.raises((ValueError, RuntimeError)):
            view.restore_setup(incoming)
        assert view.export_setup() == before
    finally:
        if conflict == 'worker':
            view._busy(False)


def test_first_handoff_uses_displayed_source_instead_of_main_input_defaults(panel, tmp_path):
    _, view = panel
    main = tmp_path / 'different-main-input.pdb'
    main.write_text(FIXTURE.read_text())
    view.vars['structure'].set(str(main))
    view.restore_setup(_handoff(('selected', 'A:G1W')), initial=True)
    assert view.vars['structure'].get() == str(FIXTURE)
    assert view.mutation_context.get() == ''
    assert view.task.get() == 'Mutation effects (ddG)'


def test_cross_task_navigation_keeps_an_existing_draft_for_another_source(panel, tmp_path):
    _, view = panel
    view.vars['structure'].set(str(FIXTURE))
    view.variants.insert('1.0', 'draft\tA:G1W')
    other = tmp_path / 'other.pdb'
    other.write_text(FIXTURE.read_text())
    before = view.export_setup()
    with pytest.raises(ValueError, match='different structure'):
        view.restore_setup({'prediction_task': 'mutation', 'structure': str(other), 'chains_a': 'A', 'chains_b': 'B'})
    assert view.export_setup() == before


def test_closing_and_reopening_retains_session_draft_and_rejects_a_different_source(panel, tmp_path):
    from topoppi_affinity.gui import open_panel
    from topoppi_affinity.integration import AffinityController

    root, view = panel
    view.host.root = root
    view.host._affinity_panel = view
    view.vars['structure'].set(str(FIXTURE))
    view.vars['model'].set(str(write_model(tmp_path / 'unseen.joblib', target='ddg', kind='unseen')))
    view.inspect_local_model()
    view.variants.insert('1.0', 'single\tA:G1W\njoint\tA:G1W,B:G1A')
    view.vars['temperature'].set('310')
    view.vars['wt_kd'].set('1e-8')
    before = view.export_setup()
    assert view.close()
    restored = open_panel(view.host)
    try:
        assert restored.export_setup() == before
        assert 'restored from this session' in restored.status.get()
        assert restored.close()
        other = tmp_path / 'other.pdb'
        other.write_text(FIXTURE.read_text())
        view.host._open_affinity = lambda: open_panel(view.host)
        with pytest.raises(ValueError, match='different structure'):
            AffinityController(view.host).open({'prediction_task': 'mutation', 'structure': str(other),
                                               'chains_a': 'A', 'chains_b': 'B'})
        assert view.host._affinity_panel.export_setup() == before
    finally:
        view.host._affinity_panel.close()


def test_csv_append_preserves_ids_and_replace_is_explicit_and_undoable(panel, tmp_path, monkeypatch):
    root, view = panel
    view.restore_setup({'prediction_task': 'mutation'})
    root.update()
    original = 'existing\tA:G1W'
    view.variants.insert('1.0', original)
    source = tmp_path / 'library.csv'
    source.write_text('variant_id,mutations\nnew_row,B:G1A\n')
    monkeypatch.setattr('topoppi_affinity.gui.filedialog.askopenfilename', lambda **_: str(source))
    view.import_variants()
    appended = view.variants.get('1.0', 'end-1c')
    assert variant_lines(appended).variant_id.tolist() == ['existing', 'new_row']
    view.import_variants()
    assert view.variants.get('1.0', 'end-1c') == appended
    assert 'already exist' in view.errors.get() and 'new_row' in view.errors.get()
    view.import_mode.set('replace')
    view.import_variants()
    assert view.variants.get('1.0', 'end-1c') == 'new_row\tB:G1A'
    assert not view.errors.get()
    view.variants.focus_force()
    view.variants.event_generate('<Control-z>')
    root.update()
    assert view.variants.get('1.0', 'end-1c') == appended


@pytest.mark.parametrize('field,value,label', [('temperature', '310 K', 'Reference temperature'),
                                               ('wt_kd', '10 nM', 'Measured WT Kd'),
                                               ('layout_timeout', '1 minute', 'Layout timeout')])
def test_numeric_input_errors_name_the_field_and_required_number(panel, field, value, label):
    _, view = panel
    view.restore_setup({'prediction_task': 'mutation', 'mutation_context': 'unseen'})
    view.vars['structure'].set(str(FIXTURE))
    view.variants.insert('1.0', 'candidate\tA:G1W')
    view.layout.set(True)
    view.vars[field].set(value)
    with pytest.raises(ValueError, match=label + ' must be a number'):
        view._values()
    view.vars[field].set('310' if field == 'temperature' else '1e-8' if field == 'wt_kd' else '60')
    assert view._values()[field] > 0


def test_successful_export_retry_removes_export_error_but_keeps_result_warning(panel, tmp_path, monkeypatch):
    _, view = panel
    view.document = _failed_document(tmp_path)
    view.errors.set('Scientific result remains unavailable.')
    destination = tmp_path / 'missing-directory' / 'export.csv'
    monkeypatch.setattr('topoppi_affinity.gui.filedialog.asksaveasfilename', lambda **_: str(destination))
    view.export_predictions()
    assert view.errors.get().startswith('Cannot export predictions:')
    destination = tmp_path / 'export.csv'
    view.export_predictions()
    assert destination.is_file()
    assert view.errors.get() == 'Scientific result remains unavailable.'


@pytest.mark.parametrize('kind,saved_context,expected', [('descriptor', None, ''),
                                                       ('descriptor', 'unseen', 'unseen'),
                                                       ('known', None, 'known')])
def test_reuse_inputs_restores_recorded_applicability_without_inheriting_other_work(panel, tmp_path, kind, saved_context, expected):
    _, view = panel
    document = _failed_document(tmp_path)
    document.report.update(model={'kind': kind, 'target': 'ddg'},
                           input={'path': str(FIXTURE), 'chains_a': 'A', 'chains_b': 'B'})
    if saved_context:
        document.report['scope_assignment'] = {'context': saved_context, 'user_basis': 'Saved assessment'}
    view.document = document
    view.mutation_context.set('known')
    view.reuse_inputs()
    assert view.mutation_context.get() == expected


def test_explicit_setup_file_keeps_incomplete_draft_without_model_or_old_result(panel, tmp_path, monkeypatch):
    import json

    _, view = panel
    view.restore_setup({'prediction_task': 'mutation', 'mutation_context': ''})
    view.model_source.set('Local model file')
    view.vars['structure'].set('')
    view.vars['model'].set('')
    view.vars['temperature'].set('310 K')
    view.variants.insert('1.0', 'unfinished\tB:')
    view.document = SimpleNamespace(path=tmp_path / 'unavailable-result' / 'analysis.json')
    before = view.export_setup()
    target = tmp_path / 'unfinished.topoppi-setup.json'
    monkeypatch.setattr('topoppi_affinity.gui.filedialog.asksaveasfilename', lambda **_: str(target))
    monkeypatch.setattr('joblib.load', lambda *_: pytest.fail('Saving/opening a setup must not load model weights'))
    monkeypatch.setattr('topoppi_affinity.workflow.load_analysis', lambda *_: pytest.fail('A setup must not need the old displayed result'))
    view.save_setup()
    stored = json.loads(target.read_text())
    assert stored['format'] == 'topoppi-affinity-setup'
    assert 'analysis' not in stored['setup']
    assert stored['setup']['variants'] == 'unfinished\tB:'
    view.vars['structure'].set(str(FIXTURE))
    view.variants.delete('1.0', 'end')
    view.variants.insert('1.0', 'other\tA:G1W')
    monkeypatch.setattr('topoppi_affinity.gui.filedialog.askopenfilename', lambda **_: str(target))
    view.open_setup()
    assert view.export_setup() == before
    assert 'Setup opened' in view.status.get()


@pytest.mark.parametrize('contents', ['not JSON', '{"format":"topoppi-affinity-analysis","version":1}',
                                    '{"format":"topoppi-affinity-setup","version":1,"setup":{"fields":[]}}'])
def test_wrong_setup_file_reports_recovery_without_changing_edited_fields(panel, tmp_path, monkeypatch, contents):
    _, view = panel
    view.vars['structure'].set(str(FIXTURE))
    view.variants.insert('1.0', 'existing\tA:G1W')
    before = view.export_setup()
    target = tmp_path / 'wrong-file.json'
    target.write_text(contents)
    monkeypatch.setattr('topoppi_affinity.gui.filedialog.askopenfilename', lambda **_: str(target))
    view.open_setup()
    assert view.errors.get().startswith('Cannot open prediction setup:')
    assert view.export_setup() == before


def test_open_setup_is_disabled_during_prediction_and_file_actions_remain_visible(panel):
    root, view = panel
    view.window.geometry('760x540')
    root.update()
    for button in (view.open_setup_button, view.save_setup_button):
        assert button.winfo_ismapped()
        assert button.winfo_rooty() >= view.window.winfo_rooty()
        assert button.winfo_rooty() + button.winfo_height() < view.window.winfo_rooty() + view.window.winfo_height()
    view._busy(True)
    try:
        assert 'disabled' in view.open_setup_button.state()
        assert 'disabled' not in view.save_setup_button.state()
    finally:
        view._busy(False)
