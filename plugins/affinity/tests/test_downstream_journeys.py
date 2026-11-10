"""Actual bundled-model desktop journeys through successive prediction tasks."""
from __future__ import annotations

import socket
import threading
import time
import tkinter as tk

import pytest
from topoppi.examples import copy_example
from topoppi.gui_app.app import InterfaceApp


def _button(parent, label):
    for child in parent.winfo_children():
        if child.winfo_class() in {"TButton", "Button"} and child.cget("text") == label:
            return child
        if found := _button(child, label):
            return found


def _idle(root, panel):
    deadline = time.monotonic() + 60
    while panel._working and time.monotonic() < deadline:
        root.update()
        time.sleep(.01)
    root.update()
    assert not panel._working, "The local CPU prediction did not finish"


def _finish(root, panel):
    _idle(root, panel)
    assert panel.document is not None, panel.errors.get()
    assert panel.document.predictions.status.eq("ok").all(), panel.document.predictions.to_dict("records")


@pytest.fixture
def desktop(tmp_path, monkeypatch):
    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("A graphical display is required")
    monkeypatch.setattr(InterfaceApp, "_recent_store_path", lambda _: tmp_path / "recent.json")
    monkeypatch.setenv("TOPOPPI_AFFINITY_CACHE", str(tmp_path / "resources"))
    callbacks = []
    root.report_callback_exception = lambda *error: callbacks.append(error)
    app = InterfaceApp(root)
    assert {"absolute", "mutation"} <= set(app._binding_tasks), app._extension_errors
    source = copy_example(tmp_path / "input", synthetic=True)
    app.var_input_path.set(str(source))
    app.var_chain_a.set("A")
    app.var_chain_b.set("B")
    app.var_output_dir.set(str(tmp_path))

    def offline(*args, **kwargs):
        raise AssertionError("A bundled CPU prediction attempted a network connection")

    monkeypatch.setattr(socket.socket, "connect", offline)
    try:
        yield app
        assert not callbacks
    finally:
        app.close()


@pytest.mark.parametrize("first_task", ["Binding affinity", "Mutation effects"])
def test_reusing_bundled_prediction_keeps_task_model_for_the_next_task(desktop, tmp_path, first_task):
    app = desktop
    app.btn_binding_analysis.invoke()
    app.root.update()
    _button(app._binding_window, first_task).invoke()
    panel = app._affinity_panel
    if first_task == "Mutation effects":
        panel.variants.insert("1.0", "candidate\tA:F2A")
    panel.vars["output"].set(str(tmp_path / "first-result"))
    panel.run_button.invoke()
    _finish(app.root, panel)
    original = panel.document
    _button(panel.window, "Reuse inputs").invoke()
    app.root.update()
    assert panel.model_source.get() == "Prepared task model"
    assert panel.document is original

    # Reusing a default prediction must leave the other default task usable.
    next_task = "Mutation effects" if first_task == "Binding affinity" else "Binding affinity"
    _button(panel.window, next_task).invoke()
    if next_task == "Mutation effects":
        panel.variants.insert("1.0", "candidate\tA:F2A")
    panel.vars["output"].set(str(tmp_path / "next-result"))
    panel.run_button.invoke()
    _finish(app.root, panel)
    expected = "ddg" if next_task == "Mutation effects" else "pkd"
    assert panel.document.report["model"]["target"] == expected
    assert original.path.is_file()


@pytest.mark.parametrize("problem", ["file", "existing_analysis"])
def test_invalid_results_folder_keeps_the_previous_result_and_allows_correction(desktop, tmp_path, monkeypatch, problem):
    app = desktop
    panel = app.open_binding_task("absolute")
    panel.vars["output"].set(str(tmp_path / "first-result"))
    panel.run_button.invoke()
    _finish(app.root, panel)
    previous = panel.document
    previous_bytes = previous.path.read_bytes()
    worker = panel.worker
    _button(panel.window, "Reuse inputs").invoke()
    if problem == "file":
        path = tmp_path / "export.csv"
        path.write_text("existing export\n")
        panel.vars["output"].set(str(path))
    else:
        panel.retry.set(False)
    alerts = []
    monkeypatch.setattr("topoppi_affinity.gui.messagebox.showerror", lambda title, text, **_: alerts.append(text))
    panel.run_button.invoke()
    _idle(app.root, panel)
    assert panel.document is previous
    assert panel.worker is worker
    assert panel.table.get_children()
    assert alerts and "folder" in alerts[0].lower()
    if problem == "existing_analysis":
        assert "Keep previous runs in separate directories" in alerts[0]
    panel.vars["output"].set(str(tmp_path / "corrected-result"))
    panel.retry.set(True)
    panel.run_button.invoke()
    _finish(app.root, panel)
    assert panel.document.path.parent.name == "corrected-result"
    assert previous.path.read_bytes() == previous_bytes


def test_cancelled_mutation_run_keeps_error_rows_and_recovers_through_reuse_inputs(desktop, tmp_path, monkeypatch):
    from topoppi_affinity import workflow

    app = desktop
    panel = app.open_binding_task("mutation")
    panel.variants.insert("1.0", "candidate\tA:F2A\nmissing\tA:F999A")
    panel.vars["output"].set(str(tmp_path / "cancelled-result"))
    entered, release = threading.Event(), threading.Event()
    original_prepare = workflow.prepare_interface

    def pause_before_geometry(*args, **kwargs):
        entered.set()
        assert release.wait(10), "The test did not release the paused geometry step"
        return original_prepare(*args, **kwargs)

    monkeypatch.setattr(workflow, "prepare_interface", pause_before_geometry)
    saved_setup = tmp_path / "unfinished.topoppi-setup.json"
    monkeypatch.setattr("topoppi_affinity.gui.filedialog.asksaveasfilename", lambda **_: str(saved_setup))
    panel.run_button.invoke()
    try:
        deadline = time.monotonic() + 10
        while not entered.is_set() and time.monotonic() < deadline:
            app.root.update()
            time.sleep(.01)
        assert entered.is_set(), panel.errors.get()
        assert "disabled" in panel.open_setup_button.state()
        assert "disabled" not in panel.save_setup_button.state()
        panel.save_setup_button.invoke()
        assert saved_setup.is_file()
        panel.cancel_button.invoke()
        assert panel.cancel.is_set()
    finally:
        release.set()
    _idle(app.root, panel)
    assert panel.document is not None, panel.errors.get()
    assert panel.document.predictions.status.tolist() == ["cancelled", "error"]
    assert panel.table.item("candidate", "values")[-1] == "Cancelled"
    assert panel.table.item("missing", "values")[-1] == "Failed"
    cancelled_path = panel.document.path
    _button(panel.window, "Reuse inputs").invoke()
    assert panel.model_source.get() == "Prepared task model"
    panel.variants.delete("1.0", tk.END)
    panel.variants.insert("1.0", "candidate\tA:F2A")
    panel.vars["output"].set(str(tmp_path / "corrected-result"))
    panel.run_button.invoke()
    _finish(app.root, panel)
    assert panel.document.report["model"]["target"] == "ddg"
    assert workflow.load_analysis(cancelled_path).predictions.status.tolist() == ["cancelled", "error"]
