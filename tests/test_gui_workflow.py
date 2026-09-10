import json
import os
import sys
import tempfile
import threading
import tkinter as tk
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from topoppi.config import DEFAULT_RUN_CONFIG
from topoppi.gui_app.app import ProtSurfApp
from topoppi.gui_app.plot_mixin import PlotMixin
from topoppi.gui_app.ui_mixin import UIMixin
from topoppi.gui_app.workflow_mixin import WorkflowMixin, _GuiLogAdapter


class _ValueHolder:
    def __init__(self):
        self.value = ""

    def set(self, value):
        self.value = value


class _WorkflowHarness(WorkflowMixin):
    def __init__(self):
        self.messages = []
        self.var_prolif_path = _ValueHolder()

    def log(self, message):
        self.messages.append(str(message))

    @staticmethod
    def post_to_ui(callback, *args):
        callback(*args)


class _BenchmarkPipelineHarness(WorkflowMixin):
    def __init__(self):
        self._cancel_event = threading.Event()
        self.progress = mock.Mock()
        self.messages = []
        self.errors = []
        self.events = []

    def log(self, message):
        self.messages.append(str(message))

    def set_stage_progress(self, *_args):
        return None

    def _clear_benchmark_outputs(self, _config):
        self.events.append("clear")

    def _finish_task(self):
        self.events.append("finish")

    def show_error(self, message):
        self.errors.append(str(message))

    def finish_cancelled(self, message):
        self.errors.append(str(message))

    @staticmethod
    def post_to_ui(callback, *args):
        callback(*args)


class GUIWorkflowTests(unittest.TestCase):
    def test_gui_log_adapter_routes_debug_details_to_the_module_logger(self):
        messages = []
        adapter = _GuiLogAdapter(messages.append)
        with self.assertLogs("topoppi.gui", level="DEBUG") as captured:
            adapter.debug("prepared %d atoms", 12)

        self.assertEqual(messages, [])
        self.assertIn("prepared 12 atoms", captured.output[0])

    @unittest.skipUnless(
        sys.platform in {"win32", "darwin"} or bool(os.environ.get("DISPLAY")),
        "A graphical display is required for the Tk startup smoke test",
    )
    def test_desktop_app_starts_and_closes(self):
        root = tk.Tk()
        root.withdraw()
        try:
            app = ProtSurfApp(root)
            root.update_idletasks()
            self.assertEqual(str(app.btn_run.cget("text")), "Create Interface Map")
            self.assertTrue(hasattr(app, "btn_redraw"))
            self.assertTrue(hasattr(app, "entry_output_dir_basic"))
            app.close()
        finally:
            try:
                root.destroy()
            except tk.TclError:
                pass

    def test_gui_auto_generation_keeps_the_user_override_empty(self):
        harness = _WorkflowHarness()
        with mock.patch(
            "topoppi.gui_app.workflow_mixin.generate_prolif_interactions",
            return_value="complex.A-B.prolif.json",
        ) as generate:
            result = harness.resolve_prolif_interactions(
                "complex.pdb",
                "A",
                "B",
                source_sha256="input-digest",
                output_dir="results",
            )

        self.assertEqual(result, "complex.A-B.prolif.json")
        generate.assert_called_once_with(
            "complex.pdb",
            "A",
            "B",
            log=mock.ANY,
            source_sha256="input-digest",
            output_dir="results",
        )
        self.assertEqual(harness.var_prolif_path.value, "")

    def test_gui_does_not_silently_replace_failed_prolif_generation(self):
        harness = _WorkflowHarness()
        with mock.patch(
            "topoppi.gui_app.workflow_mixin.generate_prolif_interactions",
            side_effect=RuntimeError("missing interaction stack"),
        ):
            with self.assertRaisesRegex(
                RuntimeError,
                "Could not generate ProLIF annotations",
            ):
                harness.resolve_prolif_interactions("complex.pdb", "A", "B")

    def test_gui_geometric_fallback_is_explicit_opt_in(self):
        harness = _WorkflowHarness()
        with mock.patch(
            "topoppi.gui_app.workflow_mixin.generate_prolif_interactions",
            side_effect=RuntimeError("missing interaction stack"),
        ):
            result = harness.resolve_prolif_interactions(
                "complex.pdb",
                "A",
                "B",
                allow_geometric_fallback=True,
            )

        self.assertIsNone(result)
        self.assertTrue(any("selected geometric fallback" in item for item in harness.messages))

    def test_benchmark_overwrite_clears_after_other_preflight_checks_pass(self):
        harness = _BenchmarkPipelineHarness()
        runner = mock.Mock()
        output_reason = "Benchmark output_root must be empty or contain a matching resume checkpoint."

        def preflight():
            harness.events.append("preflight")
            if harness.events.count("preflight") == 1:
                return {
                    "ready": False,
                    "accepted_job_count": 1,
                    "planned_worker_process_count": 1,
                    "blockers": [output_reason],
                    "output_state": {
                        "acceptable": False,
                        "state": "nonempty_unmatched",
                        "reason": output_reason,
                    },
                }
            return {
                "ready": True,
                "accepted_job_count": 1,
                "planned_worker_process_count": 1,
                "blockers": [],
                "output_state": {"acceptable": True, "state": "empty"},
            }

        runner.preflight.side_effect = preflight
        runner.run.side_effect = lambda: harness.events.append("run") or {"summary": {}}
        params = {"run_mode": "overwrite", "output_root": "results"}

        with (
            mock.patch("topoppi.gui_app.workflow_mixin.BenchmarkRunner", return_value=runner),
            mock.patch("topoppi.gui_app.workflow_mixin.messagebox.showinfo"),
        ):
            harness.run_benchmark_pipeline(params, mock.Mock())

        self.assertEqual(harness.events, ["preflight", "clear", "preflight", "run", "finish"])
        self.assertEqual(harness.errors, [])

    def test_benchmark_overwrite_preserves_outputs_when_another_check_fails(self):
        harness = _BenchmarkPipelineHarness()
        runner = mock.Mock()
        output_reason = "Benchmark output_root must be empty or contain a matching resume checkpoint."
        runner.preflight.return_value = {
            "ready": False,
            "accepted_job_count": 1,
            "planned_worker_process_count": 1,
            "blockers": [output_reason, "OptCuts binary could not be resolved."],
            "output_state": {
                "acceptable": False,
                "state": "nonempty_unmatched",
                "reason": output_reason,
            },
        }

        with mock.patch("topoppi.gui_app.workflow_mixin.BenchmarkRunner", return_value=runner):
            harness.run_benchmark_pipeline(
                {"run_mode": "overwrite", "output_root": "results"},
                mock.Mock(),
            )

        self.assertNotIn("clear", harness.events)
        runner.run.assert_not_called()
        self.assertIn("OptCuts binary could not be resolved", harness.errors[0])


class _PipelineHarness(WorkflowMixin):
    def __init__(self, *, previous_result=False):
        self._cancel_event = threading.Event()
        self.messages = []
        self.errors = []
        self.cancellations = []
        self.progress_updates = []
        self.current_fig = object() if previous_result else None
        self._successful_single_run = {"figure": self.current_fig} if previous_result else None

    def log(self, message):
        self.messages.append(str(message))

    def set_stage_progress(self, *args):
        self.progress_updates.append(args)

    def show_error(self, message):
        self.errors.append(str(message))

    def finish_cancelled(self, message):
        self.cancellations.append(str(message))

    @staticmethod
    def post_to_ui(callback, *args):
        callback(*args)


class _PlotFailureHarness(_PipelineHarness, PlotMixin):
    pass


class SingleRunPreflightTests(unittest.TestCase):
    @staticmethod
    def _params(**changes):
        params = {
            "path": "complex.pdb",
            "chain_a": "A",
            "chain_b": "B",
            "prolif": "",
            "output_dir": "",
            "auto_save": False,
        }
        params.update(changes)
        return params

    def test_missing_chain_is_reported_before_optcuts_hash_or_prolif(self):
        harness = _PipelineHarness(previous_result=True)
        loader = mock.Mock()
        loader.get_protein_chain_ids.return_value = ["A", "B"]

        with (
            mock.patch(
                "topoppi.gui_app.workflow_mixin.PDBLoader",
                return_value=loader,
            ) as loader_class,
            mock.patch("topoppi.gui_app.workflow_mixin.OptCutsUVOptimizer") as optimizer_class,
            mock.patch("topoppi.gui_app.workflow_mixin.sha256_file") as hash_file,
            mock.patch.object(
                harness,
                "resolve_prolif_interactions",
            ) as resolve_prolif,
        ):
            harness.run_pipeline(
                self._params(chain_b="Z"),
                DEFAULT_RUN_CONFIG,
            )

        loader_class.assert_called_once_with("complex.pdb")
        loader.get_chain_data.assert_not_called()
        optimizer_class.assert_not_called()
        hash_file.assert_not_called()
        resolve_prolif.assert_not_called()
        self.assertIn("Selected chain(s) Z were not found", harness.errors[0])
        self.assertIn("Available protein chains: A, B", harness.errors[0])
        self.assertIn("previous successful result", harness.errors[0])

    def test_structure_is_parsed_once_before_optcuts_and_prolif(self):
        harness = _PipelineHarness()
        events = []
        loader = mock.Mock()
        loader.get_protein_chain_ids.side_effect = lambda: (
            events.append("chain list")
            or [
                "A",
                "B",
            ]
        )
        loader.get_chain_data.side_effect = lambda chain: events.append(f"chain {chain}") or ([1.0], [object()])
        optimizer = mock.Mock()
        optimizer.preflight_binary.side_effect = lambda: events.append("OptCuts preflight") or {}

        def build_loader(_path):
            events.append("parse structure")
            return loader

        def build_optimizer(*_args, **_kwargs):
            events.append("create optimizer")
            return optimizer

        def stop_after_prolif(*_args, **_kwargs):
            events.append("generate ProLIF")
            raise RuntimeError("stop after preflight")

        with (
            mock.patch(
                "topoppi.gui_app.workflow_mixin.PDBLoader",
                side_effect=build_loader,
            ) as loader_class,
            mock.patch(
                "topoppi.gui_app.workflow_mixin.OptCutsUVOptimizer",
                side_effect=build_optimizer,
            ),
            mock.patch(
                "topoppi.gui_app.workflow_mixin.sha256_file",
                side_effect=lambda _path: events.append("hash input") or "digest",
            ),
            mock.patch.object(
                harness,
                "resolve_prolif_interactions",
                side_effect=stop_after_prolif,
            ),
        ):
            harness.run_pipeline(self._params(), DEFAULT_RUN_CONFIG)

        loader_class.assert_called_once_with("complex.pdb")
        self.assertEqual(
            events,
            [
                "parse structure",
                "chain list",
                "chain A",
                "chain B",
                "create optimizer",
                "OptCuts preflight",
                "hash input",
                "generate ProLIF",
            ],
        )
        self.assertEqual(harness.errors, ["stop after preflight"])

    def test_auto_save_folder_is_created_before_structure_work(self):
        harness = _PipelineHarness()
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "new" / "figures"
            params = self._params(auto_save=True, output_dir=str(output_dir))

            def stop_after_output_preparation(_path):
                self.assertTrue(output_dir.is_dir())
                raise RuntimeError("stop after output preparation")

            with mock.patch(
                "topoppi.gui_app.workflow_mixin.PDBLoader",
                side_effect=stop_after_output_preparation,
            ) as loader_class:
                harness.run_pipeline(params, DEFAULT_RUN_CONFIG)

            loader_class.assert_called_once_with("complex.pdb")
            self.assertTrue(output_dir.is_dir())
            self.assertEqual(params["output_dir"], str(output_dir))
            self.assertIn("stop after output preparation", harness.errors[0])

    def test_invalid_auto_save_folder_stops_before_structure_work(self):
        harness = _PipelineHarness()
        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "figure-target"
            output_path.write_text("already a file", encoding="utf-8")
            with (
                mock.patch("topoppi.gui_app.workflow_mixin.PDBLoader") as loader_class,
                mock.patch("topoppi.gui_app.workflow_mixin.OptCutsUVOptimizer") as optimizer_class,
            ):
                harness.run_pipeline(
                    self._params(auto_save=True, output_dir=str(output_path)),
                    DEFAULT_RUN_CONFIG,
                )

        loader_class.assert_not_called()
        optimizer_class.assert_not_called()
        self.assertIn("Could not prepare the auto-save folder", harness.errors[0])
        self.assertIn("Choose another output folder", harness.errors[0])

    def test_cancelled_run_labels_a_retained_previous_result(self):
        harness = _PipelineHarness(previous_result=True)
        harness._cancel_event.set()
        harness.run_pipeline(self._params(), DEFAULT_RUN_CONFIG)

        self.assertIn("Run cancelled by user", harness.cancellations[0])
        self.assertIn("previous successful result", harness.cancellations[0])

    def test_render_failure_labels_a_retained_previous_result(self):
        harness = _PlotFailureHarness(previous_result=True)
        viz = mock.Mock()
        viz.plot_patches.side_effect = RuntimeError("render failed")

        harness.update_plot(
            viz,
            [],
            {},
            complete_task=True,
            run_params={},
            run_manifest={},
        )

        self.assertIn("Failed to generate plot: render failed", harness.errors[0])
        self.assertIn("previous successful result", harness.errors[0])


class _ContextHarness(UIMixin):
    def __init__(self):
        interpreter = tk.Tcl()
        self._interpreter = interpreter
        self.var_input_path = tk.StringVar(master=interpreter, value="complex.pdb")
        self.var_chain_a = tk.StringVar(master=interpreter, value="A")
        self.var_chain_b = tk.StringVar(master=interpreter, value="B")
        self.var_prolif_path = tk.StringVar(master=interpreter, value="")
        self._track_prolif_override_context()


class ProLIFOverrideContextTests(unittest.TestCase):
    def test_basic_chain_or_input_change_clears_a_previous_override(self):
        harness = _ContextHarness()

        harness.var_prolif_path.set("complex.A-B.prolif.json")
        harness.var_chain_a.set("B")
        self.assertEqual(harness.var_prolif_path.get(), "")

        harness.var_prolif_path.set("complex.B-A.prolif.json")
        harness.var_input_path.set("other_complex.pdb")
        self.assertEqual(harness.var_prolif_path.get(), "")

    def test_run_readiness_waits_for_a_path_and_two_distinct_chains(self):
        self.assertFalse(UIMixin._required_inputs_ready("single", "", "A", "B"))
        self.assertFalse(UIMixin._required_inputs_ready("single", "complex.pdb", "A", "A"))
        self.assertTrue(UIMixin._required_inputs_ready("single", "complex.pdb", "A", "B"))
        self.assertTrue(UIMixin._required_inputs_ready("benchmark", "structures", "", ""))


class _Figure:
    def __init__(self):
        self.saved = []

    def savefig(self, path, **kwargs):
        self.saved.append((path, kwargs))


class _SaveHarness(UIMixin, PlotMixin):
    def __init__(self, successful_run):
        self._successful_single_run = successful_run
        self.current_fig = successful_run["figure"] if successful_run else None
        self.config = SimpleNamespace(figure_dpi=300)
        self.messages = []
        self.current_run_log = []

    def get_style_config(self):
        return {"preset": "Exploration"}

    def _remember_recent_output_dir(self, _path):
        return None

    def log(self, message):
        self.messages.append(str(message))


class SuccessfulRunSnapshotTests(unittest.TestCase):
    @staticmethod
    def _successful_run(output_dir, figure):
        return {
            "figure": figure,
            "params": {
                "path": "/data/success.pdb",
                "chain_a": "A",
                "chain_b": "B",
                "cutoff": 4.0,
                "res": 1.0,
                "sigma": 1.0,
                "min_points": 3,
                "output_dir": output_dir,
            },
            "manifest": {
                "input_file": "/data/success.pdb",
                "chain_a": "A",
                "chain_b": "B",
                "prolif_source": "generated",
                "run_id": "successful-run",
                "config": {"output_file": "old.png"},
            },
            "log": ["successful run log"],
            "patches": [object()],
            "style": {"preset": "Successful style"},
            "viz": object(),
        }

    def test_save_uses_one_successful_snapshot_after_later_task_state_changes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            figure = _Figure()
            successful_run = self._successful_run(temp_dir, figure)
            harness = _SaveHarness(successful_run)
            harness.current_run_log = ["failed or benchmark run log"]
            saved_path = str(Path(temp_dir) / "saved.png")

            with (
                mock.patch(
                    "topoppi.gui_app.ui_mixin.filedialog.asksaveasfilename",
                    return_value=saved_path,
                ) as choose_path,
                mock.patch("topoppi.gui_app.ui_mixin.messagebox.showinfo"),
                mock.patch("topoppi.gui_app.ui_mixin.messagebox.showerror") as show_error,
            ):
                harness.save_figure()

            expected_default = "success_A-B_cutoff4_res1_sigma1_min3_generated_successful-run.png"
            self.assertEqual(choose_path.call_args.kwargs["initialfile"], expected_default)
            self.assertEqual(figure.saved[0][0], saved_path)
            with open(Path(temp_dir) / "saved.topoppi.json", encoding="utf-8") as handle:
                payload = json.load(handle)
            self.assertEqual(payload["run"]["input_file"], "/data/success.pdb")
            self.assertEqual(payload["run"]["run_id"], "successful-run")
            self.assertEqual(payload["run"]["chain_a"], "A")
            self.assertEqual(payload["run"]["chain_b"], "B")
            self.assertEqual(payload["log"], ["successful run log"])
            self.assertEqual(payload["style"], {"preset": "Successful style"})
            show_error.assert_not_called()

    def test_save_is_unavailable_without_a_successful_single_run(self):
        harness = _SaveHarness(None)
        with mock.patch("topoppi.gui_app.ui_mixin.filedialog.asksaveasfilename") as choose_path:
            harness.save_figure()
        choose_path.assert_not_called()


if __name__ == "__main__":
    unittest.main()
