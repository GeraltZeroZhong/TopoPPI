import logging
import math
import os
import platform
import shutil
import sys
import threading
import time
from dataclasses import asdict, replace
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from tkinter import messagebox

from topoppi import __version__
from topoppi.atlas.footprints import (
    contact_partner_degrees,
    geometric_contact_partner_map,
    residue_aware_residue_weights,
    residue_fragmentation_report,
    source_atom_residue_labels,
)
from topoppi.atlas.uv import set_uv_layout
from topoppi.benchmarking import BenchmarkRunner
from topoppi.config import BenchmarkConfig, TopoPPIRunConfig
from topoppi.errors import ConfigurationError
from topoppi.file_utils import git_worktree_state, sha256_file
from topoppi.interactions.interaction_engine import (
    generate_prolif_interactions,
    load_prolif_partner_map,
)
from topoppi.io.io_loader import PDBLoader
from topoppi.mesh.parameterization import Parameterizer
from topoppi.mesh.surface import SurfaceGenerator
from topoppi.mesh.topology import TopologyManager
from topoppi.optimization.optcuts import OptCutsUVOptimizer
from topoppi.visualization.visualizer import InterfaceVisualizer

logger = logging.getLogger("topoppi.gui")


class PipelineCancelled(Exception):
    """Raised when a GUI run is cancelled cooperatively."""


class _GuiLogAdapter:
    def __init__(self, log_fn):
        self._log_fn = log_fn

    def info(self, message, *args):
        self._emit(message, *args)

    def debug(self, message, *args):
        logger.debug(message, *args)

    def warning(self, message, *args):
        self._emit("Warning: " + str(message), *args)

    def error(self, message, *args):
        self._emit("Error: " + str(message), *args)

    def _emit(self, message, *args):
        if args:
            message = str(message) % args
        self._log_fn(str(message))


class WorkflowMixin:
    def _check_cancelled(self):
        if self._cancel_event.is_set():
            raise PipelineCancelled("Run cancelled by user.")

    @staticmethod
    def _default_optcuts_frame_dir(input_path: str) -> str:
        base_dir = os.path.dirname(input_path) or os.getcwd()
        stem = os.path.splitext(os.path.basename(input_path))[0]
        ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        return os.path.join(base_dir, f"{stem}_optcuts_frames_{ts}")

    def start_benchmark(self):
        try:
            form = self.read_benchmark_form()
        except ConfigurationError as exc:
            messagebox.showerror("Invalid Input", str(exc))
            return

        params = form.to_params()
        if params["formal_mode"] and not messagebox.askyesno(
            "Confirm Formal Benchmark",
            "This runs the full formal benchmark with the selected inputs and output folder. "
            "The run starts after its preflight checks pass.\n\n"
            "Continue?",
        ):
            self.log("Formal benchmark cancelled before preflight; no output was changed.")
            return
        params["run_id"] = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        if params["run_mode"] == "new":
            params["output_root"] = self._timestamped_output_root(params["output_root"], params["run_id"])
        config = form.to_config(output_root=params["output_root"])
        self._remember_recent_file(params["folder"])
        self._remember_recent_output_dir(
            params["output_root"] if os.path.isdir(params["output_root"]) else os.path.dirname(params["output_root"])
        )
        self._update_run_summary()
        self._begin_task("Starting benchmark pipeline...", progress_mode="determinate")
        self._worker_thread = threading.Thread(target=self.run_benchmark_pipeline, args=(params, config), daemon=True)
        self._worker_thread.start()

    def run_benchmark_pipeline(self, params, config: BenchmarkConfig):
        try:
            self.set_stage_progress("Benchmark", 0, "Preparing jobs")
            self._check_cancelled()
            runner = BenchmarkRunner(
                config=config, log_fn=self.log, progress_fn=self._on_benchmark_progress, cancel_event=self._cancel_event
            )
            preflight = runner.preflight()
            output_state = preflight["output_state"]
            if params["run_mode"] == "overwrite" and output_state["state"] == "nonempty_unmatched":
                output_blocker = str(output_state["reason"])
                other_blockers = [blocker for blocker in preflight["blockers"] if blocker != output_blocker]
                if not other_blockers:
                    self._check_cancelled()
                    self._clear_benchmark_outputs(config)
                    preflight = runner.preflight()
            self.log(
                "Benchmark preflight: ready={}, accepted={}, planned isolated processes={}.".format(
                    preflight["ready"],
                    preflight["accepted_job_count"],
                    preflight["planned_worker_process_count"],
                )
            )
            if not preflight["ready"]:
                raise ConfigurationError("Benchmark preflight failed: " + "; ".join(preflight["blockers"]))
            self.log("Benchmark OptCuts runs in headless mode (viewer disabled).")
            output = runner.run()
            self._check_cancelled()
            self.post_to_ui(lambda: self.progress.configure(value=100))
            summary = output.get("summary", {})
            method_distributions = summary.get("method_distributions", {})
            self.log(
                "Benchmark done. attempted={}, valid={}, complete={}, failed={}.".format(
                    int(summary.get("attempted_structure_count", 0)),
                    int(summary.get("valid_structure_count", 0)),
                    int(summary.get("complete_comparison_structure_count", 0)),
                    int(summary.get("failed_structure_count", 0)),
                )
            )
            method_labels = {
                "lscm": "LSCM",
                "harmonic": "harmonic",
                "slim": "SLIM",
                "spherical": "spherical",
                "cylindrical": "cylindrical",
                "optcuts_automatic": "geometry-only OptCuts",
                "optcuts_lscm_initialized": "LSCM-initialized OptCuts",
                "residue_aware_optcuts": "TopoPPI",
            }
            distortion_parts = []
            for method, distribution in method_distributions.items():
                mean = float((distribution.get("distortion_mean") or {}).get("mean", float("nan")))
                if math.isfinite(mean):
                    distortion_parts.append(f"{method_labels.get(method, method)}={mean:.4f}")
            if distortion_parts:
                self.log("Mean distortion: " + ", ".join(distortion_parts))
            runtime = summary.get("isolated_end_to_end_wall_sec", {})
            runtime_median = float(runtime.get("median", float("nan")))
            if math.isfinite(runtime_median):
                self.log(
                    "Isolated end-to-end runtime: median={:.2f}s, p95={:.2f}s; right-censored runs={}.".format(
                        runtime_median,
                        float(runtime.get("p95", float("nan"))),
                        int(summary.get("performance_right_censored_run_count", 0)),
                    )
                )
            self.post_to_ui(
                lambda: messagebox.showinfo(
                    "Benchmark Completed",
                    f"Results saved to:\n{params['output_root']}\n\n"
                    "Generated evidence includes report, summary, frozen manifest, failures, "
                    "per-patch retention, per-face audit sample, worker logs, and checksums.",
                ),
            )
            self.post_to_ui(self._finish_task)
        except PipelineCancelled as exc:
            self.post_to_ui(self.finish_cancelled, str(exc))
        except Exception as e:
            logger.exception("GUI benchmark failed")
            if self._cancel_event.is_set():
                self.post_to_ui(self.finish_cancelled, "Benchmark cancellation requested.")
            else:
                self.post_to_ui(self.show_error, f"Benchmark failed: {e}")

    def _on_benchmark_progress(self, completed: int, total: int, message: str):
        self.post_to_ui(self._set_benchmark_progress_ui, completed, total, message)

    def _set_benchmark_progress_ui(self, completed: int, total: int, message: str):
        total_safe = max(1, int(total))
        completed_safe = max(0, min(int(completed), total_safe))
        percent = int((completed_safe / total_safe) * 100.0)
        self.progress.configure(mode="determinate", maximum=100, value=percent)
        self.stage_status_var.set(f"Benchmark: {percent}% - {message}")
        self.log(f"[Benchmark][Progress] {completed_safe}/{total_safe} ({percent}%) - {message}")

    def resolve_prolif_interactions(
        self,
        pdb_path,
        chain_a,
        chain_b,
        *,
        source_sha256=None,
        output_dir=None,
        allow_geometric_fallback=False,
    ):
        self.log("Generating ProLIF interaction annotations...")
        try:
            output_json = generate_prolif_interactions(
                pdb_path,
                chain_a,
                chain_b,
                log=_GuiLogAdapter(self.log),
                source_sha256=source_sha256,
                output_dir=output_dir,
            )
        except Exception as exc:
            if allow_geometric_fallback:
                self.log(f"Warning: ProLIF generation failed; using the selected geometric fallback. Details: {exc}")
                return None
            raise RuntimeError(
                "Could not generate ProLIF annotations. Check the selected chains and "
                "residue completeness, or choose an existing ProLIF JSON on the Advanced page. "
                f"Details: {exc}"
            ) from exc
        return output_json

    def run_pipeline(self, params, config: TopoPPIRunConfig):
        try:
            self.set_stage_progress("Load", 0, "Preparing inputs")
            self._check_cancelled()
            self._prepare_auto_save_output_dir(params)
            surface_config = config.surface
            topology_config = config.topology
            parameterization_config = config.parameterization
            optcuts_config = config.optcuts
            stage_timings = {}
            run_start = time.perf_counter()

            self.log("Loading protein structure...")
            self.set_stage_progress("Load", 5, "Reading structure")
            t0 = time.perf_counter()
            try:
                loader = PDBLoader(params["path"])
                available_chains = loader.get_protein_chain_ids()
            except Exception as exc:
                raise ValueError(f"Could not read the input structure: {exc}") from exc
            missing_chains = [
                chain_id for chain_id in (params["chain_a"], params["chain_b"]) if chain_id not in available_chains
            ]
            if missing_chains:
                available = ", ".join(available_chains) or "none"
                missing = ", ".join(missing_chains)
                raise ValueError(f"Selected chain(s) {missing} were not found. Available protein chains: {available}.")
            coords_A, atoms_A = loader.get_chain_data(params["chain_a"])
            coords_B, atoms_B = loader.get_chain_data(params["chain_b"])
            if len(coords_A) == 0:
                raise ValueError(f"Chain {params['chain_a']} has no standard protein atoms.")
            if len(coords_B) == 0:
                raise ValueError(f"Chain {params['chain_b']} has no standard protein atoms.")
            stage_timings["load_structure_sec"] = time.perf_counter() - t0
            self._check_cancelled()

            optimizer = OptCutsUVOptimizer(optcuts_config, cancel_event=self._cancel_event)
            optcuts_artifact = optimizer.preflight_binary()
            input_sha256 = sha256_file(params["path"])
            explicit_geometric = config.interaction_source == "geometric"
            provided_prolif = params["prolif"] if not explicit_geometric else ""
            prolif_file = provided_prolif
            if not prolif_file and not explicit_geometric:
                self.log("No ProLIF JSON selected. Generating the required ProLIF annotations automatically.")
                self.set_stage_progress("Load", 8, "Resolving interactions")
                prolif_file = self.resolve_prolif_interactions(
                    params["path"],
                    params["chain_a"],
                    params["chain_b"],
                    source_sha256=input_sha256,
                    output_dir=params["output_dir"] or os.path.dirname(params["path"]) or os.getcwd(),
                    allow_geometric_fallback=(config.visualization.use_geometric_interaction_fallback),
                )
                self._check_cancelled()
            prolif_source = (
                "geometric"
                if explicit_geometric
                else "provided"
                if provided_prolif
                else "generated"
                if prolif_file
                else "geometric_fallback"
            )
            t0 = time.perf_counter()
            if prolif_file:
                interaction_partners = load_prolif_partner_map(
                    prolif_file,
                    atoms_A,
                    atoms_B,
                    expected_chain_a=params["chain_a"],
                    expected_chain_b=params["chain_b"],
                    expected_source_sha256=input_sha256,
                )
                interaction_source = "prolif"
                if not interaction_partners:
                    raise ValueError(
                        "ProLIF did not yield any Chain-A/Chain-B interaction "
                        "residue pairs that resolve against the selected structure."
                    )
                self.log(
                    "ProLIF interactions: {} Chain-A residues and {} residue pairs.".format(
                        len(interaction_partners),
                        sum(len(values) for values in interaction_partners.values()),
                    )
                )
            else:
                interaction_partners = geometric_contact_partner_map(
                    coords_A,
                    atoms_A,
                    coords_B,
                    atoms_B,
                    distance_cutoff=float(params["contact_distance_angstrom"]),
                )
                interaction_source = "geometric" if explicit_geometric else "geometric_fallback"
                self.log(
                    "Geometric contacts: {} Chain-A residues and {} residue pairs at <= {:g} Angstrom.".format(
                        len(interaction_partners),
                        sum(len(values) for values in interaction_partners.values()),
                        float(params["contact_distance_angstrom"]),
                    )
                )
            stage_timings["interaction_assignment_sec"] = time.perf_counter() - t0
            self.set_stage_progress("Load", 20, "Structure loaded")
            self._check_cancelled()

            self.log("Generating molecular surface...")
            self.set_stage_progress("Surface", 25, "Generating molecular surface")
            t0 = time.perf_counter()
            surf_gen = SurfaceGenerator(coords_A, config=surface_config)
            mesh_A = surf_gen.generate_mesh()
            if mesh_A is None:
                raise ValueError(f"Surface generation failed: {surf_gen.last_report}")
            stage_timings["surface_sec"] = time.perf_counter() - t0
            self.set_stage_progress("Surface", 40, "Surface ready")
            self._check_cancelled()

            self.log("Extracting interface patches...")
            self.set_stage_progress("Patch", 44, "Extracting interface patches")
            t0 = time.perf_counter()
            param = Parameterizer(config=parameterization_config)
            if params["save_optcuts_frames"]:
                self.log(
                    "OptCuts frame export enabled: stride={}, min_size={}px, output={}".format(
                        optcuts_config.optcuts_frame_stride,
                        optcuts_config.optcuts_min_frame_long_edge,
                        optcuts_config.optcuts_frames_dir,
                    )
                )
            viz = InterfaceVisualizer(
                chain_A_atoms=atoms_A,
                chain_A_coords=coords_A,
                chain_B_coords=coords_B,
                chain_B_atoms=atoms_B,
                chain_a_id=params["chain_a"],
                chain_b_id=params["chain_b"],
                structure_label=Path(params["path"]).stem,
                prolif_file=prolif_file,
                config=replace(
                    config.visualization,
                    use_geometric_interaction_fallback=(
                        explicit_geometric or config.visualization.use_geometric_interaction_fallback
                    ),
                ),
                interaction_partner_map=interaction_partners,
                contact_distance_angstrom=params["contact_distance_angstrom"],
            )
            viz.interaction_residue_source = interaction_source

            topo = TopologyManager(mesh_A, coords_B, config=topology_config)
            patches = topo.get_interface_patches()
            if not patches:
                raise ValueError(f"No interface found with cutoff {params['cutoff']:.2f}.")
            stage_timings["patch_extraction_sec"] = time.perf_counter() - t0
            self.set_stage_progress("Patch", 50, f"{len(patches)} patch(es) found")
            self._check_cancelled()

            self.log(f"Preparing {len(patches)} patches for OptCuts...")
            self.set_stage_progress("Patch", 54, "Preparing patches")
            t0 = time.perf_counter()
            parameterized_patches = self._parameterize_patches(
                patches,
                param,
                use_input_uv=optcuts_config.use_input_uv,
            )
            if not parameterized_patches:
                raise ValueError("Topology preparation or requested UV initialization failed for all patches.")
            stage_timings["parameterization_sec"] = time.perf_counter() - t0
            self._check_cancelled()

            self.log(f"Running OptCuts on every parameterized interface ({len(parameterized_patches)} patches)...")
            self.set_stage_progress("OptCuts", 64, "Optimizing UVs")
            t0 = time.perf_counter()
            optimized_patches, optimizer_report = self._optimize_patches(
                parameterized_patches,
                optimizer,
                atoms_A,
                contact_distance_angstrom=params["contact_distance_angstrom"],
                interaction_partner_map=interaction_partners,
                interaction_source=interaction_source,
            )
            stage_timings["optcuts_sec"] = time.perf_counter() - t0
            self.set_stage_progress("OptCuts", 85, "Optimization complete")
            self._check_cancelled()

            selected_patches, valid_count, invalid_count = self._split_interfaces_by_interaction_count(
                optimized_patches,
                viz,
                params["min_points"],
            )
            if config.visualization.map_style == "footprints" or config.visualization.view == "surface":
                selected_patches = optimized_patches
            self.log(
                f"Display filter summary (minimum interaction residues = {params['min_points']}): "
                f"above threshold={valid_count}, below threshold={invalid_count}; "
                f"displayed={len(selected_patches)}; all were optimized."
            )
            self.log("Rendering visualization...")
            self.set_stage_progress("Render", 92, "Preparing figure")
            stage_timings["total_pipeline_sec"] = time.perf_counter() - run_start
            manifest = {
                "schema_name": "topoppi_gui_run",
                "schema_version": "2.1",
                "topoppi_version": __version__,
                "run_id": params["run_id"],
                "created_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
                "input_file": os.path.abspath(params["path"]),
                "input_sha256": input_sha256,
                "output_dir": os.path.abspath(params["output_dir"] or os.path.dirname(params["path"]) or os.getcwd()),
                "chain_a": params["chain_a"],
                "chain_b": params["chain_b"],
                "prolif_file": os.path.abspath(prolif_file) if prolif_file else None,
                "prolif_sha256": sha256_file(prolif_file) if prolif_file else None,
                "prolif_source": prolif_source,
                "environment": self._environment_block(),
                "git_commit": self._git_commit(),
                "optcuts_resolved": optcuts_artifact,
                "config": config.to_dict(),
                "surface": asdict(surface_config),
                "surface_generation": dict(surf_gen.last_report),
                "topology": asdict(topology_config),
                "topology_extraction": dict(topo.last_report),
                "parameterization": asdict(parameterization_config),
                "optcuts": asdict(optcuts_config),
                "interaction_summary": {
                    "source": interaction_source,
                    "definition": (
                        "Chain A/B residue pairs present in ProLIF interaction records"
                        if interaction_source == "prolif"
                        else (
                            "Chain A/B heavy-atom residue pair distance <= "
                            f"{float(params['contact_distance_angstrom']):g} Angstrom "
                            f"({interaction_source.replace('_', ' ')})"
                        )
                    ),
                    "chain_a_interaction_residue_count": len(interaction_partners),
                    "interaction_residue_pair_count": sum(len(values) for values in interaction_partners.values()),
                },
                "visualization": asdict(config.visualization),
                "stage_timings": stage_timings,
                "optimizer_report": optimizer_report,
                "interface_counts": {
                    "raw_patches": len(patches),
                    "parameterized_patches": len(parameterized_patches),
                    "optimized_patches": len(optimized_patches),
                    "valid_patches": valid_count,
                    "invalid_patches": invalid_count,
                    "displayed_patches": len(selected_patches),
                    "footprint_display_uses_complete_domain": config.visualization.map_style == "footprints",
                    "surface_display_uses_complete_domain": config.visualization.view == "surface",
                },
            }
            self.set_stage_progress("Render", 100, "Rendering")
            self.post_to_ui(
                self.accept_pipeline_result,
                viz,
                selected_patches,
                manifest,
                params,
                optimized_patches,
            )
        except PipelineCancelled as exc:
            self.post_to_ui(self.finish_cancelled, self._previous_result_message(str(exc)))
        except Exception as e:
            logger.exception("GUI single-run pipeline failed")
            if self._cancel_event.is_set():
                self.post_to_ui(
                    self.finish_cancelled,
                    self._previous_result_message("Run cancelled by user."),
                )
            else:
                self.post_to_ui(self.show_error, self._previous_result_message(str(e)))

    @staticmethod
    def _prepare_auto_save_output_dir(params):
        if not params["auto_save"]:
            return
        output_dir = Path(params["output_dir"] or os.path.dirname(params["path"]) or os.getcwd()).expanduser()
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise ValueError(
                f"Could not prepare the auto-save folder '{output_dir}': {exc}. "
                "Choose another output folder or create it first."
            ) from exc
        params["output_dir"] = str(output_dir)

    def _previous_result_message(self, message):
        successful_run = self._successful_single_run
        if successful_run is not None and successful_run["figure"] is self.current_fig:
            return f"{message}\n\nThe canvas still shows the previous successful result."
        return message

    def accept_pipeline_result(self, viz, patches, manifest, params, all_patches=None):
        self._pending_single_run = {
            "viz": viz,
            "patches": list(patches),
            "all_patches": list(patches if all_patches is None else all_patches),
            "manifest": dict(manifest),
            "params": dict(params),
            "style": dict(self._run_style or {}),
        }
        if self._run_style is not None:
            self._restore_atlas_style(self._run_style)
        self.var_annotation_target.set("Current map")
        self._render_pending_result(style=self._run_style)

    def _parameterize_patches(self, patches, parameterizer, *, use_input_uv):
        valid_patches = []
        for p in patches:
            self._check_cancelled()
            if parameterizer.prepare_patch(p) is None:
                continue
            if use_input_uv:
                uv = parameterizer.flatten_patch(p)
                if uv is None:
                    continue
                set_uv_layout(p, uv, key="uv")
            valid_patches.append(p)
        return valid_patches

    def _optimize_patches(
        self,
        patches,
        optimizer,
        atoms_a,
        *,
        contact_distance_angstrom,
        interaction_partner_map,
        interaction_source,
    ):
        self._check_cancelled()
        interaction_counts = contact_partner_degrees(interaction_partner_map)
        source_labels = source_atom_residue_labels(atoms_a)
        objective_weights = residue_aware_residue_weights(source_labels, interaction_counts)
        result = optimizer.optimize_patches(
            patches,
            source_residue_labels=source_labels,
            residue_weights=objective_weights,
        )
        self._check_cancelled()
        report = optimizer.get_last_report()
        weight_definition = (
            "distinct Chain-B residues paired with each Chain-A residue in ProLIF records"
            if interaction_source == "prolif"
            else (
                "distinct Chain-B residues with any heavy-atom pair at distance <= "
                f"{float(contact_distance_angstrom):g} Angstrom "
                f"({interaction_source.replace('_', ' ')})"
            )
        )
        report["residue_footprint_fragmentation"] = {
            "interaction_residue_source": interaction_source,
            "interaction_weight_definition": weight_definition,
            **residue_fragmentation_report(
                result,
                source_labels,
                uv_key="uv_optcuts",
                interaction_weights=interaction_counts,
                objective_weights=objective_weights,
            ),
        }
        self._log_joint_report(optimizer)
        return result, report

    def _log_joint_report(self, optimizer):
        report = optimizer.get_last_report()
        pq = report["parameterization_quality"]
        tc = report["topology_complexity"]
        au = report["atlas_usability"]
        se = report["stability_efficiency"]
        self.log(
            "[JointReport] flip={:.4f}, dist(mean/max/p95)=({:.4f}/{:.4f}/{:.4f})".format(
                float(pq["flip_rate_mean"]),
                float(pq["distortion"]["mean"]),
                float(pq["distortion"]["max"]),
                float(pq["distortion"]["p95"]),
            )
        )
        self.log(
            "[JointReport] internal_seam_3d={:.3f}, charts={}, polygonal_overlap={:.4f}, padding_viol={}, util={:.4f}".format(
                float(tc["seam_length_3d"]),
                int(tc["chart_count"]),
                float(au["overlap_area"]),
                int(au["padding_violations"]),
                float(au["utilization"]),
            )
        )
        objective_drop = se["objective_drop"]
        objective_drop_text = "n/a" if objective_drop is None else f"{float(objective_drop):.4f}"
        self.log(
            "[JointReport] obj_drop={}, total_time={:.3f}s, failure_rate={:.3f}".format(
                objective_drop_text,
                float(se["total_time_sec"]),
                float(se["failure_rate"]),
            )
        )

    def _split_interfaces_by_interaction_count(self, patches, viz, min_points):
        valid = []
        invalid = 0
        for p in patches:
            self._check_cancelled()
            interaction_residue_count = viz.count_patch_interaction_residues(p)
            p.metadata["interaction_residue_count"] = interaction_residue_count
            if interaction_residue_count >= min_points:
                valid.append(p)
            else:
                invalid += 1
        return valid, len(valid), invalid

    def _timestamped_output_root(self, output_root: str, run_id: str) -> str:
        base = Path(output_root).expanduser()
        if base.name:
            return str(base.with_name(f"{base.name}_{run_id}"))
        return str(base / f"benchmark_results_{run_id}")

    def _clear_benchmark_outputs(self, config: BenchmarkConfig) -> None:
        os.makedirs(config.output_root, exist_ok=True)
        for filename in (
            config.checkpoint_filename,
            config.report_filename,
            config.summary_filename,
            config.manifest_filename,
            config.failures_filename,
            config.per_patch_filename,
            config.per_face_sample_filename,
            config.per_residue_filename,
            config.provenance_filename,
            config.optcuts_execution_filename,
            config.artifact_checksums_filename,
        ):
            path = os.path.join(config.output_root, filename)
            if os.path.exists(path):
                os.remove(path)
                self.log(f"Removed previous benchmark output: {path}")
        worker_logs = os.path.join(config.output_root, config.worker_log_folder)
        if os.path.isdir(worker_logs):
            shutil.rmtree(worker_logs)
            self.log(f"Removed previous benchmark worker logs: {worker_logs}")

    def _git_commit(self) -> str | None:
        repo_root = Path(__file__).resolve().parents[3]
        commit, _dirty = git_worktree_state(repo_root)
        return commit

    def _environment_block(self) -> dict[str, object]:
        names = ["numpy", "scipy", "biopython", "matplotlib", "trimesh", "scikit-image", "topoppi"]
        versions = {}
        for name in names:
            try:
                versions[name] = metadata.version(name)
            except metadata.PackageNotFoundError:
                versions[name] = None
        return {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "packages": versions,
        }
