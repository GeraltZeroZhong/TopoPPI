from __future__ import annotations

import hashlib
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

import numpy as np
import trimesh
from PIL import Image
from scipy.spatial import cKDTree

from topoppi.atlas.footprints import (
    residue_fragmentation_report,
    write_residue_footprint_sidecar,
)
from topoppi.atlas.metrics import UVAtlasMetrics, weighted_stats
from topoppi.atlas.packing import apply_packed_uv, pack_mesh_charts, resolved_chart_gap
from topoppi.atlas.uv import as_corner_uv, set_uv_layout, uv_checksum
from topoppi.config import OptCutsConfig
from topoppi.file_utils import sha256_file
from topoppi.install_optcuts import (
    LINUX_X86_64_SHA256,
    OPTCUTS_AUDITED_UPSTREAM_COMMIT,
    OPTCUTS_UPSTREAM_URL,
    InstallError,
    default_install_dir,
)
from topoppi.mesh.provenance import OPTCUTS_GEOMETRY_VERTEX_IDS

logger = logging.getLogger("UVOOptimizer")
_RESIDUE_AWARE_MARKER = b"residue footprint energy enabled"


@dataclass(frozen=True)
class ParsedOBJUV:
    vertices: np.ndarray
    faces: np.ndarray
    texcoords: np.ndarray
    face_texcoord_indices: np.ndarray
    corner_uv: np.ndarray


def _requested_optcuts_binary(config: OptCutsConfig) -> str:
    configured = str(config.optcuts_bin)
    if configured not in {"OptCuts_bin", "OptCuts_bin.exe"}:
        return configured
    return os.environ.get(config.optcuts_env_var) or configured


def resolve_optcuts_binary(config: OptCutsConfig) -> Optional[str]:
    """Resolve an explicit executable, or discover a compatible installation/build."""

    requested = _requested_optcuts_binary(config)
    bin_path = os.path.expanduser(requested)
    local_bin = os.path.abspath(bin_path)
    resolved_bin = local_bin if os.path.isfile(local_bin) else shutil.which(bin_path)
    # An explicit selection must be corrected by the user, not silently replaced.
    if os.environ.get(config.optcuts_env_var) or requested not in {"OptCuts_bin", "OptCuts_bin.exe"}:
        return os.path.abspath(resolved_bin) if resolved_bin and os.path.isfile(resolved_bin) else None
    filename = "OptCuts_bin.exe" if os.name == "nt" else "OptCuts_bin"
    candidates = [Path(resolved_bin)] if resolved_bin else []
    candidates.append(Path(sys.executable).resolve().parent / filename)
    try:
        candidates.append(default_install_dir() / filename)
    except InstallError:
        pass  # A local source build can still be available on other platforms.
    for start in (Path(__file__).resolve().parent, Path.cwd()):
        candidates.extend(parent / "tools" / "OptCuts" / filename for parent in (start, *start.parents))
    available = list(dict.fromkeys(str(candidate.resolve()) for candidate in candidates if candidate.is_file()))
    if config.residue_fragmentation_weight > 0.0:
        for binary in available:
            if supports_residue_footprint_energy(binary):
                return binary
    return available[0] if available else None


def supports_residue_footprint_energy(binary_path: str) -> bool:
    """Return whether an OptCuts executable advertises the optional energy."""

    with open(binary_path, "rb") as handle:
        return _RESIDUE_AWARE_MARKER in handle.read()


class OptCutsUVOptimizer:
    """OptCuts UV optimizer with explicit initialization and seam-safe output."""

    def __init__(self, config: Optional[OptCutsConfig] = None, cancel_event=None):
        self.config = config or OptCutsConfig()
        self.last_report: Dict[str, object] = {}
        self.cancel_event = cancel_event
        self._binary: tuple[str, str] | None = None

    def _check_cancelled(self) -> None:
        if self.cancel_event is not None and self.cancel_event.is_set():
            raise RuntimeError("OptCuts cancelled by user.")

    def preflight_binary(self) -> Dict[str, str]:
        """Resolve and verify the executable once for this optimizer instance."""

        resolved, digest = self._resolved_binary()
        return {
            "requested": _requested_optcuts_binary(self.config),
            "resolved": resolved,
            "sha256": digest,
            "env_var": self.config.optcuts_env_var,
        }

    def optimize_patches(
        self,
        patches: List[trimesh.Trimesh],
        *,
        initialization: Optional[str] = None,
        pack: bool = True,
        build_report: bool = True,
        source_residue_labels: Sequence[str] | None = None,
        residue_weights: Mapping[str, float] | None = None,
        timeout_sec: float | None = None,
    ) -> List[trimesh.Trimesh]:
        start_ts = time.perf_counter()
        if not patches:
            self.last_report = {"status": "empty_input"}
            return patches

        init_mode = initialization or ("provided" if self.config.use_input_uv else "automatic")
        if init_mode not in {"provided", "automatic"}:
            raise ValueError("OptCuts initialization must be 'provided' or 'automatic'.")
        residue_aware = float(self.config.residue_fragmentation_weight) > 0.0
        if residue_aware and source_residue_labels is None:
            raise ValueError("source_residue_labels are required when residue_fragmentation_weight is positive.")
        effective_timeout = float(self.config.timeout_sec if timeout_sec is None else timeout_sec)
        if not math.isfinite(effective_timeout) or effective_timeout <= 0.0:
            raise ValueError("OptCuts timeout_sec must be finite and positive.")

        for idx, patch in enumerate(patches):
            self._check_cancelled()
            reference_uv = None
            if init_mode == "provided":
                try:
                    reference_uv = as_corner_uv(patch, key="uv")
                except ValueError as exc:
                    raise RuntimeError(f"Patch {idx} is missing initial UV before OptCuts.") from exc
            patch_vertex_count = int(len(patch.vertices))
            patch_face_count = int(len(patch.faces))

            t0 = time.perf_counter()
            opt_uv, execution = self._run_optcuts_for_patch(
                patch,
                reference_uv,
                patch_index=idx,
                source_residue_labels=source_residue_labels,
                residue_weights=residue_weights,
                timeout_sec=effective_timeout,
            )
            self._check_cancelled()
            elapsed_patch = time.perf_counter() - t0
            set_uv_layout(patch, opt_uv, key="uv_optcuts")
            set_uv_layout(patch, opt_uv, key="uv")
            if residue_aware:
                footprint_report = residue_fragmentation_report(
                    [patch],
                    source_residue_labels,
                    uv_key="uv_optcuts",
                    objective_weights=residue_weights,
                )
                execution["residue_aware_objective"].update(
                    {
                        "final_objective_weighted_fragmentation": footprint_report["objective_weighted_fragmentation"],
                        "final_area_weighted_fragmentation": footprint_report["area_weighted_fragmentation"],
                        "final_nonseparating_seam_crossing_edge_count": footprint_report["nonlocality_audit"][
                            "nonseparating_seam_crossing_edge_count"
                        ],
                    }
                )
            patch.metadata["optcuts_runtime_sec"] = float(elapsed_patch)
            patch.metadata["optcuts_execution"] = execution
            patch.metadata["optcuts_initialization"] = init_mode
            logger.info(
                "OptCuts patch %d done in %.3fs (verts=%d, faces=%d, quick=%s, initialization=%s).",
                idx,
                elapsed_patch,
                patch_vertex_count,
                patch_face_count,
                bool(self.config.optcuts_quick_mode),
                init_mode,
            )

        packing_report = {"status": "disabled", "chart_count": len(patches)}
        if pack:
            packed_uv, transforms, packing_report = pack_mesh_charts(
                patches,
                key="uv_optcuts",
                gap=self.config.patch_gap,
            )
            apply_packed_uv(patches, packed_uv, transforms, key="uv_global")
        elif build_report:
            for patch in patches:
                set_uv_layout(patch, as_corner_uv(patch, key="uv_optcuts"), key="uv_global")

        if not build_report:
            self.last_report = {}
            return patches

        self.last_report = self._build_report(
            patches=patches,
            iteration_time=time.perf_counter() - start_ts,
            packing_report=packing_report,
        )
        self.last_report["status"] = "ok"
        self.last_report["treatment"] = (
            "provided_uv_initialized_optcuts" if init_mode == "provided" else "optcuts_automatic_initialization"
        )
        self.last_report["initialization"] = init_mode
        self.last_report["optcuts_runtime"] = {
            "quick_mode": bool(self.config.optcuts_quick_mode),
            "total_patch_count": int(len(patches)),
        }
        self.last_report["residue_aware_objective"] = {
            "enabled": residue_aware,
            "residue_fragmentation_weight": float(self.config.residue_fragmentation_weight),
        }
        for p in patches:
            p.metadata["joint_opt_report"] = self.last_report
        return patches

    def _build_report(
        self,
        patches: List[trimesh.Trimesh],
        iteration_time: float,
        packing_report: Dict[str, object],
    ) -> Dict[str, object]:
        flip_vals = []
        dist_samples = []
        angle_samples = []
        area_samples = []
        patch_weights = []
        seam_stats = []
        for p in patches:
            uv = as_corner_uv(p, key="uv_optcuts")
            patch_weights.append(float(p.area))
            flip_vals.append(UVAtlasMetrics.flip_rate(p, uv))
            dist_samples.append(UVAtlasMetrics.distortion_samples(p, uv))
            angle_samples.append(UVAtlasMetrics.angle_distortion_samples(p, uv))
            area_samples.append(UVAtlasMetrics.area_distortion_samples(p, uv))
            seam_stats.append(UVAtlasMetrics.seam_stats(p, uv))

        def _aggregate_samples(samples):
            values = np.concatenate([np.asarray(item[0], dtype=np.float64) for item in samples])
            weights = np.concatenate([np.asarray(item[1], dtype=np.float64) for item in samples])
            return weighted_stats(values, weights)

        atlas_stats = UVAtlasMetrics.atlas_geometry_stats(
            patches,
            key="uv_global",
            padding=resolved_chart_gap(patches, self.config.patch_gap),
        )
        total_area = float(np.sum(patch_weights))

        return {
            "parameterization_quality": {
                "flip_rate_mean": float(np.average(flip_vals, weights=patch_weights)),
                "distortion": _aggregate_samples(dist_samples),
                "angle_distortion": {**_aggregate_samples(angle_samples), "unit": "radian"},
                "area_distortion": _aggregate_samples(area_samples),
                "aggregation": "original_3d_face_area_weighted_across_all_patches",
            },
            "topology_complexity": {
                "seam_edge_count": int(sum(int(item["seam_edge_count"]) for item in seam_stats)),
                "seam_length_3d": float(sum(float(item["seam_length_3d"]) for item in seam_stats)),
                "seam_length_3d_normalized": float(
                    sum(float(item["seam_length_3d"]) for item in seam_stats) / np.sqrt(total_area)
                ),
                "boundary_length_3d": float(sum(float(item["boundary_length_3d"]) for item in seam_stats)),
                "chart_count": int(len(patches)),
            },
            "atlas_usability": {
                **atlas_stats,
                "packing": packing_report,
            },
            "stability_efficiency": {
                "objective_history": None,
                "objective_drop": None,
                "objective_trace_status": "not_collected_from_optcuts",
                "total_time_sec": float(iteration_time),
                "failure_rate": 0.0,
            },
        }

    def get_last_report(self) -> Dict[str, object]:
        return dict(self.last_report)

    def _run_optcuts_for_patch(
        self,
        patch: trimesh.Trimesh,
        reference_uv: Optional[np.ndarray],
        patch_index: int,
        timeout_sec: float,
        source_residue_labels: Sequence[str] | None = None,
        residue_weights: Mapping[str, float] | None = None,
    ) -> tuple[np.ndarray, Dict[str, object]]:
        effective_timeout = float(timeout_sec)
        bijectivity_enabled = (
            self.config.optcuts_quick_use_bijectivity
            if self.config.optcuts_quick_mode
            else self.config.optcuts_use_bijectivity
        )
        initial_injectivity = None
        source_initial_injectivity = None
        source_initial_uv_checksum = None
        provided_uv_transform = "not_applicable"
        if reference_uv is not None:
            source_initial_uv_checksum = uv_checksum(patch, reference_uv)
            source_initial_injectivity = UVAtlasMetrics.parameterization_injectivity_stats(
                patch,
                reference_uv,
            )
            reference_uv = as_corner_uv(patch, reference_uv).copy()
            if source_initial_injectivity["global_reflection_required_for_positive_orientation"]:
                reference_uv[..., 0] *= -1.0
                provided_uv_transform = "global_u_reflection_for_optcuts_positive_orientation"
                initial_injectivity = UVAtlasMetrics.parameterization_injectivity_stats(
                    patch,
                    reference_uv,
                )
            else:
                provided_uv_transform = "identity"
                initial_injectivity = source_initial_injectivity
            if bijectivity_enabled and not initial_injectivity["globally_injective"]:
                raise RuntimeError(
                    "Provided UV is not globally injective and cannot initialize "
                    "OptCuts with bijectivity enabled "
                    f"(flipped_faces={initial_injectivity['flip_face_count']}, "
                    f"overdraw_ratio={initial_injectivity['overdraw_ratio']:.6g})."
                )
        resolved_bin, binary_sha256 = self._resolved_binary()

        try:
            with tempfile.TemporaryDirectory(prefix="optcuts_") as tmpdir:
                in_obj = os.path.join(tmpdir, "patch_in.obj")
                input_geometry = self._write_obj_with_uv(patch, in_obj, reference_uv)
                fragmentation_weight = float(self.config.residue_fragmentation_weight)
                residue_aware = fragmentation_weight > 0.0
                footprint_metadata: dict[str, object] | None = None
                footprint_sidecar: str | None = None
                if residue_aware:
                    footprint_sidecar = os.path.join(tmpdir, "residue_footprints.txt")
                    initial_uv = (
                        reference_uv
                        if reference_uv is not None
                        else np.zeros((len(patch.vertices), 2), dtype=np.float64)
                    )
                    footprint_metadata = write_residue_footprint_sidecar(
                        patch,
                        initial_uv,
                        source_residue_labels,
                        footprint_sidecar,
                        residue_weights=residue_weights,
                        input_source_vertices=input_geometry["footprint_topology_vertex_ids"],
                    )

                # The bundled binary is invoked via positional parameters (see tools/OptCuts/install_optcuts.sh).
                # Keep the output inside the temporary directory by setting cwd.
                lambda_init = (
                    self.config.optcuts_quick_lambda_init
                    if self.config.optcuts_quick_mode
                    else self.config.optcuts_lambda_init
                )
                distortion_bound = (
                    self.config.optcuts_quick_distortion_bound
                    if self.config.optcuts_quick_mode
                    else self.config.optcuts_distortion_bound
                )
                cmd = [
                    resolved_bin,
                    str(int(self.config.optcuts_mode)),  # mode
                    os.path.basename(in_obj),  # input mesh in tmpdir
                    f"{float(lambda_init):.17g}",  # lambda_init
                    str(int(self.config.optcuts_prog_mode)),  # testID
                    str(int(self.config.optcuts_method_type)),  # methodType
                    f"{float(distortion_bound):.17g}",  # distortionBound
                    str(int(bijectivity_enabled)),  # useBijectivity
                    str(int(self.config.optcuts_initial_cut_option)),
                    self.config.optcuts_output_tag,  # output tag
                ]
                proc_env = os.environ.copy()
                proc_env.pop("TOPOPPI_FOOTPRINT_SIDECAR", None)
                proc_env.pop("TOPOPPI_FRAGMENTATION_WEIGHT", None)
                if residue_aware:
                    proc_env["TOPOPPI_FOOTPRINT_SIDECAR"] = str(footprint_sidecar)
                    proc_env["TOPOPPI_FRAGMENTATION_WEIGHT"] = f"{fragmentation_weight:.17g}"
                if int(self.config.optcuts_mode) >= int(self.config.optcuts_headless_mode):
                    # Headless mode must not inherit an interactive display backend.
                    proc_env.pop("DISPLAY", None)
                    proc_env.pop("WAYLAND_DISPLAY", None)
                proc = subprocess.Popen(
                    cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=tmpdir, env=proc_env
                )
                child_cpu_affinity = (
                    sorted(os.sched_getaffinity(proc.pid)) if hasattr(os, "sched_getaffinity") else None
                )
                process_started = time.perf_counter()
                while True:
                    try:
                        stdout, stderr = proc.communicate(timeout=0.1)
                        break
                    except subprocess.TimeoutExpired:
                        if time.perf_counter() - process_started > effective_timeout:
                            self._stop_process(proc)
                            raise RuntimeError(f"OptCuts timed out after {effective_timeout:.1f}s.") from None
                        if self.cancel_event is not None and self.cancel_event.is_set():
                            self._stop_process(proc)
                            raise RuntimeError("OptCuts cancelled by user.") from None
                if proc.returncode != 0:
                    raise RuntimeError(f"OptCuts failed (code={proc.returncode}): {stderr.strip()}")
                out_obj = self._locate_optcuts_output_obj(tmpdir)
                if not os.path.exists(out_obj):
                    raise RuntimeError(f"OptCuts output OBJ not found: {out_obj}")

                parsed = self._parse_obj_uv(out_obj)
                if parsed is None:
                    raise RuntimeError(f"Failed to parse UV from OptCuts output: {out_obj}")
                uv = self._align_output_corners(patch, parsed)
                output_injectivity = UVAtlasMetrics.parameterization_injectivity_stats(patch, uv)
                if bijectivity_enabled and not output_injectivity["globally_injective"]:
                    raise RuntimeError(
                        "OptCuts returned a UV map that violates the enabled global "
                        "bijectivity constraint "
                        f"(flipped_faces={output_injectivity['flip_face_count']}, "
                        f"overdraw_ratio={output_injectivity['overdraw_ratio']:.6g})."
                    )
                output_constraint_energy = UVAtlasMetrics.optcuts_constraint_energy(
                    patch,
                    uv,
                )
                constraint_tolerance = max(1.0e-6, 1.0e-6 * float(distortion_bound))
                constraint_satisfied = bool(
                    np.isfinite(output_constraint_energy)
                    and output_constraint_energy <= float(distortion_bound) + constraint_tolerance
                )
                if not constraint_satisfied:
                    raise RuntimeError(
                        "OptCuts returned a UV map outside the requested distortion "
                        "constraint "
                        f"({output_constraint_energy:.9g} > {float(distortion_bound):.9g} "
                        f"+ {constraint_tolerance:.3g} tolerance)."
                    )

                self._maybe_export_optcuts_frames(tmpdir=tmpdir, patch_index=patch_index)
                input_source_vertices = np.asarray(
                    input_geometry["footprint_topology_vertex_ids"],
                    dtype=np.int64,
                )
                input_geometry_report = {
                    key: value
                    for key, value in input_geometry.items()
                    if key
                    not in {
                        "source_vertex_ids",
                        "footprint_topology_vertex_ids",
                    }
                }
                input_geometry_report.update(
                    {
                        "footprint_topology_vertex_id_count": int(len(input_source_vertices)),
                        "unique_footprint_topology_vertex_id_count": int(len(np.unique(input_source_vertices))),
                        "footprint_topology_vertex_ids_sha256": hashlib.sha256(
                            np.ascontiguousarray(input_source_vertices).tobytes()
                        ).hexdigest(),
                    }
                )
                execution = {
                    "status": "ok",
                    "initialization": "provided_uv" if reference_uv is not None else "optcuts_automatic",
                    "initial_uv_checksum": uv_checksum(patch, reference_uv) if reference_uv is not None else None,
                    "source_initial_uv_checksum": source_initial_uv_checksum,
                    "provided_uv_transform": provided_uv_transform,
                    "source_initial_uv_injectivity": source_initial_injectivity,
                    "initial_uv_injectivity": initial_injectivity,
                    "input_obj_sha256": sha256_file(in_obj),
                    "input_geometry": input_geometry_report,
                    "output_obj_sha256": sha256_file(out_obj),
                    "output_uv_checksum": uv_checksum(patch, uv),
                    "output_uv_injectivity": output_injectivity,
                    "output_distortion_constraint": {
                        "satisfied": constraint_satisfied,
                        "energy": float(output_constraint_energy),
                        "bound": float(distortion_bound),
                        "numeric_tolerance": float(constraint_tolerance),
                        "identity_value": 4.0,
                        "scale_alignment": "raw_optcuts_output_uv",
                        "aggregation": "original_3d_face_area_weighted_mean",
                        "formula": "||J||_F^2 + ||J^-1||_F^2",
                    },
                    "binary_path": resolved_bin,
                    "binary_sha256": binary_sha256,
                    "cpu_affinity": child_cpu_affinity,
                    "upstream_reference": {
                        "repository": OPTCUTS_UPSTREAM_URL,
                        "audited_commit": OPTCUTS_AUDITED_UPSTREAM_COMMIT,
                        "matches_packaged_linux_binary": binary_sha256 == LINUX_X86_64_SHA256,
                        "note": "Custom binaries require independent source/build provenance.",
                    },
                    "command": [os.path.basename(resolved_bin), *cmd[1:2], "<input.obj>", *cmd[3:]],
                    "lambda_init": float(lambda_init),
                    "distortion_bound": float(distortion_bound),
                    "timeout_sec": effective_timeout,
                    "use_bijectivity": bool(bijectivity_enabled),
                    "initial_cut_option": int(self.config.optcuts_initial_cut_option),
                    "stdout_tail": stdout[-4000:],
                    "stderr_tail": stderr[-4000:],
                    "output_face_count": int(len(parsed.faces)),
                    "output_texture_coordinate_count": int(len(parsed.texcoords)),
                    "per_corner_uv_preserved": True,
                    "residue_aware_objective": {
                        "enabled": residue_aware,
                        "capability_confirmed": bool(residue_aware),
                        "residue_fragmentation_weight": fragmentation_weight,
                        "sidecar_sha256": sha256_file(footprint_sidecar) if footprint_sidecar else None,
                        "sidecar_schema_version": (
                            int(footprint_metadata["schema_version"]) if footprint_metadata else None
                        ),
                        "residue_count": (int(footprint_metadata["residue_count"]) if footprint_metadata else 0),
                        "internal_edge_count": (
                            int(footprint_metadata["internal_edge_count"]) if footprint_metadata else 0
                        ),
                        "initial_seam_edge_count": (
                            int(footprint_metadata["initial_seam_edge_count"]) if footprint_metadata else 0
                        ),
                    },
                }
                return uv, execution
        except (OSError, ValueError) as exc:
            raise RuntimeError(f"OptCuts execution error: {exc}") from exc

    def _resolved_binary(self) -> tuple[str, str]:
        if self._binary is not None:
            return self._binary
        requested = _requested_optcuts_binary(self.config)
        resolved = resolve_optcuts_binary(self.config)
        if not resolved:
            raise RuntimeError(
                f"OptCuts binary not found: {requested}. Run 'topoppi-install-optcuts', "
                f"select a native OptCuts executable in Advanced options, or set {self.config.optcuts_env_var}."
            )
        if os.name != "nt" and not os.access(resolved, os.X_OK):
            raise RuntimeError(f"OptCuts is not executable: {resolved}. Enable executable permission or select another binary in Advanced options.")
        digest = sha256_file(resolved)
        expected = self.config.expected_binary_sha256.strip().lower()
        if expected and digest.lower() != expected:
            raise RuntimeError(f"OptCuts binary checksum mismatch: expected {expected}, got {digest}.")
        if self.config.residue_fragmentation_weight > 0.0 and not supports_residue_footprint_energy(resolved):
            raise RuntimeError(
                f"The selected OptCuts binary ({resolved}) does not support residue-aware layouts. "
                "Select a TopoPPI residue-aware OptCuts build in Advanced options, or build it from tools/OptCuts. "
                "For a geometry-only layout, set the residue fragmentation weight to 0."
            )
        self._binary = (resolved, digest)
        return self._binary

    @staticmethod
    def _stop_process(proc: subprocess.Popen) -> None:
        proc.terminate()
        try:
            proc.communicate(timeout=3)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.communicate(timeout=3)

    def _maybe_export_optcuts_frames(self, tmpdir: str, patch_index: int) -> None:
        if not self.config.save_optcuts_frames:
            return

        output_dir = self.config.optcuts_frames_dir.strip()
        if not output_dir:
            output_dir = os.path.join(os.getcwd(), "optcuts_frames")
        os.makedirs(output_dir, exist_ok=True)

        raster_paths: list[str] = []
        gif_paths: list[str] = []
        for root, _, files in os.walk(tmpdir):
            for name in files:
                lower = name.lower()
                full_path = os.path.join(root, name)
                if lower.endswith((".png", ".bmp", ".jpg", ".jpeg", ".tif", ".tiff")):
                    raster_paths.append(full_path)
                elif lower.endswith(".gif"):
                    gif_paths.append(full_path)

        patch_dir = os.path.join(output_dir, f"patch_{patch_index:03d}")
        os.makedirs(patch_dir, exist_ok=True)

        stride = int(self.config.optcuts_frame_stride)
        min_long_edge = int(self.config.optcuts_min_frame_long_edge)

        def _frame_sort_key(path: str):
            base = os.path.splitext(os.path.basename(path))[0]
            chunks = re.findall(r"\d+", base)
            nums = tuple(int(c) for c in chunks) if chunks else ()
            return (0 if chunks else 1, nums, base, path)

        ordered_rasters = sorted(raster_paths, key=_frame_sort_key)

        def _is_diagnostic(path: str) -> bool:
            base = os.path.basename(path).lower()
            stem = os.path.splitext(base)[0]
            return base == "finalresult.png" or stem.endswith("_distortion") or stem.endswith("_seam")

        def _is_timeline(path: str) -> bool:
            stem = os.path.splitext(os.path.basename(path))[0].lower()
            return not _is_diagnostic(path) and (
                stem.isdigit() or any(tag in stem for tag in ("frame", "viewer", "iter", "step", "anim"))
            )

        timeline = [path for path in ordered_rasters if _is_timeline(path)]
        non_diagnostic = [path for path in ordered_rasters if not _is_diagnostic(path)]
        frame_candidates = timeline if len(timeline) > 1 else non_diagnostic
        if len(frame_candidates) > 1:
            selected = frame_candidates[::stride]
            self._export_raster_frames(selected, patch_dir, min_long_edge)
            logger.info(
                "OptCuts raster frames exported for patch %d: %d image(s) -> %s",
                patch_index,
                len(selected),
                patch_dir,
            )
            return

        extracted_gif_frames = self._export_frames_from_gifs(
            gif_paths=gif_paths,
            patch_dir=patch_dir,
            stride=stride,
            min_long_edge=0,
        )
        if gif_paths and self.config.optcuts_copy_raw_gif:
            for gif_idx, gif_path in enumerate(sorted(gif_paths)):
                gif_name = f"source_{gif_idx:02d}_{os.path.basename(gif_path)}"
                shutil.copy2(gif_path, os.path.join(patch_dir, gif_name))
        if extracted_gif_frames > 0:
            logger.info(
                "OptCuts GIF viewer frames exported for patch %d: %d image(s) -> %s",
                patch_index,
                extracted_gif_frames,
                patch_dir,
            )
            return

        if not ordered_rasters:
            logger.info(
                "OptCuts frame export enabled, but no raster or GIF frames were found for patch %d.", patch_index
            )
            return

        selected = ordered_rasters[::stride]
        self._export_raster_frames(selected, patch_dir, min_long_edge)
        logger.info(
            "OptCuts diagnostic images exported for patch %d: %d image(s) -> %s",
            patch_index,
            len(selected),
            patch_dir,
        )

    @staticmethod
    def _export_raster_frames(paths: List[str], output_dir: str, min_dimension: int) -> None:
        for index, source in enumerate(paths):
            stem = os.path.splitext(os.path.basename(source))[0]
            OptCutsUVOptimizer._copy_png_with_min_resolution(
                src_path=source,
                dst_path=os.path.join(output_dir, f"{index:04d}_{stem}.png"),
                min_long_edge=min_dimension,
            )

    @staticmethod
    def _export_frames_from_gifs(gif_paths: List[str], patch_dir: str, stride: int, min_long_edge: int) -> int:
        if not gif_paths:
            return 0

        def _gif_sort_key(path: str):
            base = os.path.basename(path).lower()
            return (0 if base == "anim.gif" else 1, base, path)

        frame_count = 0
        ordered_gifs = sorted(gif_paths, key=_gif_sort_key)
        for gif_idx, gif_path in enumerate(ordered_gifs):
            with Image.open(gif_path) as gif:
                for frame_idx in range(0, gif.n_frames, stride):
                    gif.seek(frame_idx)
                    frame = gif.convert("RGBA")
                    out_name = f"{frame_count:04d}_viewer_g{gif_idx:02d}_f{frame_idx:05d}.png"
                    OptCutsUVOptimizer._save_png_with_min_resolution(
                        image=frame,
                        dst_path=os.path.join(patch_dir, out_name),
                        min_long_edge=min_long_edge,
                    )
                    frame_count += 1
        return frame_count

    @staticmethod
    def _resize_image_if_needed(image: Image.Image, min_long_edge: int) -> Image.Image:
        if min_long_edge <= 0:
            return image
        w, h = image.size
        scale = max(1.0, float(min_long_edge) / float(max(w, h)))
        if scale <= 1.0:
            return image
        new_size = (int(round(w * scale)), int(round(h * scale)))
        return image.resize(new_size, Image.Resampling.LANCZOS)

    @staticmethod
    def _copy_png_with_min_resolution(src_path: str, dst_path: str, min_long_edge: int) -> None:
        with Image.open(src_path) as img:
            OptCutsUVOptimizer._save_png_with_min_resolution(
                image=img.convert("RGBA"),
                dst_path=dst_path,
                min_long_edge=min_long_edge,
            )

    @staticmethod
    def _save_png_with_min_resolution(image: Image.Image, dst_path: str, min_long_edge: int) -> None:
        out = OptCutsUVOptimizer._resize_image_if_needed(image, min_long_edge=min_long_edge)
        out.save(dst_path)

    @staticmethod
    def _locate_optcuts_output_obj(tmpdir: str) -> str:
        candidate_paths = [
            os.path.join(tmpdir, "output", "finalResult_mesh.obj"),
            os.path.join(tmpdir, "finalResult_mesh.obj"),
        ]
        for path in candidate_paths:
            if os.path.exists(path):
                return path

        for root, _, files in os.walk(tmpdir):
            if "finalResult_mesh.obj" in files:
                return os.path.join(root, "finalResult_mesh.obj")
        return os.path.join(tmpdir, "output", "finalResult_mesh.obj")

    @staticmethod
    def _write_obj_with_uv(
        mesh: trimesh.Trimesh,
        obj_path: str,
        initial_uv: Optional[np.ndarray],
    ) -> dict[str, object]:
        """Write repaired 3-D topology with diskification represented only in UV.

        Parameterization may duplicate vertices solely to turn a manifold patch
        into a disk.  Those copies must share one OBJ geometry vertex so OptCuts
        can price, move, and merge the corresponding seam.  Copies introduced
        earlier to repair disconnected vertex fans retain distinct geometry IDs.
        """

        mesh_vertices = np.asarray(mesh.vertices, dtype=np.float64)
        mesh_faces = np.asarray(mesh.faces, dtype=np.int64)
        source_vertices = np.asarray(
            mesh.metadata.get("source_vertex_ids", np.arange(len(mesh_vertices))),
            dtype=np.int64,
        )
        if source_vertices.shape != (len(mesh_vertices),):
            raise ValueError("source_vertex_ids must contain one ID per mesh vertex.")
        geometry_vertex_ids = np.asarray(
            mesh.metadata.get(
                OPTCUTS_GEOMETRY_VERTEX_IDS,
                np.arange(len(mesh_vertices), dtype=np.int64),
            ),
            dtype=np.int64,
        )
        if geometry_vertex_ids.shape != (len(mesh_vertices),):
            raise ValueError("optcuts_geometry_vertex_ids must contain one ID per mesh vertex.")
        footprint_topology_vertex_ids = (
            geometry_vertex_ids if OPTCUTS_GEOMETRY_VERTEX_IDS in mesh.metadata else source_vertices
        )

        geometry_to_input: dict[int, int] = {}
        representatives: list[int] = []
        input_sources: list[int] = []
        input_footprint_topology_vertices: list[int] = []
        mesh_to_input = np.empty(len(mesh_vertices), dtype=np.int64)
        for mesh_vertex, geometry_vertex in enumerate(geometry_vertex_ids):
            geometry_id = int(geometry_vertex)
            input_vertex = geometry_to_input.get(geometry_id)
            if input_vertex is None:
                input_vertex = len(representatives)
                geometry_to_input[geometry_id] = input_vertex
                representatives.append(mesh_vertex)
                input_sources.append(int(source_vertices[mesh_vertex]))
                input_footprint_topology_vertices.append(int(footprint_topology_vertex_ids[mesh_vertex]))
            elif int(source_vertices[mesh_vertex]) != input_sources[input_vertex]:
                raise ValueError("One OptCuts geometry vertex maps to inconsistent root source vertices.")
            mesh_to_input[mesh_vertex] = input_vertex

        representative_indices = np.asarray(representatives, dtype=np.int64)
        vertices = mesh_vertices[representative_indices]
        faces = mesh_to_input[mesh_faces]
        if np.any(np.diff(np.sort(faces, axis=1), axis=1) == 0):
            raise ValueError("Collapsing diskification copies creates a degenerate OptCuts face.")
        geometry_scale = max(float(np.max(np.ptp(mesh_vertices, axis=0))), 1.0)
        coordinate_tolerance = 1e-8 * geometry_scale
        deviations = np.linalg.norm(mesh_vertices - vertices[mesh_to_input], axis=1)
        if np.any(deviations > coordinate_tolerance):
            raise ValueError(
                "Vertices sharing one OptCuts geometry ID disagree in 3-D position "
                f"({float(np.max(deviations)):.3g} > {coordinate_tolerance:.3g})."
            )
        corners = None if initial_uv is None else as_corner_uv(mesh, initial_uv)
        texcoords: list[tuple[float, float] | None] = [] if corners is None else [None] * len(vertices)
        face_texcoord_indices = np.empty_like(faces)
        if corners is not None:
            texture_index: dict[tuple[int, float, float], int] = {}
            for face_index, face in enumerate(faces):
                for corner_index, input_vertex in enumerate(face):
                    u, v = (float(value) for value in corners[face_index, corner_index])
                    key = (int(input_vertex), u, v)
                    index = texture_index.get(key)
                    if index is None:
                        vertex = int(input_vertex)
                        index = vertex if texcoords[vertex] is None else len(texcoords)
                        texture_index[key] = index
                        if index == vertex:
                            texcoords[vertex] = (u, v)
                        else:
                            texcoords.append((u, v))
                    face_texcoord_indices[face_index, corner_index] = index
        with open(obj_path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write("# TopoPPI seam-preserving OptCuts input\n")
            for x, y, z in vertices:
                handle.write(f"v {x:.17g} {y:.17g} {z:.17g}\n")
            for texcoord in texcoords:
                u, v = texcoord if texcoord is not None else (0.0, 0.0)
                handle.write(f"vt {u:.17g} {v:.17g}\n")
            for face_index, tri in enumerate(faces):
                if corners is None:
                    handle.write("f " + " ".join(str(int(vertex) + 1) for vertex in tri) + "\n")
                else:
                    tokens = [
                        f"{int(vertex) + 1}/{int(face_texcoord_indices[face_index, corner_index]) + 1}"
                        for corner_index, vertex in enumerate(tri)
                    ]
                    handle.write("f " + " ".join(tokens) + "\n")
        return {
            "vertex_count": int(len(vertices)),
            "mesh_vertex_count": int(len(mesh_vertices)),
            "collapsed_vertex_copy_count": int(len(mesh_vertices) - len(vertices)),
            "collapsed_diskification_vertex_copy_count": int(len(mesh_vertices) - len(vertices)),
            "preserved_topology_vertex_copy_count": int(len(input_sources) - len(np.unique(input_sources))),
            "source_vertex_ids": input_sources,
            "footprint_topology_vertex_ids": input_footprint_topology_vertices,
        }

    @staticmethod
    def _parse_obj_uv(obj_path: str) -> Optional[ParsedOBJUV]:
        vertices = []
        texcoords = []
        raw_faces = []
        raw_face_texcoords = []
        with open(obj_path, "r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if line.startswith("v "):
                    parts = line.split()
                    if len(parts) >= 4:
                        vertices.append((float(parts[1]), float(parts[2]), float(parts[3])))
                elif line.startswith("vt "):
                    parts = line.split()
                    if len(parts) >= 3:
                        texcoords.append((float(parts[1]), float(parts[2])))
                elif line.startswith("f "):
                    vertex_indices = []
                    texture_indices = []
                    for token in line.split()[1:]:
                        chunks = token.split("/")
                        if len(chunks) < 2 or not chunks[0] or not chunks[1]:
                            vertex_indices = []
                            break
                        vertex_indices.append(int(chunks[0]))
                        texture_indices.append(int(chunks[1]))
                    if len(vertex_indices) < 3:
                        continue
                    for corner in range(1, len(vertex_indices) - 1):
                        raw_faces.append([vertex_indices[0], vertex_indices[corner], vertex_indices[corner + 1]])
                        raw_face_texcoords.append(
                            [texture_indices[0], texture_indices[corner], texture_indices[corner + 1]]
                        )

        if not vertices or not texcoords or not raw_faces:
            return None
        vertex_count = len(vertices)
        texture_count = len(texcoords)

        def _resolve(raw_index: int, count: int) -> int:
            return count + raw_index if raw_index < 0 else raw_index - 1

        faces = np.asarray(
            [[_resolve(index, vertex_count) for index in face] for face in raw_faces],
            dtype=np.int64,
        )
        face_texcoords = np.asarray(
            [[_resolve(index, texture_count) for index in face] for face in raw_face_texcoords],
            dtype=np.int64,
        )
        if np.any(faces < 0) or np.any(faces >= vertex_count):
            return None
        if np.any(face_texcoords < 0) or np.any(face_texcoords >= texture_count):
            return None
        texture_array = np.asarray(texcoords, dtype=np.float64)
        vertex_array = np.asarray(vertices, dtype=np.float64)
        if not np.isfinite(vertex_array).all() or not np.isfinite(texture_array).all():
            return None
        return ParsedOBJUV(
            vertices=vertex_array,
            faces=faces,
            texcoords=texture_array,
            face_texcoord_indices=face_texcoords,
            corner_uv=texture_array[face_texcoords],
        )

    @staticmethod
    def _align_output_corners(mesh: trimesh.Trimesh, parsed: ParsedOBJUV) -> np.ndarray:
        """Align OptCuts face-corner UV to the unchanged input face order."""

        input_faces = np.asarray(mesh.faces, dtype=np.int64)
        if len(parsed.faces) != len(input_faces):
            raise RuntimeError(
                f"OptCuts changed face count ({len(parsed.faces)} vs {len(input_faces)}); "
                "the output cannot be scored on the frozen domain."
            )
        input_vertices = np.asarray(mesh.vertices, dtype=np.float64)
        geometry_scale = max(float(np.max(np.ptp(input_vertices, axis=0))), 1.0)
        coordinate_tolerance = max(
            1e-12 * geometry_scale,
            256.0 * np.finfo(np.float64).eps * max(float(np.max(np.abs(input_vertices))), 1.0),
        )
        canonical_vertices, input_to_canonical = np.unique(input_vertices, axis=0, return_inverse=True)
        neighbor_count = min(2, len(canonical_vertices))
        distances, output_to_canonical = cKDTree(canonical_vertices).query(
            parsed.vertices,
            k=neighbor_count,
        )
        if neighbor_count == 2:
            if np.any(np.asarray(distances)[:, 1] <= coordinate_tolerance):
                raise RuntimeError("OptCuts output vertex is geometrically ambiguous at the OBJ round-trip tolerance.")
            distances = np.asarray(distances)[:, 0]
            output_to_canonical = np.asarray(output_to_canonical)[:, 0]
        if np.any(distances > coordinate_tolerance):
            raise RuntimeError(
                "OptCuts changed 3-D geometry beyond the OBJ round-trip tolerance "
                f"({float(np.max(distances)):.3g} > {coordinate_tolerance:.3g})."
            )
        mapped_input_faces = np.asarray(input_to_canonical, dtype=np.int64)[input_faces]
        mapped_output_faces = np.asarray(output_to_canonical, dtype=np.int64)[parsed.faces]
        if np.array_equal(mapped_output_faces, mapped_input_faces):
            return np.ascontiguousarray(parsed.corner_uv)

        # A topology cut deliberately duplicates seam vertices at identical 3-D
        # coordinates. Collapse exact input-coordinate copies to stable geometric
        # IDs, then map rounded OBJ output coordinates to those IDs.
        def face_key(ids: np.ndarray) -> tuple[int, int, int]:
            return tuple(sorted(int(value) for value in ids))

        output_face_map: dict[tuple[int, int, int], list[int]] = {}
        for output_index, ids in enumerate(mapped_output_faces):
            key = face_key(ids)
            output_face_map.setdefault(key, []).append(output_index)

        aligned = np.empty((len(input_faces), 3, 2), dtype=np.float64)
        for input_index, input_ids in enumerate(mapped_input_faces):
            candidates = output_face_map.get(face_key(input_ids))
            if not candidates:
                raise RuntimeError(f"OptCuts output face {input_index} does not match the frozen input geometry.")
            output_index = candidates.pop(0)
            output_ids = mapped_output_faces[output_index]
            for input_corner, canonical_id in enumerate(input_ids):
                matches = np.flatnonzero(output_ids == canonical_id)
                if len(matches) != 1:
                    raise RuntimeError("OptCuts output contains a degenerate or ambiguous geometric face.")
                output_corner = int(matches[0])
                aligned[input_index, input_corner] = parsed.corner_uv[output_index, output_corner]
        return aligned
