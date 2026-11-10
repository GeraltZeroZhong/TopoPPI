"""Pure GUI form parsing helpers.

These functions intentionally do not import Tkinter so they can be unit-tested
without a display server.
"""

from __future__ import annotations

import os
from dataclasses import asdict, dataclass, replace
from typing import Mapping

from topoppi.config import (
    DEFAULT_BENCHMARK_CONFIG,
    DEFAULT_GUI_CONFIG,
    DEFAULT_RUN_CONFIG,
    BenchmarkConfig,
    TopoPPIRunConfig,
)
from topoppi.errors import ConfigurationError

DEFAULT_GUI_BENCHMARK_VARIANTS = (
    "optcuts_automatic",
    "optcuts_lscm_initialized",
    "residue_aware_optcuts",
)


@dataclass(frozen=True)
class SingleRunForm:
    path: str
    chain_a: str
    chain_b: str
    prolif: str
    interaction_source: str
    map_style: str
    cutoff: float
    contact_distance_angstrom: float
    res: float
    sigma: float
    surface_level: float
    surface_padding: float
    max_voxels: int
    adaptive_resolution: bool
    max_adaptive_resolution: float
    parameterization_method: str
    slim_iterations: int
    slim_boundary_constraint_weight: float
    min_points: int
    optcuts_bin: str
    expected_optcuts_sha256: str
    patch_gap: float
    optcuts_lambda: float
    optcuts_distortion_bound: float
    optcuts_initial_cut_option: int
    optcuts_use_bijectivity: bool
    optcuts_initialization: str
    optcuts_timeout: float
    residue_fragmentation_weight: float
    save_optcuts_frames: bool
    optcuts_frame_stride: int
    optcuts_min_frame_long_edge: int
    optcuts_frames_dir: str
    output_dir: str
    auto_save: bool

    def to_params(self) -> dict[str, object]:
        return asdict(self)

    def to_config(self, *, output_file: str | None = None) -> TopoPPIRunConfig:
        """Build the exact runtime configuration represented by this form."""

        target_output = output_file or os.path.join(
            self.output_dir or os.path.dirname(self.path) or os.getcwd(),
            "interface_map.png",
        )
        config = TopoPPIRunConfig(
            pdb_file=self.path,
            chain_a=self.chain_a,
            chain_b=self.chain_b,
            output_file=target_output,
            prolif_file=self.prolif or None,
            interaction_source=self.interaction_source,
            contact_distance_angstrom=self.contact_distance_angstrom,
            surface=_surface_config(self),
            topology=replace(DEFAULT_RUN_CONFIG.topology, distance_cutoff=self.cutoff),
            parameterization=_parameterization_config(self),
            optcuts=replace(
                _optcuts_config(self),
                use_input_uv=self.optcuts_initialization == "provided",
                save_optcuts_frames=self.save_optcuts_frames,
                optcuts_frame_stride=self.optcuts_frame_stride,
                optcuts_min_frame_long_edge=self.optcuts_min_frame_long_edge,
                optcuts_frames_dir=self.optcuts_frames_dir,
            ).for_headless(),
            visualization=replace(
                DEFAULT_RUN_CONFIG.visualization,
                min_points=self.min_points,
                map_style=self.map_style,
            ),
        )
        return config


@dataclass(frozen=True)
class BenchmarkForm:
    folder: str
    chain_a: str
    chain_b: str
    cutoff: float
    contact_distance_angstrom: float
    res: float
    sigma: float
    surface_level: float
    surface_padding: float
    max_voxels: int
    adaptive_resolution: bool
    max_adaptive_resolution: float
    parameterization_method: str
    slim_iterations: int
    slim_boundary_constraint_weight: float
    patch_gap: float
    optcuts_bin: str
    expected_optcuts_sha256: str
    optcuts_lambda: float
    optcuts_distortion_bound: float
    optcuts_initial_cut_option: int
    optcuts_use_bijectivity: bool
    optcuts_timeout: float
    residue_fragmentation_weight: float
    output_root: str
    resume: bool
    run_mode: str
    chain_selection_mode: str
    manifest_path: str
    repetitions: int
    warmup_runs: int
    formal_mode: bool
    max_workers: int | None
    benchmark_purpose: str
    execution_profile: str
    optcuts_variants: tuple[str, ...]
    include_topology_ablation: bool
    threads_per_worker: int
    worker_timeout_sec: float
    worker_memory_limit_mb: float | None
    raster_size: int
    min_chain_residues: int
    per_face_sample_size_per_patch: int
    bootstrap_iterations: int
    random_seed: int
    expected_git_commit: str
    coordinate_audit_path: str
    expected_coordinate_audit_sha256: str

    def to_params(self) -> dict[str, object]:
        return asdict(self)

    def to_config(self, *, output_root: str | None = None) -> BenchmarkConfig:
        """Build the benchmark runner configuration represented by this form."""

        config = BenchmarkConfig(
            input_folder=self.folder,
            output_root=output_root or self.output_root,
            chain_a=self.chain_a,
            chain_b=self.chain_b,
            chain_selection_mode=self.chain_selection_mode,
            manifest_path=self.manifest_path,
            surface=_surface_config(self),
            topology=replace(DEFAULT_RUN_CONFIG.topology, distance_cutoff=self.cutoff),
            parameterization=_parameterization_config(self),
            optcuts=_optcuts_config(self).for_headless(),
            raster_size=self.raster_size,
            max_workers=self.max_workers,
            repetitions=self.repetitions,
            warmup_runs=self.warmup_runs,
            formal_mode=self.formal_mode,
            expected_git_commit=self.expected_git_commit,
            coordinate_audit_path=self.coordinate_audit_path,
            expected_coordinate_audit_sha256=self.expected_coordinate_audit_sha256,
            benchmark_purpose=self.benchmark_purpose,
            execution_profile=self.execution_profile,
            optcuts_variants=self.optcuts_variants,
            include_topology_ablation=self.include_topology_ablation,
            random_seed=self.random_seed,
            bootstrap_iterations=self.bootstrap_iterations,
            threads_per_worker=self.threads_per_worker,
            contact_distance_angstrom=self.contact_distance_angstrom,
            per_face_sample_size_per_patch=self.per_face_sample_size_per_patch,
            worker_timeout_sec=self.worker_timeout_sec,
            worker_memory_limit_mb=self.worker_memory_limit_mb,
            min_chain_residues=self.min_chain_residues,
            show_tqdm=False,
            resume=self.resume,
        )
        return config


def _surface_config(form: SingleRunForm | BenchmarkForm):
    return replace(
        DEFAULT_RUN_CONFIG.surface,
        grid_resolution=form.res,
        sigma=form.sigma,
        level=form.surface_level,
        padding=form.surface_padding,
        max_voxels=form.max_voxels,
        adaptive_resolution=form.adaptive_resolution,
        max_adaptive_resolution=form.max_adaptive_resolution,
    )


def _parameterization_config(form: SingleRunForm | BenchmarkForm):
    return replace(
        DEFAULT_RUN_CONFIG.parameterization,
        method=form.parameterization_method,
        slim_iterations=form.slim_iterations,
        slim_boundary_constraint_weight=form.slim_boundary_constraint_weight,
    )


def _optcuts_config(form: SingleRunForm | BenchmarkForm):
    return replace(
        DEFAULT_RUN_CONFIG.optcuts,
        optcuts_bin=form.optcuts_bin,
        expected_binary_sha256=form.expected_optcuts_sha256,
        patch_gap=form.patch_gap,
        optcuts_lambda_init=form.optcuts_lambda,
        optcuts_distortion_bound=form.optcuts_distortion_bound,
        optcuts_initial_cut_option=form.optcuts_initial_cut_option,
        optcuts_use_bijectivity=form.optcuts_use_bijectivity,
        timeout_sec=form.optcuts_timeout,
        residue_fragmentation_weight=form.residue_fragmentation_weight,
    )


def parse_single_run_form(raw: Mapping[str, object]) -> SingleRunForm:
    path = _optional_text(raw.get("path"))

    prolif = _optional_text(raw.get("prolif"))
    output_dir = _optional_output_root(raw.get("output_dir"), "Output directory")

    chain_a = _optional_text(raw.get("chain_a"))
    chain_b = _optional_text(raw.get("chain_b"))

    form = SingleRunForm(
        path=path,
        chain_a=chain_a,
        chain_b=chain_b,
        prolif=prolif,
        interaction_source=_optional_text(raw.get("interaction_source", "prolif")),
        map_style=_choice(raw.get("map_style", "markers"), "Map style", {"markers", "footprints"}),
        cutoff=_parse_float(
            raw.get("cutoff", DEFAULT_GUI_CONFIG.default_patch_cutoff),
            "Interface cutoff",
        ),
        contact_distance_angstrom=_parse_float(
            raw.get("contact_distance_angstrom", DEFAULT_RUN_CONFIG.contact_distance_angstrom),
            "Geometric fallback distance",
        ),
        res=_parse_float(raw.get("res"), "Grid resolution"),
        sigma=_parse_float(raw.get("sigma"), "Surface sigma"),
        surface_level=_parse_float(
            raw.get("surface_level", DEFAULT_RUN_CONFIG.surface.level),
            "Surface level",
        ),
        surface_padding=_parse_float(
            raw.get("surface_padding", DEFAULT_RUN_CONFIG.surface.padding),
            "Surface padding",
        ),
        max_voxels=_parse_int(
            raw.get("max_voxels"),
            "Maximum voxels",
            default=DEFAULT_RUN_CONFIG.surface.max_voxels,
        ),
        adaptive_resolution=_boolean(raw.get("adaptive_resolution"), DEFAULT_RUN_CONFIG.surface.adaptive_resolution),
        max_adaptive_resolution=_parse_float(
            raw.get("max_adaptive_resolution", DEFAULT_RUN_CONFIG.surface.max_adaptive_resolution),
            "Maximum adaptive resolution",
        ),
        parameterization_method=_choice(
            raw.get("parameterization_method", DEFAULT_RUN_CONFIG.parameterization.method),
            "Parameterization method",
            {"auto", "lscm", "harmonic", "slim", "spherical", "cylindrical"},
        ),
        slim_iterations=_parse_int(
            raw.get("slim_iterations"),
            "SLIM iterations",
            default=DEFAULT_RUN_CONFIG.parameterization.slim_iterations,
        ),
        slim_boundary_constraint_weight=_parse_float(
            raw.get(
                "slim_boundary_constraint_weight",
                DEFAULT_RUN_CONFIG.parameterization.slim_boundary_constraint_weight,
            ),
            "SLIM boundary constraint weight",
        ),
        min_points=_parse_int(
            raw.get("min_points"),
            "Minimum interaction residues",
            default=DEFAULT_GUI_CONFIG.default_min_points,
        ),
        optcuts_bin=_text_or_default(raw.get("optcuts_bin"), DEFAULT_RUN_CONFIG.optcuts.optcuts_bin),
        expected_optcuts_sha256=_optional_text(raw.get("expected_optcuts_sha256")).lower(),
        patch_gap=_parse_float(raw.get("patch_gap", DEFAULT_RUN_CONFIG.optcuts.patch_gap), "Patch gap"),
        optcuts_lambda=_parse_float(
            raw.get("optcuts_lambda", DEFAULT_RUN_CONFIG.optcuts.optcuts_lambda_init),
            "OptCuts lambda",
        ),
        optcuts_distortion_bound=_parse_float(
            raw.get("optcuts_distortion_bound", DEFAULT_RUN_CONFIG.optcuts.optcuts_distortion_bound),
            "OptCuts distortion bound",
        ),
        optcuts_initial_cut_option=_parse_int(
            raw.get("optcuts_initial_cut_option"),
            "OptCuts initial cut option",
            default=DEFAULT_RUN_CONFIG.optcuts.optcuts_initial_cut_option,
        ),
        optcuts_use_bijectivity=_boolean(
            raw.get("optcuts_use_bijectivity"),
            DEFAULT_RUN_CONFIG.optcuts.optcuts_use_bijectivity,
        ),
        optcuts_initialization=_choice(
            raw.get(
                "optcuts_initialization",
                "provided" if DEFAULT_RUN_CONFIG.optcuts.use_input_uv else "automatic",
            ),
            "OptCuts initialization",
            {"provided", "automatic"},
        ),
        optcuts_timeout=_parse_float(
            raw.get("optcuts_timeout", DEFAULT_RUN_CONFIG.optcuts.timeout_sec),
            "OptCuts timeout",
        ),
        residue_fragmentation_weight=_parse_float(
            raw.get(
                "residue_fragmentation_weight",
                DEFAULT_RUN_CONFIG.optcuts.residue_fragmentation_weight,
            ),
            "TopoPPI objective weight",
        ),
        save_optcuts_frames=_boolean(raw.get("save_optcuts_frames"), False),
        optcuts_frame_stride=_parse_int(
            raw.get("optcuts_frame_stride"),
            "OptCuts frame stride",
            default=DEFAULT_RUN_CONFIG.optcuts.optcuts_frame_stride,
        ),
        optcuts_min_frame_long_edge=_parse_int(
            raw.get("optcuts_min_frame_long_edge"),
            "OptCuts minimum frame size",
            default=DEFAULT_RUN_CONFIG.optcuts.optcuts_min_frame_long_edge,
        ),
        optcuts_frames_dir=str(raw.get("optcuts_frames_dir") or "").strip(),
        output_dir=output_dir,
        auto_save=_boolean(raw.get("auto_save"), DEFAULT_GUI_CONFIG.auto_save_single_run),
    )
    form.to_config().validate()
    return form


def parse_benchmark_form(raw: Mapping[str, object]) -> BenchmarkForm:
    folder = _optional_text(raw.get("folder"))

    output_root = _optional_output_root(raw.get("output_root"), "Benchmark output folder")
    if not output_root:
        output_root = os.path.join(folder, DEFAULT_GUI_CONFIG.benchmark_output_folder)

    selection_mode = _choice(
        raw.get("chain_selection_mode", "configured"),
        "Chain selection mode",
        {"configured", "auto_contact", "manifest"},
    )
    chain_a = _text_or_default(raw.get("chain_a"), DEFAULT_RUN_CONFIG.chain_a)
    chain_b = _text_or_default(raw.get("chain_b"), DEFAULT_RUN_CONFIG.chain_b)
    manifest_path = _optional_text(raw.get("manifest_path"))

    run_mode = str(raw.get("run_mode") or "resume").strip().lower()
    if run_mode not in {"resume", "new", "overwrite"}:
        raise ConfigurationError("Benchmark run mode must be resume, new, or overwrite.")

    formal_mode = _boolean(raw.get("formal_mode"), False)
    repetitions = _parse_int(raw.get("repetitions"), "Measured repetitions", default=3)
    warmup_runs = _parse_int(raw.get("warmup_runs"), "Warm-up runs", default=0)
    max_workers = _optional_int(raw.get("max_workers"), "Benchmark workers")
    expected_sha = _optional_text(raw.get("expected_optcuts_sha256")).lower()
    benchmark_purpose = _choice(
        raw.get("benchmark_purpose", DEFAULT_BENCHMARK_CONFIG.benchmark_purpose),
        "Benchmark purpose",
        {"quality", "performance"},
    )
    execution_profile = _choice(
        raw.get("execution_profile", DEFAULT_BENCHMARK_CONFIG.execution_profile),
        "Execution profile",
        {"comparative", "operational_optcuts"},
    )
    residue_fragmentation_weight = _parse_float(
        raw.get(
            "residue_fragmentation_weight",
            DEFAULT_RUN_CONFIG.optcuts.residue_fragmentation_weight,
        ),
        "TopoPPI objective weight",
    )
    optcuts_variants = _method_tuple(
        raw.get("optcuts_variants"),
        default=(
            DEFAULT_GUI_BENCHMARK_VARIANTS
            if residue_fragmentation_weight > 0.0
            else tuple(DEFAULT_BENCHMARK_CONFIG.resolved_optcuts_variants())
        ),
    )
    coordinate_audit_path = _optional_text(raw.get("coordinate_audit_path"))
    expected_coordinate_audit_sha256 = _optional_text(raw.get("expected_coordinate_audit_sha256")).lower()
    if formal_mode and run_mode == "overwrite":
        raise ConfigurationError(
            "Formal benchmark cannot use overwrite mode; use a new empty output folder "
            "or resume from a matching checkpoint."
        )

    form = BenchmarkForm(
        folder=folder,
        chain_a=chain_a,
        chain_b=chain_b,
        cutoff=_parse_float(
            raw.get("cutoff", DEFAULT_GUI_CONFIG.default_patch_cutoff),
            "Interface cutoff",
        ),
        contact_distance_angstrom=_parse_float(
            raw.get("contact_distance_angstrom", DEFAULT_RUN_CONFIG.contact_distance_angstrom),
            "Geometric contact distance",
        ),
        res=_parse_float(raw.get("res"), "Grid resolution"),
        sigma=_parse_float(raw.get("sigma"), "Surface sigma"),
        surface_level=_parse_float(
            raw.get("surface_level", DEFAULT_RUN_CONFIG.surface.level),
            "Surface level",
        ),
        surface_padding=_parse_float(
            raw.get("surface_padding", DEFAULT_RUN_CONFIG.surface.padding),
            "Surface padding",
        ),
        max_voxels=_parse_int(
            raw.get("max_voxels"),
            "Maximum voxels",
            default=DEFAULT_RUN_CONFIG.surface.max_voxels,
        ),
        adaptive_resolution=_boolean(raw.get("adaptive_resolution"), DEFAULT_RUN_CONFIG.surface.adaptive_resolution),
        max_adaptive_resolution=_parse_float(
            raw.get("max_adaptive_resolution", DEFAULT_RUN_CONFIG.surface.max_adaptive_resolution),
            "Maximum adaptive resolution",
        ),
        parameterization_method=_choice(
            raw.get("parameterization_method", DEFAULT_RUN_CONFIG.parameterization.method),
            "Parameterization method",
            {"auto", "lscm", "harmonic", "slim", "spherical", "cylindrical"},
        ),
        slim_iterations=_parse_int(
            raw.get("slim_iterations"),
            "SLIM iterations",
            default=DEFAULT_RUN_CONFIG.parameterization.slim_iterations,
        ),
        slim_boundary_constraint_weight=_parse_float(
            raw.get(
                "slim_boundary_constraint_weight",
                DEFAULT_RUN_CONFIG.parameterization.slim_boundary_constraint_weight,
            ),
            "SLIM boundary constraint weight",
        ),
        patch_gap=_parse_float(raw.get("patch_gap", DEFAULT_RUN_CONFIG.optcuts.patch_gap), "Patch gap"),
        optcuts_bin=_text_or_default(raw.get("optcuts_bin"), DEFAULT_RUN_CONFIG.optcuts.optcuts_bin),
        expected_optcuts_sha256=expected_sha,
        optcuts_lambda=_parse_float(
            raw.get("optcuts_lambda", DEFAULT_RUN_CONFIG.optcuts.optcuts_lambda_init),
            "OptCuts lambda",
        ),
        optcuts_distortion_bound=_parse_float(
            raw.get("optcuts_distortion_bound", DEFAULT_RUN_CONFIG.optcuts.optcuts_distortion_bound),
            "OptCuts distortion bound",
        ),
        optcuts_initial_cut_option=_parse_int(
            raw.get("optcuts_initial_cut_option"),
            "OptCuts initial cut option",
            default=DEFAULT_RUN_CONFIG.optcuts.optcuts_initial_cut_option,
        ),
        optcuts_use_bijectivity=_boolean(
            raw.get("optcuts_use_bijectivity"),
            DEFAULT_RUN_CONFIG.optcuts.optcuts_use_bijectivity,
        ),
        optcuts_timeout=_parse_float(
            raw.get("optcuts_timeout", DEFAULT_RUN_CONFIG.optcuts.timeout_sec),
            "OptCuts timeout",
        ),
        residue_fragmentation_weight=residue_fragmentation_weight,
        output_root=output_root,
        resume=(run_mode == "resume"),
        run_mode=run_mode,
        chain_selection_mode=selection_mode,
        manifest_path=manifest_path,
        repetitions=repetitions,
        warmup_runs=warmup_runs,
        formal_mode=formal_mode,
        max_workers=max_workers,
        benchmark_purpose=benchmark_purpose,
        execution_profile=execution_profile,
        optcuts_variants=optcuts_variants,
        include_topology_ablation=_boolean(
            raw.get("include_topology_ablation"),
            DEFAULT_BENCHMARK_CONFIG.include_topology_ablation,
        ),
        threads_per_worker=_parse_int(
            raw.get("threads_per_worker"),
            "Threads per worker",
            default=DEFAULT_BENCHMARK_CONFIG.threads_per_worker,
        ),
        worker_timeout_sec=_parse_float(
            raw.get("worker_timeout_sec", DEFAULT_BENCHMARK_CONFIG.worker_timeout_sec),
            "Worker timeout",
        ),
        worker_memory_limit_mb=_optional_float(
            raw.get("worker_memory_limit_mb"),
            "Worker memory limit",
        ),
        raster_size=_parse_int(
            raw.get("raster_size"),
            "Raster size",
            default=DEFAULT_BENCHMARK_CONFIG.raster_size,
        ),
        min_chain_residues=_parse_int(
            raw.get("min_chain_residues"),
            "Minimum chain residues",
            default=DEFAULT_BENCHMARK_CONFIG.min_chain_residues,
        ),
        per_face_sample_size_per_patch=_parse_int(
            raw.get("per_face_sample_size_per_patch"),
            "Per-face sample size",
            default=DEFAULT_BENCHMARK_CONFIG.per_face_sample_size_per_patch,
        ),
        bootstrap_iterations=_parse_int(
            raw.get("bootstrap_iterations"),
            "Bootstrap iterations",
            default=DEFAULT_BENCHMARK_CONFIG.bootstrap_iterations,
        ),
        random_seed=_parse_int(
            raw.get("random_seed"),
            "Random seed",
            default=DEFAULT_BENCHMARK_CONFIG.random_seed,
        ),
        expected_git_commit=_optional_text(raw.get("expected_git_commit")).lower(),
        coordinate_audit_path=coordinate_audit_path,
        expected_coordinate_audit_sha256=expected_coordinate_audit_sha256,
    )
    form.to_config().validate()
    return form


def _text_or_default(value: object, default: str) -> str:
    text = str(value or "").strip()
    return text or str(default)


def _optional_text(value: object) -> str:
    return str(value or "").strip()


def _optional_output_root(value: object, label: str) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    if os.path.exists(text) and not os.path.isdir(text):
        raise ConfigurationError(f"{label} is not a folder: {text}")
    return text


def _parse_float(value: object, label: str) -> float:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{label} must be a number.") from exc


def _boolean(value: object, default: bool) -> bool:
    if value is None:
        return bool(default)
    if not isinstance(value, bool):
        raise ConfigurationError("Expected a boolean value.")
    return value


def _choice(value: object, label: str, choices: set[str]) -> str:
    normalized = str(value or "").strip().lower()
    if normalized not in choices:
        allowed = ", ".join(sorted(choices))
        raise ConfigurationError(f"{label} must be one of: {allowed}.")
    return normalized


def _method_tuple(value: object, *, default: tuple[str, ...]) -> tuple[str, ...]:
    methods = tuple(default if value is None else value)
    if not methods:
        raise ConfigurationError("At least one benchmark method must be selected.")
    return methods


def _optional_int(value: object, label: str) -> int | None:
    if value is None or str(value).strip() == "":
        return None
    return _parse_int(value, label)


def _optional_float(value: object, label: str) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    return _parse_float(value, label)


def _parse_int(value: object, label: str, default: int | None = None) -> int:
    if (value is None or str(value).strip() == "") and default is not None:
        return int(default)
    try:
        return int(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{label} must be an integer.") from exc
