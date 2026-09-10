"""Central configuration objects for TopoPPI."""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field, replace
from numbers import Integral, Real
from pathlib import Path
from typing import Dict, Optional

from topoppi.benchmark_methods import (
    OPTCUTS_VARIANTS,
    RESIDUE_AWARE_BASELINE,
    RESIDUE_AWARE_OPTCUTS_METHODS,
    resolved_optcuts_variants,
)
from topoppi.errors import ConfigurationError

DEFAULT_RESIDUE_FRAGMENTATION_WEIGHT = 20.0
DEFAULT_CONTACT_DISTANCE_ANGSTROM = 6.0
DEFAULT_INTERFACE_CUTOFF_ANGSTROM = 4.0
DEFAULT_MIN_INTERACTION_RESIDUES = 1


@dataclass(frozen=True)
class SurfaceConfig:
    """Parameters for molecular surface generation."""

    grid_resolution: float = 1.0
    sigma: float = 1.0
    level: float = 0.1
    padding: float = 10.0
    max_voxels: int = 40_000_000
    adaptive_resolution: bool = True
    max_adaptive_resolution: float = 2.0
    smoothing_iterations: int = 0

    def validate(self) -> None:
        _require_positive("surface.grid_resolution", self.grid_resolution)
        _require_positive("surface.sigma", self.sigma)
        _require_positive("surface.level", self.level)
        _require_positive("surface.padding", self.padding)
        _require_positive_integer("surface.max_voxels", self.max_voxels)
        _require_boolean("surface.adaptive_resolution", self.adaptive_resolution)
        _require_positive("surface.max_adaptive_resolution", self.max_adaptive_resolution)
        _require_non_negative_integer("surface.smoothing_iterations", self.smoothing_iterations)
        if self.max_adaptive_resolution < self.grid_resolution:
            raise ConfigurationError("surface.max_adaptive_resolution must be >= surface.grid_resolution.")


@dataclass(frozen=True)
class TopologyConfig:
    """Parameters for interface patch extraction and sanitation."""

    distance_cutoff: float = DEFAULT_INTERFACE_CUTOFF_ANGSTROM
    min_patch_area_angstrom2: float = 10.0
    min_patch_vertices: int = 3
    degenerate_face_area: float = 1e-9
    max_edge_face_incidence: int = 2

    def validate(self) -> None:
        _require_positive("topology.distance_cutoff", self.distance_cutoff)
        _require_non_negative("topology.min_patch_area_angstrom2", self.min_patch_area_angstrom2)
        _require_positive_integer("topology.min_patch_vertices", self.min_patch_vertices)
        _require_positive("topology.degenerate_face_area", self.degenerate_face_area)
        _require_positive_integer("topology.max_edge_face_incidence", self.max_edge_face_incidence)
        if self.max_edge_face_incidence != 2:
            raise ConfigurationError("topology.max_edge_face_incidence must be 2 for a two-manifold surface.")


@dataclass(frozen=True)
class ParameterizationConfig:
    """Parameters for mesh cleanup and UV parameterization."""

    method: str = "auto"
    min_vertices: int = 3
    min_face_area: float = 1e-12
    min_angle_deg: float = 1e-6
    max_aspect_ratio: float = 1e12
    uv_epsilon: float = 1e-6
    expected_euler_characteristic: int = 1
    expected_boundary_loops: int = 1
    lscm_pin_a: tuple[float, float] = (0.0, 0.0)
    lscm_pin_b: tuple[float, float] = (1.0, 0.0)
    slim_iterations: int = 20
    slim_boundary_constraint_weight: float = 1e11

    def validate(self) -> None:
        method = self.method.strip().lower()
        if method not in {"auto", "lscm", "harmonic", "slim", "spherical", "cylindrical"}:
            raise ConfigurationError(f"Unsupported parameterization method: {self.method}")
        _require_positive_integer("parameterization.min_vertices", self.min_vertices)
        _require_positive("parameterization.min_face_area", self.min_face_area)
        _require_positive("parameterization.min_angle_deg", self.min_angle_deg)
        _require_positive("parameterization.max_aspect_ratio", self.max_aspect_ratio)
        _require_positive("parameterization.uv_epsilon", self.uv_epsilon)
        _require_integer("parameterization.expected_euler_characteristic", self.expected_euler_characteristic)
        _require_positive_integer("parameterization.expected_boundary_loops", self.expected_boundary_loops)
        if self.expected_euler_characteristic != 1 or self.expected_boundary_loops != 1:
            raise ConfigurationError(
                "Parameterization domains must be connected disks (Euler characteristic 1, one boundary loop)."
            )
        _require_positive_integer("parameterization.slim_iterations", self.slim_iterations)
        _require_positive(
            "parameterization.slim_boundary_constraint_weight",
            self.slim_boundary_constraint_weight,
        )


@dataclass(frozen=True)
class OptCutsConfig:
    """Parameters for invoking OptCuts."""

    optcuts_bin: str = "OptCuts_bin"
    patch_gap: float = 0.08
    optcuts_mode: int = 10
    optcuts_headless_mode: int = 100
    optcuts_prog_mode: int = 1
    optcuts_method_type: int = 0
    optcuts_initial_cut_option: int = 0
    optcuts_lambda_init: float = 0.999
    optcuts_distortion_bound: float = 4.1
    optcuts_use_bijectivity: bool = True
    optcuts_quick_mode: bool = False
    optcuts_quick_distortion_bound: float = 6.0
    optcuts_quick_lambda_init: float = 0.95
    optcuts_quick_use_bijectivity: bool = False
    optcuts_output_tag: str = "patch"
    save_optcuts_frames: bool = False
    optcuts_frame_stride: int = 1
    optcuts_min_frame_long_edge: int = 3200
    optcuts_frames_dir: str = ""
    optcuts_copy_raw_gif: bool = True
    optcuts_env_var: str = "TOPOPPI_OPTCUTS_BIN"
    expected_binary_sha256: str = ""
    use_input_uv: bool = False
    residue_fragmentation_weight: float = 0.0
    timeout_sec: float = 600.0

    def validate(self) -> None:
        if not self.optcuts_bin.strip():
            raise ConfigurationError("optcuts.optcuts_bin is required.")
        _require_non_negative("optcuts.patch_gap", self.patch_gap)
        _require_non_negative(
            "optcuts.residue_fragmentation_weight",
            self.residue_fragmentation_weight,
        )
        _require_positive_integer("optcuts.optcuts_mode", self.optcuts_mode)
        _require_positive_integer("optcuts.optcuts_headless_mode", self.optcuts_headless_mode)
        _require_non_negative_integer("optcuts.optcuts_prog_mode", self.optcuts_prog_mode)
        _require_non_negative_integer("optcuts.optcuts_method_type", self.optcuts_method_type)
        if self.optcuts_method_type != 0:
            raise ConfigurationError(
                "optcuts.optcuts_method_type must be 0; other upstream methods are not TopoPPI OptCuts treatments."
            )
        _require_non_negative_integer("optcuts.optcuts_initial_cut_option", self.optcuts_initial_cut_option)
        if self.optcuts_initial_cut_option not in {0, 1}:
            raise ConfigurationError("optcuts.optcuts_initial_cut_option must be 0 or 1.")
        _require_positive("optcuts.optcuts_lambda_init", self.optcuts_lambda_init)
        _require_positive("optcuts.optcuts_distortion_bound", self.optcuts_distortion_bound)
        _require_positive("optcuts.optcuts_quick_distortion_bound", self.optcuts_quick_distortion_bound)
        _require_positive("optcuts.optcuts_quick_lambda_init", self.optcuts_quick_lambda_init)
        _require_positive_integer("optcuts.optcuts_frame_stride", self.optcuts_frame_stride)
        _require_non_negative_integer("optcuts.optcuts_min_frame_long_edge", self.optcuts_min_frame_long_edge)
        for name in (
            "optcuts_use_bijectivity",
            "optcuts_quick_mode",
            "optcuts_quick_use_bijectivity",
            "save_optcuts_frames",
            "optcuts_copy_raw_gif",
            "use_input_uv",
        ):
            _require_boolean(f"optcuts.{name}", getattr(self, name))
        if not self.optcuts_output_tag.strip():
            raise ConfigurationError("optcuts.optcuts_output_tag is required.")
        for name, value in (
            ("optcuts.optcuts_lambda_init", self.optcuts_lambda_init),
            ("optcuts.optcuts_quick_lambda_init", self.optcuts_quick_lambda_init),
        ):
            if value >= 1.0:
                raise ConfigurationError(f"{name} must be in (0, 1).")
        for name, value in (
            ("optcuts.optcuts_distortion_bound", self.optcuts_distortion_bound),
            ("optcuts.optcuts_quick_distortion_bound", self.optcuts_quick_distortion_bound),
        ):
            if value <= 4.0:
                raise ConfigurationError(
                    f"{name} must be > 4.0, the lower limit of OptCuts' symmetric Dirichlet energy."
                )
        expected_sha = self.expected_binary_sha256.strip().lower()
        if expected_sha and (len(expected_sha) != 64 or any(char not in "0123456789abcdef" for char in expected_sha)):
            raise ConfigurationError("optcuts.expected_binary_sha256 must be a 64-character SHA-256 digest.")
        _require_positive("optcuts.timeout_sec", self.timeout_sec)

    def for_headless(self) -> "OptCutsConfig":
        return replace(self, optcuts_mode=self.optcuts_headless_mode)


@dataclass(frozen=True)
class VisualizationConfig:
    """Parameters for 2D/3D plotting and residue annotation."""

    show_plot: bool = False
    min_points: int = DEFAULT_MIN_INTERACTION_RESIDUES
    residue_scope: str = "interaction"
    color_by_interaction_type: bool = True
    use_geometric_interaction_fallback: bool = False
    vdw_distance: float = 4.5
    ionic_distance: float = 5.0
    polar_contact_distance: float = 3.8
    mesh_fill_alpha: float = 0.20
    mesh_line_alpha: float = 0.60
    label_offset: float = 0.04
    map_style: str = "markers"
    view: str = "atlas"
    show_mesh: bool = True
    surface_projection: str = "orthographic"
    surface_elevation: float = 73.0
    surface_azimuth: float = -90.0
    surface_zoom: float = 1.0
    highlight_residues: tuple[str, ...] = ()
    annotation_file: str = ""
    annotation_label: str = "Value"
    value_min: Optional[float] = None
    value_max: Optional[float] = None
    footprint_labels: str = "all"
    show_seams: bool = True
    show_residue_borders: bool = True
    footprint_color: str = "#DCE8EF"
    highlight_color: str = "#A64D79"
    missing_color: str = "#D9D9D9"

    def validate(self) -> None:
        _require_positive_integer("visualization.min_points", self.min_points)
        if self.map_style not in {"markers", "footprints"}:
            raise ConfigurationError("visualization.map_style must be 'markers' or 'footprints'.")
        if self.view not in {"atlas", "surface"}:
            raise ConfigurationError("visualization.view must be 'atlas' or 'surface'.")
        if self.surface_projection not in {"orthographic", "perspective"}:
            raise ConfigurationError("visualization.surface_projection must be 'orthographic' or 'perspective'.")
        for name in ("surface_elevation", "surface_azimuth"):
            _require_finite_number(f"visualization.{name}", getattr(self, name))
        _require_positive("visualization.surface_zoom", self.surface_zoom)
        if self.footprint_labels not in {"all", "highlighted", "none"}:
            raise ConfigurationError("visualization.footprint_labels must be 'all', 'highlighted', or 'none'.")
        for name in ("show_seams", "show_residue_borders", "show_mesh"):
            _require_boolean(f"visualization.{name}", getattr(self, name))
        for name in ("value_min", "value_max"):
            if getattr(self, name) is not None:
                _require_finite_number(f"visualization.{name}", getattr(self, name))
        if self.value_min is not None and self.value_max is not None and self.value_min >= self.value_max:
            raise ConfigurationError("visualization.value_min must be smaller than value_max.")
        if isinstance(self.highlight_residues, str):
            raise ConfigurationError("visualization.highlight_residues must be a sequence of residue keys.")
        if self.map_style == "footprints" or self.view == "surface":
            from matplotlib.colors import is_color_like

            for name in ("footprint_color", "highlight_color", "missing_color"):
                if not is_color_like(getattr(self, name)):
                    raise ConfigurationError(f"visualization.{name} must be a valid colour.")
        if str(self.residue_scope).strip().lower() not in {"interaction", "patch"}:
            raise ConfigurationError("visualization.residue_scope must be 'interaction' or 'patch'.")
        _require_boolean(
            "visualization.color_by_interaction_type",
            self.color_by_interaction_type,
        )
        _require_boolean(
            "visualization.use_geometric_interaction_fallback",
            self.use_geometric_interaction_fallback,
        )
        for name in (
            "vdw_distance",
            "ionic_distance",
            "polar_contact_distance",
            "label_offset",
        ):
            _require_positive(f"visualization.{name}", getattr(self, name))
        _require_alpha("visualization.mesh_fill_alpha", self.mesh_fill_alpha)
        _require_alpha("visualization.mesh_line_alpha", self.mesh_line_alpha)


@dataclass(frozen=True)
class GUIConfig:
    """Defaults for the Tkinter GUI."""

    window_width: int = 1400
    window_height: int = 900
    min_window_width: int = 960
    min_window_height: int = 600
    sidebar_width: int = 460
    ttk_theme: str = "clam"
    tk_scaling: float = 1.2
    font_family: str = "Segoe UI"
    font_fallbacks: tuple[str, ...] = ("Segoe UI", "Arial", "DejaVu Sans", "Noto Sans", "Liberation Sans")
    font_size: int = 10
    header_font_size: int = 16
    log_visible_lines: int = 7
    ui_poll_interval_ms: int = 80
    default_patch_cutoff: float = DEFAULT_INTERFACE_CUTOFF_ANGSTROM
    default_min_points: int = DEFAULT_MIN_INTERACTION_RESIDUES
    default_residue_color: str = "#d62728"
    auto_save_single_run: bool = True
    label_font_size: int = 9
    label_font_min_size: int = 5
    label_font_max_size: int = 20
    figure_dpi: int = 300
    benchmark_output_folder: str = "benchmark_results_resume"
    default_benchmark_run_mode: str = "resume"


@dataclass(frozen=True)
class TopoPPIRunConfig:
    """Configuration for one protein-protein interface map run."""

    pdb_file: str = ""
    chain_a: str = "A"
    chain_b: str = "B"
    output_file: str = "interface_map.png"
    prolif_file: Optional[str] = None
    interaction_source: str = "prolif"
    atlas_output: Optional[str] = None
    contact_distance_angstrom: float = DEFAULT_CONTACT_DISTANCE_ANGSTROM
    surface: SurfaceConfig = field(default_factory=SurfaceConfig)
    topology: TopologyConfig = field(default_factory=TopologyConfig)
    parameterization: ParameterizationConfig = field(default_factory=ParameterizationConfig)
    optcuts: OptCutsConfig = field(
        default_factory=lambda: OptCutsConfig(residue_fragmentation_weight=DEFAULT_RESIDUE_FRAGMENTATION_WEIGHT)
    )
    visualization: VisualizationConfig = field(default_factory=VisualizationConfig)

    def validate(self) -> None:
        input_path = Path(self.pdb_file)
        if not self.pdb_file:
            raise ConfigurationError("Choose an input PDB or mmCIF structure.")
        if not input_path.is_file():
            raise ConfigurationError(f"Input structure was not found: {self.pdb_file}")
        if not str(self.chain_a).strip():
            raise ConfigurationError("Choose the surface chain (Chain A).")
        if not str(self.chain_b).strip():
            raise ConfigurationError("Choose the partner chain (Chain B).")
        if str(self.chain_a).strip() == str(self.chain_b).strip():
            raise ConfigurationError("Choose two different chains for Chain A and Chain B.")
        if self.interaction_source not in {"prolif", "geometric"}:
            raise ConfigurationError("interaction_source must be 'prolif' or 'geometric'.")
        if self.interaction_source == "geometric" and self.prolif_file:
            raise ConfigurationError("Choose geometric interactions or a ProLIF file; remove the other setting.")
        if self.prolif_file and not Path(self.prolif_file).is_file():
            raise ConfigurationError(f"ProLIF JSON was not found: {self.prolif_file}")
        if not str(self.output_file).strip():
            raise ConfigurationError("Choose an output image path.")
        if Path(self.output_file).expanduser().resolve() == input_path.expanduser().resolve():
            raise ConfigurationError("The output image path matches the input structure path.")
        if self.atlas_output:
            atlas_path = Path(self.atlas_output).expanduser().resolve()
            if atlas_path in {input_path.expanduser().resolve(), Path(self.output_file).expanduser().resolve()}:
                raise ConfigurationError("Choose an atlas output path separate from the input structure and output image.")
        _require_positive("contact_distance_angstrom", self.contact_distance_angstrom)
        self.surface.validate()
        self.topology.validate()
        self.parameterization.validate()
        self.optcuts.validate()
        self.visualization.validate()

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class BenchmarkConfig:
    """Configuration for benchmark runs."""

    input_folder: str
    output_root: str
    chain_a: str = "A"
    chain_b: str = "B"
    chain_selection_mode: str = "configured"
    manifest_path: str = ""
    surface: SurfaceConfig = field(default_factory=SurfaceConfig)
    topology: TopologyConfig = field(default_factory=TopologyConfig)
    parameterization: ParameterizationConfig = field(default_factory=ParameterizationConfig)
    optcuts: OptCutsConfig = field(default_factory=lambda: OptCutsConfig().for_headless())
    raster_size: int = 256
    max_workers: Optional[int] = 1
    repetitions: int = 3
    warmup_runs: int = 0
    formal_mode: bool = False
    expected_git_commit: str = ""
    coordinate_audit_path: str = ""
    expected_coordinate_audit_sha256: str = ""
    benchmark_purpose: str = "performance"
    execution_profile: str = "comparative"
    optcuts_variants: tuple[str, ...] | None = None
    include_topology_ablation: bool = True
    random_seed: int = 20260817
    bootstrap_iterations: int = 2000
    threads_per_worker: int = 1
    contact_distance_angstrom: float = DEFAULT_CONTACT_DISTANCE_ANGSTROM
    per_face_sample_size_per_patch: int = 128
    worker_timeout_sec: float = 7200.0
    worker_memory_limit_mb: float | None = None
    worker_poll_interval_sec: float = 0.05
    show_tqdm: bool = True
    resume: bool = True
    checkpoint_interval_structures: int = 1
    min_chain_residues: int = 11
    checkpoint_filename: str = "benchmark_checkpoint.json"
    report_filename: str = "benchmark_report.json"
    summary_filename: str = "benchmark_summary.csv"
    manifest_filename: str = "benchmark_manifest.csv"
    failures_filename: str = "benchmark_failures.csv"
    per_patch_filename: str = "benchmark_per_patch.csv"
    per_face_sample_filename: str = "benchmark_per_face_sample.csv"
    per_residue_filename: str = "benchmark_per_residue.csv.gz"
    provenance_filename: str = "benchmark_provenance.csv.gz"
    optcuts_execution_filename: str = "benchmark_optcuts_executions.jsonl.gz"
    artifact_checksums_filename: str = "benchmark_artifact_checksums.json"
    worker_log_folder: str = "worker_logs"

    def validate(self) -> None:
        input_path = Path(self.input_folder)
        if not input_path.is_dir():
            raise ConfigurationError(f"Benchmark input folder does not exist: {self.input_folder}")
        if not self.output_root.strip():
            raise ConfigurationError("output_root is required.")
        selection_mode = self.chain_selection_mode.strip().lower()
        self._validate_selection(selection_mode)
        self.surface.validate()
        self.topology.validate()
        self.parameterization.validate()
        self.optcuts.validate()
        self._validate_runtime_limits()
        self._validate_artifact_names()
        if self.formal_mode:
            self._validate_formal_mode(selection_mode)

    def _validate_selection(self, selection_mode: str) -> None:
        if selection_mode not in {"configured", "auto_contact", "manifest"}:
            raise ConfigurationError("benchmark.chain_selection_mode must be configured, auto_contact, or manifest.")
        if selection_mode == "configured":
            if not self.chain_a.strip() or not self.chain_b.strip():
                raise ConfigurationError("benchmark chain IDs are required in configured mode.")
            if self.chain_a.strip() == self.chain_b.strip():
                raise ConfigurationError("benchmark chain IDs must be different in configured mode.")
        if selection_mode == "manifest":
            manifest = Path(self.manifest_path)
            if not self.manifest_path or not manifest.is_file():
                raise ConfigurationError("benchmark.manifest_path must be an existing file in manifest mode.")

    def _validate_runtime_limits(self) -> None:
        purpose = self.benchmark_purpose.strip().lower()
        if purpose not in {"quality", "performance"}:
            raise ConfigurationError("benchmark.benchmark_purpose must be quality or performance.")
        execution_profile = self.execution_profile.strip().lower()
        if execution_profile not in {"comparative", "operational_optcuts"}:
            raise ConfigurationError("benchmark.execution_profile must be comparative or operational_optcuts.")
        variants = self.resolved_optcuts_variants()
        if not variants:
            raise ConfigurationError("benchmark.optcuts_variants must enable at least one OptCuts variant.")
        if len(variants) != len(set(variants)):
            raise ConfigurationError("benchmark.optcuts_variants must not contain duplicates.")
        unsupported = sorted(set(variants) - set(OPTCUTS_VARIANTS))
        if unsupported:
            raise ConfigurationError("Unsupported benchmark.optcuts_variants: " + ", ".join(unsupported))
        for residue_aware_method in set(variants) & set(RESIDUE_AWARE_OPTCUTS_METHODS):
            if self.optcuts.residue_fragmentation_weight <= 0.0:
                raise ConfigurationError(f"{residue_aware_method} requires optcuts.residue_fragmentation_weight > 0.")
            baseline = RESIDUE_AWARE_BASELINE[residue_aware_method]
            if execution_profile != "operational_optcuts" and baseline not in variants:
                raise ConfigurationError(f"{residue_aware_method} requires its matched baseline {baseline}.")
        if "optcuts_lscm_initialized" in variants and "optcuts_automatic" not in variants:
            raise ConfigurationError(
                "optcuts_lscm_initialized requires optcuts_automatic for the matched initialization comparison."
            )
        if self.optcuts.residue_fragmentation_weight > 0.0 and not set(variants) & set(RESIDUE_AWARE_OPTCUTS_METHODS):
            raise ConfigurationError(
                "A positive residue_fragmentation_weight requires the residue_aware_optcuts variant."
            )
        if self.include_topology_ablation and "optcuts_automatic" not in variants:
            raise ConfigurationError("include_topology_ablation requires optcuts_automatic.")
        if execution_profile == "operational_optcuts":
            if purpose != "performance":
                raise ConfigurationError("operational_optcuts is a performance-only execution profile.")
            if len(variants) != 1 or variants[0] not in {
                "optcuts_automatic",
                "residue_aware_optcuts",
            }:
                raise ConfigurationError("operational_optcuts requires exactly one automatic OptCuts variant.")
            if self.include_topology_ablation:
                raise ConfigurationError("operational_optcuts excludes comparison-only topology ablation.")
            if variants == ("optcuts_automatic",) and self.optcuts.residue_fragmentation_weight != 0.0:
                raise ConfigurationError("operational optcuts_automatic requires residue_fragmentation_weight=0.")
        _require_positive_integer("benchmark.raster_size", self.raster_size)
        _require_positive_integer("benchmark.min_chain_residues", self.min_chain_residues)
        _require_positive_integer("benchmark.repetitions", self.repetitions)
        _require_non_negative_integer("benchmark.warmup_runs", self.warmup_runs)
        _require_positive_integer("benchmark.bootstrap_iterations", self.bootstrap_iterations)
        _require_positive_integer("benchmark.threads_per_worker", self.threads_per_worker)
        _require_positive("benchmark.contact_distance_angstrom", self.contact_distance_angstrom)
        _require_positive_integer(
            "benchmark.per_face_sample_size_per_patch",
            self.per_face_sample_size_per_patch,
        )
        _require_positive("benchmark.worker_timeout_sec", self.worker_timeout_sec)
        if self.worker_memory_limit_mb is not None:
            _require_positive("benchmark.worker_memory_limit_mb", self.worker_memory_limit_mb)
        _require_positive("benchmark.worker_poll_interval_sec", self.worker_poll_interval_sec)
        _require_positive_integer(
            "benchmark.checkpoint_interval_structures",
            self.checkpoint_interval_structures,
        )
        if self.max_workers is not None:
            _require_positive_integer("benchmark.max_workers", self.max_workers)
        _require_non_negative_integer("benchmark.random_seed", self.random_seed)
        expected_commit = str(self.expected_git_commit).strip().lower()
        if expected_commit and (
            len(expected_commit) not in {40, 64}
            or any(character not in "0123456789abcdef" for character in expected_commit)
        ):
            raise ConfigurationError("benchmark.expected_git_commit must be a 40- or 64-character Git object ID.")
        audit_path = str(self.coordinate_audit_path).strip()
        audit_sha256 = str(self.expected_coordinate_audit_sha256).strip().lower()
        if bool(audit_path) != bool(audit_sha256):
            raise ConfigurationError(
                "benchmark.coordinate_audit_path and "
                "benchmark.expected_coordinate_audit_sha256 must be provided together."
            )
        if audit_sha256 and (
            len(audit_sha256) != 64 or any(character not in "0123456789abcdef" for character in audit_sha256)
        ):
            raise ConfigurationError("benchmark.expected_coordinate_audit_sha256 must be a 64-character SHA-256.")
        for name in ("formal_mode", "include_topology_ablation", "show_tqdm", "resume"):
            _require_boolean(f"benchmark.{name}", getattr(self, name))

    def _validate_artifact_names(self) -> None:
        artifact_names = (
            self.checkpoint_filename,
            self.report_filename,
            self.summary_filename,
            self.manifest_filename,
            self.failures_filename,
            self.per_patch_filename,
            self.per_face_sample_filename,
            self.per_residue_filename,
            self.provenance_filename,
            self.optcuts_execution_filename,
            self.artifact_checksums_filename,
        )
        for name in artifact_names:
            _require_safe_output_name("benchmark artifact filename", name)
        if len(set(artifact_names)) != len(artifact_names):
            raise ConfigurationError("Benchmark artifact filenames must be distinct.")
        _require_safe_output_name("benchmark.worker_log_folder", self.worker_log_folder)
        if self.worker_log_folder in artifact_names:
            raise ConfigurationError("benchmark.worker_log_folder must not collide with an artifact filename.")

    def _validate_formal_mode(self, selection_mode: str) -> None:
        if selection_mode != "manifest":
            raise ConfigurationError("Formal benchmark mode requires an explicit dataset manifest.")
        purpose = self.benchmark_purpose.strip().lower()
        if purpose == "performance":
            if self.repetitions < 3:
                raise ConfigurationError("Formal performance benchmarks require at least three measured repetitions.")
            if self.warmup_runs < 1:
                raise ConfigurationError("Formal performance benchmarks require at least one warm-up repetition.")
            if self.max_workers != 1:
                raise ConfigurationError("Formal performance benchmarks require max_workers=1 for uncontended timing.")
        else:
            if self.repetitions != 1 or self.warmup_runs != 0:
                raise ConfigurationError(
                    "Formal quality benchmarks require repetitions=1 and warmup_runs=0; "
                    "use a separate performance subset for repeated timing."
                )
        scheduled_optcuts_arms = len(self.resolved_optcuts_variants()) + int(self.include_topology_ablation)
        scheduled_solver_budget_sec = scheduled_optcuts_arms * float(self.optcuts.timeout_sec)
        if float(self.worker_timeout_sec) <= scheduled_solver_budget_sec:
            raise ConfigurationError(
                "Formal worker_timeout_sec must exceed the sum of all scheduled "
                "OptCuts method-arm budgets; otherwise method order can censor a later arm."
            )
        if not self.optcuts.expected_binary_sha256.strip():
            raise ConfigurationError("Formal benchmark mode requires optcuts.expected_binary_sha256.")

    def resolved_optcuts_variants(self) -> tuple[str, ...]:
        return resolved_optcuts_variants(
            self.optcuts_variants,
            residue_fragmentation_weight=self.optcuts.residue_fragmentation_weight,
        )

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


DEFAULT_RUN_CONFIG = TopoPPIRunConfig()
DEFAULT_BENCHMARK_CONFIG = BenchmarkConfig(input_folder=".", output_root="benchmark_results")
DEFAULT_GUI_CONFIG = GUIConfig()


def benchmark_config_from_dict(payload: Dict[str, object]) -> BenchmarkConfig:
    """Reconstruct nested benchmark configuration in isolated worker processes."""

    data = dict(payload)
    data["surface"] = SurfaceConfig(**dict(data.get("surface", {})))
    data["topology"] = TopologyConfig(**dict(data.get("topology", {})))
    data["parameterization"] = ParameterizationConfig(**dict(data.get("parameterization", {})))
    data["optcuts"] = OptCutsConfig(**dict(data.get("optcuts", {})))
    if data.get("optcuts_variants") is not None:
        raw_variants = data["optcuts_variants"]
        if isinstance(raw_variants, (str, bytes)) or not isinstance(raw_variants, (list, tuple)):
            raise ConfigurationError("benchmark.optcuts_variants must be a list of method names.")
        data["optcuts_variants"] = tuple(str(value) for value in raw_variants)
    return BenchmarkConfig(**data)


def _require_positive(name: str, value: float) -> None:
    numeric = _require_finite_number(name, value)
    if numeric <= 0:
        raise ConfigurationError(f"{name} must be > 0.")


def _require_non_negative(name: str, value: float) -> None:
    numeric = _require_finite_number(name, value)
    if numeric < 0:
        raise ConfigurationError(f"{name} must be >= 0.")


def _require_finite_number(name: str, value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    return float(value)


def _require_integer(name: str, value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ConfigurationError(f"{name} must be an integer.")


def _require_positive_integer(name: str, value: int) -> None:
    _require_integer(name, value)
    if int(value) <= 0:
        raise ConfigurationError(f"{name} must be > 0.")


def _require_non_negative_integer(name: str, value: int) -> None:
    _require_integer(name, value)
    if int(value) < 0:
        raise ConfigurationError(f"{name} must be >= 0.")


def _require_boolean(name: str, value: bool) -> None:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")


def _require_safe_output_name(name: str, value: str) -> None:
    text = str(value or "").strip()
    if not text or text in {".", ".."} or "/" in text or "\\" in text:
        raise ConfigurationError(f"{name} must be one non-empty filename without path separators.")


def _require_alpha(name: str, value: float) -> None:
    numeric = _require_finite_number(name, value)
    if not 0.0 <= numeric <= 1.0:
        raise ConfigurationError(f"{name} must be in [0, 1].")
