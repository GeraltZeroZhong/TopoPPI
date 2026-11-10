"""Reusable TopoPPI pipeline API."""

from __future__ import annotations

import logging
import platform
import sys
import time
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Mapping, Optional

import matplotlib.pyplot as plt
import trimesh

from topoppi import __version__
from topoppi.atlas.footprints import (
    contact_partner_degrees,
    geometric_contact_partner_map,
    residue_aware_residue_weights,
    residue_fragmentation_report,
    source_atom_residue_labels,
)
from topoppi.atlas.uv import set_uv_layout
from topoppi.config import TopoPPIRunConfig
from topoppi.errors import InputDataError, PipelineError
from topoppi.file_utils import sha256_file
from topoppi.interactions.interaction_engine import (
    generate_interactions,
    generate_prolif_interactions,
    load_interaction_document,
    load_interaction_partner_map,
)
from topoppi.interface import prepare_interface
from topoppi.io.io_loader import PDBLoader
from topoppi.json_utils import dump_json_atomic
from topoppi.mesh.parameterization import Parameterizer
from topoppi.mesh.surface import SurfaceGenerator
from topoppi.mesh.topology import TopologyManager
from topoppi.optimization.optcuts import OptCutsUVOptimizer
from topoppi.visualization.visualizer import InterfaceVisualizer, select_patches_for_display

logger = logging.getLogger("topoppi.pipeline")

_SUPPORTED_OUTPUT_SUFFIXES = frozenset({".png", ".tif", ".tiff", ".pdf", ".svg"})


def layout_interface(prepared, *, side=0, optcuts=None, parameterization=None, cancel_event=None, log=None):
    """Lay out one prepared side, leaving the shared native geometry untouched."""
    from topoppi.config import DEFAULT_RUN_CONFIG, ParameterizationConfig

    log = log or logger
    config = TopoPPIRunConfig(pdb_file=prepared.source["path"],
        chain_a=",".join(prepared.partners[side].chains), chain_b=",".join(prepared.partners[1-side].chains),
        optcuts=optcuts or DEFAULT_RUN_CONFIG.optcuts, parameterization=parameterization or ParameterizationConfig())
    optimizer = OptCutsUVOptimizer(config.optcuts, cancel_event=cancel_event)
    artifact = optimizer.preflight_binary()
    patches = _parameterize_patches(prepared.partners[side].layout_patches(), config, log)
    report = _optimize_patches(patches, optimizer, prepared.chain_data(side), config, log,
        interaction_partner_map=prepared.contacts(side), interaction_source="geometric")
    report["optcuts_resolved"] = artifact
    return patches, report


@dataclass
class TopoPPIRunResult:
    """Result metadata for one TopoPPI pipeline run."""

    output_file: str
    prolif_file: Optional[str]
    patch_count: int
    valid_patch_count: int
    elapsed_sec: float
    optimizer_report: Dict[str, object]
    manifest_file: str
    input_sha256: str
    surface_generation: Dict[str, object]
    topology_extraction: Dict[str, object]
    visualization: Dict[str, object]

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


def run_interface_mapping(config: TopoPPIRunConfig, log: Optional[logging.Logger] = None) -> TopoPPIRunResult:
    """Run the TopoPPI single-structure interface mapping pipeline."""

    log = log or logger
    config.validate()
    start_time = time.perf_counter()
    _prepare_output_path(config.output_file)
    chain_data = _load_chain_data(config, log)
    if config.visualization.map_style == "footprints":
        from topoppi.visualization.footprint_rendering import read_residue_annotations, resolve_residue_keys

        try:
            source_keys = source_atom_residue_labels(chain_data[1])
            resolve_residue_keys(config.visualization.highlight_residues, source_keys)
            if config.visualization.annotation_file:
                read_residue_annotations(config.visualization.annotation_file, source_keys)
        except (OSError, ValueError) as exc:
            raise PipelineError(f"Invalid footprint annotations: {exc}") from exc

    input_sha256 = sha256_file(config.pdb_file)
    prolif_file = _resolve_prolif_file(config, log, input_sha256=input_sha256, chain_data=chain_data)
    coords_a, _atoms_a, coords_b, _atoms_b = chain_data
    interaction_partners, interaction_source = _build_interaction_partner_map(
        chain_data,
        prolif_file,
        config,
        log,
        input_sha256=input_sha256,
    )
    prepared = prepare_interface(config.pdb_file, config.chain_a, config.chain_b,
        surface=config.surface, topology=config.topology, sides=(0,), chain_data=chain_data, log=log.info)
    partner = prepared.partners[0]
    surface_report, topology_report = partner.surface_report, partner.topology_report
    patches = partner.layout_patches()
    if not patches:
        raise PipelineError("No interface patches retained; check the chain pair and --cutoff.")
    optimizer = OptCutsUVOptimizer(config.optcuts)
    try:
        optcuts_artifact = optimizer.preflight_binary()
    except RuntimeError as exc:
        raise PipelineError(str(exc)) from exc

    valid_patches = _parameterize_patches(patches, config, log)
    optimizer_report = _optimize_patches(
        valid_patches,
        optimizer,
        chain_data,
        config,
        log,
        interaction_partner_map=interaction_partners,
        interaction_source=interaction_source,
    )
    optimizer_report["optcuts_resolved"] = optcuts_artifact
    visualization_report = _render_output(
        valid_patches,
        chain_data,
        prolif_file,
        config,
        log,
        interaction_partner_map=interaction_partners,
        interaction_source=interaction_source,
        run_metadata={"optimizer_report": optimizer_report, "surface_generation": surface_report,
                      "topology_extraction": topology_report, "interface": prepared.metadata()},
    )

    elapsed = time.perf_counter() - start_time
    log.info("Pipeline finished in %.2fs. Saved to %s", elapsed, config.output_file)
    manifest_path = str(Path(config.output_file).with_suffix(".topoppi.json"))
    result = TopoPPIRunResult(
        output_file=str(config.output_file),
        prolif_file=prolif_file,
        patch_count=len(patches),
        valid_patch_count=len(valid_patches),
        elapsed_sec=float(elapsed),
        optimizer_report=optimizer_report,
        manifest_file=manifest_path,
        input_sha256=input_sha256,
        surface_generation=surface_report,
        topology_extraction=topology_report,
        visualization=visualization_report,
    )
    _write_run_manifest(result, config)
    log.info("Run manifest saved to %s", manifest_path)
    return result


def _prepare_output_path(output_file: str) -> None:
    output_path = Path(output_file)
    suffix = output_path.suffix.lower()
    if suffix not in _SUPPORTED_OUTPUT_SUFFIXES:
        found = suffix or "no extension"
        raise PipelineError(f"Unsupported output image extension ({found}). Choose a .png, .tif, .tiff, .svg, or .pdf file.")
    if output_path.exists() and output_path.is_dir():
        raise PipelineError(f"Output image path is a directory: {output_path}")
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise PipelineError(
            f"Could not create output directory {output_path.parent}: {exc}. "
            "Choose another output path or create the directory first."
        ) from exc


def _resolve_prolif_file(
    config: TopoPPIRunConfig,
    log: logging.Logger,
    *,
    input_sha256: str,
    chain_data=None,
) -> Optional[str]:
    if config.interaction_source == "geometric":
        return None
    if config.prolif_file:
        return config.prolif_file
    try:
        generate = generate_prolif_interactions if config.interaction_source == "prolif" else generate_interactions
        options = {} if config.interaction_source == "prolif" else {"backend": config.interaction_source}
        if chain_data is not None:
            options["chain_data"] = chain_data
        generated = generate(
            config.pdb_file,
            config.chain_a,
            config.chain_b,
            log,
            source_sha256=input_sha256,
            output_dir=str(Path(config.output_file).parent),
            **options,
        )
        return generated
    except Exception as exc:
        backend_label = "ProLIF" if config.interaction_source == "prolif" else config.interaction_source
        if config.visualization.use_geometric_interaction_fallback:
            log.warning(
                "%s interaction generation failed (%s); using the selected geometric fallback.",
                backend_label,
                exc,
            )
            return None
        raise PipelineError(f"{backend_label} interaction generation failed: {exc}") from exc


def _load_chain_data(config: TopoPPIRunConfig, log: logging.Logger):
    log.info("Loading %s...", config.pdb_file)
    try:
        loader = PDBLoader(config.pdb_file)
        available_chains = loader.get_protein_chain_ids()
    except Exception as exc:
        raise InputDataError(f"Could not read the input structure: {exc}") from exc

    missing_chains = [chain_id for chain_id in (config.chain_a, config.chain_b) if chain_id not in available_chains]
    if missing_chains:
        available = ", ".join(available_chains) or "none"
        missing = ", ".join(missing_chains)
        raise InputDataError(f"Selected chain(s) {missing} were not found. Available protein chains: {available}.")

    try:
        coords_a, atoms_a = loader.get_chain_data(config.chain_a)
        coords_b, atoms_b = loader.get_chain_data(config.chain_b)
    except Exception as exc:
        raise InputDataError(f"Could not load the selected chains: {exc}") from exc

    if len(coords_a) == 0:
        raise InputDataError(f"Chain {config.chain_a} has no standard protein atoms.")
    if len(coords_b) == 0:
        raise InputDataError(f"Chain {config.chain_b} has no standard protein atoms.")
    log.debug("Loaded Chain %s: %d atoms", config.chain_a, len(coords_a))
    log.debug("Loaded Chain %s: %d atoms", config.chain_b, len(coords_b))
    return coords_a, atoms_a, coords_b, atoms_b


def _build_interaction_partner_map(
    chain_data,
    prolif_file: Optional[str],
    config: TopoPPIRunConfig,
    log: logging.Logger,
    *,
    input_sha256: str,
) -> tuple[Dict[str, Dict[str, int]], str]:
    coords_a, atoms_a, coords_b, atoms_b = chain_data
    if prolif_file and config.interaction_source != "geometric":
        try:
            partners = load_interaction_partner_map(
                prolif_file,
                atoms_a,
                atoms_b,
                expected_chain_a=config.chain_a,
                expected_chain_b=config.chain_b,
                expected_source_sha256=input_sha256,
            )
        except (OSError, ValueError) as exc:
            raise PipelineError(f"Failed to resolve interaction residues: {exc}") from exc
        if not partners:
            raise PipelineError(
                "The interaction document did not yield any Chain-A/Chain-B residue pairs "
                "that resolve against the selected structure."
            )
        payload, _, _ = load_interaction_document(prolif_file)
        source = str(payload.get("engine") or config.interaction_source)
        log.info(
            "%s interaction definition: %d Chain-A residues, %d residue pairs.",
            source,
            len(partners),
            sum(len(values) for values in partners.values()),
        )
        return partners, source

    cutoff = float(config.contact_distance_angstrom)
    partners = geometric_contact_partner_map(
        coords_a,
        atoms_a,
        coords_b,
        atoms_b,
        distance_cutoff=cutoff,
    )
    log.info(
        "Geometric contacts (<= %.3g Angstrom): %d Chain-A residues, %d residue pairs.",
        cutoff,
        len(partners),
        sum(len(values) for values in partners.values()),
    )
    return partners, "geometric" if config.interaction_source == "geometric" else "geometric_fallback"


def _generate_surface(
    coords_a,
    config: TopoPPIRunConfig,
    log: logging.Logger,
) -> tuple[trimesh.Trimesh, Dict[str, object]]:
    log.info("Generating Gaussian-density implicit molecular surface for Chain %s...", config.chain_a)
    generator = SurfaceGenerator(coords_a, config=config.surface)
    mesh_a = generator.generate_mesh()
    if mesh_a is None or len(mesh_a.vertices) == 0:
        raise PipelineError(f"Failed to generate surface mesh: {generator.last_report}")
    return mesh_a, dict(generator.last_report)


def _extract_patches(
    mesh_a: trimesh.Trimesh,
    coords_b,
    config: TopoPPIRunConfig,
    log: logging.Logger,
) -> tuple[List[trimesh.Trimesh], Dict[str, object]]:
    log.info("Extracting interface patches...")
    manager = TopologyManager(mesh_a, coords_b, config=config.topology)
    patches = manager.get_interface_patches()
    if not patches:
        report = manager.last_report
        log.debug("Topology diagnostics: %s", report)
        if report.get("status") == "no_interface_faces":
            nearest = float(report["nearest_partner_distance_angstrom"])
            raise PipelineError(
                f"No surface faces were within {config.topology.distance_cutoff:g} Å of Chain {config.chain_b}. "
                f"The nearest face was {nearest:.2f} Å away. Check the chain pair or increase --cutoff."
            )
        raise PipelineError(
            "Interface faces were found, but topology preparation retained no usable patch. "
            "Check the chain pair and the minimum patch area or vertex settings."
        )
    return list(patches), dict(manager.last_report)


def _parameterize_patches(
    patches: List[trimesh.Trimesh], config: TopoPPIRunConfig, log: logging.Logger
) -> List[trimesh.Trimesh]:
    log.info("Preparing %d patches for OptCuts...", len(patches))
    valid_patches = []
    parameterizer = Parameterizer(config=config.parameterization)
    for idx, patch in enumerate(patches):
        prepared = parameterizer.prepare_patch(patch)
        if prepared is None:
            log.warning("Skipping patch %d because topology preparation failed.", idx + 1)
            continue
        if config.optcuts.use_input_uv:
            log.info("Flattening patch %d (%d vertices)...", idx + 1, len(patch.vertices))
            uv = parameterizer.flatten_patch(patch)
            if uv is None:
                log.warning("Skipping patch %d due to parameterization failure.", idx + 1)
                continue
            set_uv_layout(patch, uv, key="uv")
        valid_patches.append(patch)

    if not valid_patches:
        raise PipelineError("All patches failed topology preparation or requested UV initialization.")
    return valid_patches


def _optimize_patches(
    patches: List[trimesh.Trimesh],
    optimizer: OptCutsUVOptimizer,
    chain_data,
    config: TopoPPIRunConfig,
    log: logging.Logger,
    *,
    interaction_partner_map: Mapping[str, Mapping[str, int]],
    interaction_source: str,
) -> Dict[str, object]:
    log.info("Running OptCuts UV optimization...")
    _coords_a, atoms_a, _coords_b, _atoms_b = chain_data
    interaction_weights = contact_partner_degrees(interaction_partner_map)
    source_labels = source_atom_residue_labels(atoms_a)
    objective_weights = residue_aware_residue_weights(source_labels, interaction_weights)
    try:
        optimizer.optimize_patches(
            patches,
            source_residue_labels=source_labels,
            residue_weights=objective_weights,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        raise PipelineError(f"OptCuts optimization failed: {exc}") from exc
    report = optimizer.get_last_report()
    weight_definition = (
        f"distinct Chain-B residues paired with each Chain-A residue in {interaction_source} records"
        if interaction_source not in {"geometric", "geometric_fallback"}
        else (
            "distinct Chain-B residues with any heavy-atom pair at distance <= "
            f"{float(config.contact_distance_angstrom):g} Angstrom "
            "(geometric contacts)"
        )
    )
    report["residue_footprint_fragmentation"] = {
        "interaction_residue_source": interaction_source,
        "interaction_weight_definition": weight_definition,
        **residue_fragmentation_report(
            patches,
            source_labels,
            uv_key="uv_optcuts",
            interaction_weights=interaction_weights,
            objective_weights=objective_weights,
        ),
    }
    quality = report["parameterization_quality"]
    topology = report["topology_complexity"]
    fragmentation = report["residue_footprint_fragmentation"]
    log.debug(
        "TopoPPI result: distortion_mean=%.6g, flip_rate=%.6g, seams=%d, "
        "objective_weighted_fragmentation=%.6g, optcuts_time=%.3fs",
        float(quality["distortion"]["mean"]),
        float(quality["flip_rate_mean"]),
        int(topology["seam_edge_count"]),
        float(fragmentation["objective_weighted_fragmentation"]),
        float(report["stability_efficiency"]["total_time_sec"]),
    )
    return report


def _render_output(
    patches: List[trimesh.Trimesh],
    chain_data,
    prolif_file: Optional[str],
    config: TopoPPIRunConfig,
    log: logging.Logger,
    *,
    interaction_partner_map: Mapping[str, Mapping[str, int]],
    interaction_source: str,
    run_metadata: Optional[Dict[str, object]] = None,
) -> Dict[str, object]:
    log.info("Visualizing results...")
    coords_a, atoms_a, coords_b, atoms_b = chain_data
    output_path = Path(config.output_file)
    try:
        viz = InterfaceVisualizer(
            chain_A_atoms=atoms_a,
            chain_A_coords=coords_a,
            chain_B_coords=coords_b,
            chain_B_atoms=atoms_b,
            chain_a_id=config.chain_a,
            chain_b_id=config.chain_b,
            structure_label=Path(config.pdb_file).stem,
            prolif_file=prolif_file,
            config=replace(config.visualization, use_geometric_interaction_fallback=True)
            if config.interaction_source == "geometric" else config.visualization,
            interaction_partner_map=interaction_partner_map,
            contact_distance_angstrom=config.contact_distance_angstrom,
        )
        viz.interaction_residue_source = interaction_source
        displayed_patches, interaction_counts = select_patches_for_display(
            patches,
            viz,
            map_style=config.visualization.map_style,
            min_points=config.visualization.min_points,
        )
        hidden_patch_count = len(patches) - len(displayed_patches)
        display_log = log.info if hidden_patch_count else log.debug
        display_log(
            "Display filter (minimum interaction residues = %d): displayed=%d, hidden=%d; all patches were optimized.",
            config.visualization.min_points,
            len(displayed_patches),
            hidden_patch_count,
        )
        if not displayed_patches:
            raise PipelineError(
                "No interface patch meets the display threshold "
                f"({config.visualization.min_points} interaction residues)."
            )
        figure = viz.plot_patches(
            displayed_patches,
            output_file=str(output_path),
            show=config.visualization.show_plot,
        )
        if config.atlas_output:
            from topoppi.visualization.atlas_io import save_atlas

            save_atlas(config.atlas_output, patches, viz,
                       run_metadata={"config": config.to_dict(), **(run_metadata or {})})
        if not config.visualization.show_plot:
            plt.close(figure)
        report = dict(viz.last_report)
        report["interaction_residue_source"] = interaction_source
        if config.atlas_output:
            report["atlas_file"] = str(Path(config.atlas_output).resolve())
        report["display_filter"] = {
            "policy": "complete_footprints" if config.visualization.map_style == "footprints" else "interaction_threshold",
            "min_points": int(config.visualization.min_points),
            "optimized_patch_count": int(len(patches)),
            "displayed_patch_count": int(len(displayed_patches)),
            "hidden_patch_count": int(hidden_patch_count),
            "patch_interaction_residue_counts": interaction_counts,
        }
        return report
    except (OSError, ValueError) as exc:
        raise PipelineError(f"Visualization failed: {exc}") from exc


def _write_run_manifest(result: TopoPPIRunResult, config: TopoPPIRunConfig) -> None:
    payload = {
        "schema_version": "2.1",
        "created_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "topoppi_version": __version__,
        "input_file": str(Path(config.pdb_file).resolve()),
        "input_sha256": result.input_sha256,
        "config": config.to_dict(),
        "result": result.to_dict(),
        "environment": {
            "python": sys.version,
            "python_executable": sys.executable,
            "platform": platform.platform(),
        },
    }
    dump_json_atomic(payload, result.manifest_file)
