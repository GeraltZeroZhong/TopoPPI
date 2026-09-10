"""Command-line interface for TopoPPI."""

from __future__ import annotations

import argparse
import logging
import math
import sys
from dataclasses import replace
from pathlib import Path
from typing import Optional, Sequence

from topoppi import __version__
from topoppi.config import DEFAULT_RUN_CONFIG, TopoPPIRunConfig
from topoppi.errors import TopoPPIError
from topoppi.logging_utils import setup_logging
from topoppi.pipeline import run_interface_mapping


def build_parser() -> argparse.ArgumentParser:
    defaults = DEFAULT_RUN_CONFIG
    parser = argparse.ArgumentParser(
        prog="topoppi",
        description="Create an annotated 2D atlas or 3D interface view from a protein complex.",
        epilog=(
            "Create a map: topoppi complex.pdb -A A -B B -o interface_map.png\n"
            "Restyle an atlas: topoppi render interface.npz -o interface.svg\n"
            "View its 3D surface: topoppi render interface.npz --view surface -o surface.png\n"
            "Use topoppi render --help for saved-atlas options. Open the desktop app with topoppi-gui."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    structure = parser.add_argument_group("structure and interactions")
    structure.add_argument("pdb_file", help="Input PDB or mmCIF structure")
    structure.add_argument(
        "--interaction-source",
        choices=("prolif", "geometric"),
        default=defaults.interaction_source,
        help="Interaction partners for residue weights and annotations: ProLIF or heavy-atom contacts",
    )
    structure.add_argument(
        "-A",
        "--chain-a",
        dest="chain_a",
        default=defaults.chain_a,
        help="Chain used to build the molecular surface",
    )
    structure.add_argument(
        "-B",
        "--chain-b",
        dest="chain_b",
        default=defaults.chain_b,
        help="Partner chain used to locate the interface",
    )
    structure.add_argument(
        "--prolif",
        default=argparse.SUPPRESS,
        help="Existing ProLIF JSON; omit to generate one beside the output image",
    )
    structure.add_argument(
        "--geometric-fallback-distance",
        dest="contact_distance",
        type=float,
        default=defaults.contact_distance_angstrom,
        help="Heavy-atom contact distance in Angstroms for the geometric source or selected fallback",
    )
    structure.add_argument(
        "--residue-scope",
        choices=("interaction", "patch"),
        default=argparse.SUPPRESS,
        help="Annotation scope: interaction partners or all patch residues (footprints and 3D views default: patch)",
    )
    structure.add_argument(
        "--min-points",
        dest="min_points",
        type=int,
        default=defaults.visualization.min_points,
        help="Minimum interaction residues per visible 2D marker patch; footprints and 3D views show all patches",
    )
    structure.add_argument(
        "--uniform-residue-color",
        action="store_false",
        dest="color_by_interaction_type",
        default=argparse.SUPPRESS,
        help="Use one color for residue markers; marker colors normally show interaction types",
    )
    structure.add_argument(
        "--geometric-interaction-fallback",
        action="store_true",
        default=defaults.visualization.use_geometric_interaction_fallback,
        help="Use distance-based interaction assignments when ProLIF generation fails",
    )

    surface = parser.add_argument_group("surface generation")
    surface.add_argument(
        "--cutoff",
        type=float,
        default=defaults.topology.distance_cutoff,
        help="Chain B distance used to extract the interface patch in Angstroms",
    )
    surface.add_argument(
        "--res",
        type=float,
        default=defaults.surface.grid_resolution,
        help="Surface grid spacing in Angstroms",
    )
    surface.add_argument(
        "--sigma",
        type=float,
        default=defaults.surface.sigma,
        help="Gaussian smoothing scale in Angstroms",
    )
    surface.add_argument(
        "--surface-level", type=float, default=defaults.surface.level, help="Gaussian-density isosurface level"
    )
    surface.add_argument(
        "--surface-padding", type=float, default=defaults.surface.padding, help="Grid padding in Angstroms"
    )
    surface.add_argument(
        "--max-voxels", type=int, default=defaults.surface.max_voxels, help="Maximum dense-grid voxel count"
    )
    surface.add_argument(
        "--max-adaptive-resolution",
        type=float,
        default=defaults.surface.max_adaptive_resolution,
        help="Largest grid spacing allowed during automatic coarsening in Angstroms",
    )
    surface.add_argument(
        "--no-adaptive-resolution",
        action="store_false",
        dest="adaptive_resolution",
        default=argparse.SUPPRESS,
        help="Keep the requested grid spacing; report an error if the grid exceeds --max-voxels",
    )

    mapping = parser.add_argument_group("UV mapping and OptCuts")
    mapping.add_argument(
        "--parameterization",
        choices=("auto", "lscm", "harmonic", "slim", "spherical", "cylindrical"),
        default=defaults.parameterization.method,
        help="Initial UV parameterization",
    )
    mapping.add_argument(
        "--slim-iterations",
        type=int,
        default=defaults.parameterization.slim_iterations,
        help="SLIM iterations when --parameterization=slim",
    )
    mapping.add_argument(
        "--slim-boundary-constraint-weight",
        type=float,
        default=defaults.parameterization.slim_boundary_constraint_weight,
        help="Soft convex-boundary constraint weight for SLIM",
    )
    mapping.add_argument(
        "--optcuts-bin",
        default=defaults.optcuts.optcuts_bin,
        help="OptCuts executable path or command name",
    )
    mapping.add_argument(
        "--patch-gap",
        type=float,
        default=defaults.optcuts.patch_gap,
        help="Chart gap as a fraction of sqrt(total 3-D chart area)",
    )
    mapping.add_argument(
        "--optcuts-lambda", type=float, default=defaults.optcuts.optcuts_lambda_init, help="OptCuts initial lambda"
    )
    mapping.add_argument(
        "--optcuts-distortion-bound",
        type=float,
        default=defaults.optcuts.optcuts_distortion_bound,
        help="OptCuts distortion bound",
    )
    mapping.add_argument(
        "--optcuts-initial-cut-option",
        type=int,
        choices=(0, 1),
        default=defaults.optcuts.optcuts_initial_cut_option,
        help="OptCuts initial cut option (0=random two-edge, 1=farthest two-point)",
    )
    mapping.add_argument(
        "--no-optcuts-bijectivity",
        action="store_false",
        dest="optcuts_use_bijectivity",
        default=argparse.SUPPRESS,
        help="Disable OptCuts bijectivity enforcement",
    )
    mapping.add_argument(
        "--optcuts-initialization",
        choices=("provided", "automatic"),
        default="provided" if defaults.optcuts.use_input_uv else "automatic",
        help="Use the selected UV parameterization or OptCuts automatic initialization",
    )
    mapping.add_argument(
        "--optcuts-timeout",
        type=float,
        default=defaults.optcuts.timeout_sec,
        help="OptCuts timeout per patch in seconds",
    )
    mapping.add_argument(
        "--residue-fragmentation-weight",
        type=float,
        default=defaults.optcuts.residue_fragmentation_weight,
        help="Weight of TopoPPI's residue-footprint objective; 0 selects the geometry-only ablation",
    )
    mapping.add_argument(
        "--expected-optcuts-sha256",
        default=argparse.SUPPRESS,
        help="Expected OptCuts SHA-256 for a frozen run",
    )

    output = parser.add_argument_group("output and logging")
    output.add_argument("--output", "-o", default=defaults.output_file, help="Output PNG, TIFF, SVG, or PDF figure")
    output.add_argument("--show", action="store_true", help="Open the Matplotlib figure after saving")
    output.add_argument("--verbose", "-v", action="store_true", help="Show debug logging")
    output.add_argument("--export-atlas", metavar="FILE.npz", help="Save a self-contained atlas for later rendering")
    _add_view_options(parser)
    _add_footprint_options(parser)
    return parser


def _finite_float(value):
    try:
        number = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Use a finite number.") from exc
    if not math.isfinite(number):
        raise argparse.ArgumentTypeError("Use a finite number.")
    return number


def _positive_float(value):
    number = _finite_float(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("Use a positive finite number.")
    return number


def _add_view_options(parser):
    group = parser.add_argument_group("view and 3D camera")
    group.add_argument("--view", choices=("atlas", "surface"), default=argparse.SUPPRESS,
                       help="Render a 2D atlas or the 3D interface surface (new maps: atlas)")
    mesh = group.add_mutually_exclusive_group()
    mesh.add_argument("--show-mesh", dest="show_mesh", action="store_true", default=argparse.SUPPRESS,
                      help="Show the triangular mesh in the 3D view (new maps: shown)")
    mesh.add_argument("--no-mesh", dest="show_mesh", action="store_false", default=argparse.SUPPRESS,
                      help="Hide the triangular mesh in the 3D view")
    group.add_argument("--projection", dest="surface_projection", choices=("orthographic", "perspective"),
                       default=argparse.SUPPRESS, help="3D camera projection (new maps: orthographic)")
    group.add_argument("--elevation", dest="surface_elevation", type=_finite_float, default=argparse.SUPPRESS,
                       metavar="DEGREES", help="3D camera elevation in degrees (new maps: 73)")
    group.add_argument("--azimuth", dest="surface_azimuth", type=_finite_float, default=argparse.SUPPRESS,
                       metavar="DEGREES", help="3D camera azimuth in degrees (new maps: -90)")
    group.add_argument("--zoom", dest="surface_zoom", type=_positive_float, default=argparse.SUPPRESS,
                       metavar="FACTOR", help="Positive 3D zoom factor; larger values zoom in (new maps: 1)")


def _add_footprint_options(parser):
    group = parser.add_argument_group("residue footprints")
    group.add_argument("--map-style", choices=("markers", "footprints"), default=argparse.SUPPRESS,
                       help="Display residues as markers or filled surface footprints (new maps: markers)")
    group.add_argument("--highlight", nargs="+", default=argparse.SUPPRESS,
                       help="Residue keys to highlight, e.g. A:GLU:37 A:TYR:40 (commas also accepted)")
    group.add_argument("--annotation-file", metavar="CSV", default=argparse.SUPPRESS,
                       help="UTF-8 CSV with residue,value columns for footprint coloring; blank values appear as missing")
    group.add_argument("--annotation-label", default=argparse.SUPPRESS, help="Colorbar label, including units")
    group.add_argument("--vmin", dest="value_min", type=float, default=argparse.SUPPRESS,
                       help="Lower colorbar limit; values below it use the endpoint color and an extension marker")
    group.add_argument("--vmax", dest="value_max", type=float, default=argparse.SUPPRESS,
                       help="Upper colorbar limit; values above it use the endpoint color and an extension marker")
    group.add_argument("--labels", dest="footprint_labels", choices=("all", "highlighted", "none"),
                       default=argparse.SUPPRESS, help="Footprint labels within the annotation scope (new maps: all)")
    seams = group.add_mutually_exclusive_group()
    seams.add_argument("--show-seams", dest="show_seams", action="store_true", default=argparse.SUPPRESS,
                       help="Show optimized cut-seam outlines (new maps: shown)")
    seams.add_argument("--hide-seams", dest="show_seams", action="store_false", default=argparse.SUPPRESS,
                       help="Hide optimized cut-seam outlines")
    borders = group.add_mutually_exclusive_group()
    borders.add_argument("--show-residue-borders", dest="show_residue_borders", action="store_true",
                         default=argparse.SUPPRESS, help="Show borders between residue footprints (new maps: shown)")
    borders.add_argument("--hide-residue-borders", dest="show_residue_borders", action="store_false",
                         default=argparse.SUPPRESS, help="Hide borders between residue footprints")
    for name, meaning in (("footprint-color", "Base region color"),
                          ("highlight-color", "Highlighted region color"),
                          ("missing-color", "Color for missing numeric values")):
        group.add_argument(f"--{name}", default=argparse.SUPPRESS,
                           help=f"{meaning}; accepts a Matplotlib color name or hex value")


def _visualization_overrides(args):
    names = ("map_style", "view", "show_mesh", "surface_projection", "surface_elevation", "surface_azimuth",
             "surface_zoom", "annotation_file", "annotation_label", "value_min", "value_max", "footprint_labels",
             "show_seams", "show_residue_borders", "footprint_color", "highlight_color", "missing_color")
    changes = {name: getattr(args, name) for name in names if hasattr(args, name)}
    if hasattr(args, "highlight"):
        changes["highlight_residues"] = tuple(
            token for value in args.highlight for token in value.replace(",", " ").split()
        )
    return changes


def _render_atlas(argv):
    from topoppi.visualization.atlas_io import load_atlas, save_atlas
    from topoppi.visualization.visualizer import select_patches_for_display

    parser = argparse.ArgumentParser(
        prog="topoppi render",
        description="Restyle and export an atlas using its saved geometry and UV coordinates.",
        epilog=("Unspecified settings retain their saved values. Example: topoppi render interface.npz "
                "--map-style footprints --highlight A:GLU:37 -o interface.svg"),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("atlas_file", help="Atlas NPZ saved by TopoPPI")
    parser.add_argument("-o", "--output", required=True, help="PNG, TIFF, SVG, or PDF image")
    parser.add_argument("--export-atlas", help="Save the updated atlas and style")
    parser.add_argument("--clear-annotations", action="store_true", help="Remove saved external values")
    parser.add_argument("--residue-scope", choices=("interaction", "patch"), default=argparse.SUPPRESS,
                        help="Annotate interaction partners or all patch residues")
    parser.add_argument("--verbose", "-v", action="store_true", help="Show debug logging")
    _add_view_options(parser)
    _add_footprint_options(parser)
    args = parser.parse_args(argv)
    setup_logging(args.verbose)
    log = logging.getLogger("topoppi.cli")
    try:
        from topoppi.pipeline import _prepare_output_path

        if Path(args.output).resolve() == Path(args.atlas_file).resolve():
            raise ValueError("Choose an image output path separate from the saved atlas.")
        _prepare_output_path(args.output)
        document = load_atlas(args.atlas_file)
        style = dict(document.style)
        changes = _visualization_overrides(args)
        if args.clear_annotations or "annotation_file" in changes:
            style.pop("annotation_values", None)
            style["annotation_file"] = ""
        style.update(changes)
        if any(name in changes for name in ("surface_projection", "surface_elevation", "surface_azimuth", "surface_zoom")):
            style.pop("surface_camera", None)
        if hasattr(args, "residue_scope"):
            style["residue_scope"] = args.residue_scope
        elif (changes.get("map_style") == "footprints" and document.style.get("map_style") != "footprints") or (
            changes.get("view") == "surface" and document.style.get("view", "atlas") != "surface"
        ):
            style["residue_scope"] = "patch"
        displayed_patches, _counts = select_patches_for_display(
            document.patches, document.visualizer, map_style=style.get("map_style", "markers"),
            view=style.get("view", "atlas"),
            min_points=style.get("min_points", document.visualizer.config.min_points),
        )
        figure = document.visualizer.plot_patches(displayed_patches, output_file=args.output, show=False,
                                                  style_config=style)
        if args.export_atlas:
            if Path(args.export_atlas).resolve() == Path(args.output).resolve():
                raise ValueError("Choose separate image and atlas output paths.")
            save_atlas(args.export_atlas, document.patches, document.visualizer, run_metadata=document.metadata)
        import matplotlib.pyplot as plt

        plt.close(figure)
        log.info("Saved %s view to %s using the stored atlas geometry.", style.get("view", "atlas"), args.output)
    except (TopoPPIError, OSError, ValueError, KeyError) as exc:
        log.error("%s", exc)
        return 1
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "render":
        return _render_atlas(argv[1:])
    parser = build_parser()
    args = parser.parse_args(argv)
    setup_logging(args.verbose)
    log = logging.getLogger("topoppi.cli")

    config = TopoPPIRunConfig(
        pdb_file=args.pdb_file,
        chain_a=args.chain_a,
        chain_b=args.chain_b,
        output_file=args.output,
        prolif_file=getattr(args, "prolif", None),
        interaction_source=args.interaction_source,
        atlas_output=args.export_atlas,
        contact_distance_angstrom=args.contact_distance,
        surface=replace(
            DEFAULT_RUN_CONFIG.surface,
            grid_resolution=args.res,
            sigma=args.sigma,
            level=args.surface_level,
            padding=args.surface_padding,
            max_voxels=args.max_voxels,
            adaptive_resolution=getattr(
                args,
                "adaptive_resolution",
                DEFAULT_RUN_CONFIG.surface.adaptive_resolution,
            ),
            max_adaptive_resolution=args.max_adaptive_resolution,
        ),
        topology=replace(DEFAULT_RUN_CONFIG.topology, distance_cutoff=args.cutoff),
        parameterization=replace(
            DEFAULT_RUN_CONFIG.parameterization,
            method=args.parameterization,
            slim_iterations=args.slim_iterations,
            slim_boundary_constraint_weight=args.slim_boundary_constraint_weight,
        ),
        optcuts=replace(
            DEFAULT_RUN_CONFIG.optcuts,
            optcuts_bin=args.optcuts_bin,
            patch_gap=args.patch_gap,
            optcuts_lambda_init=args.optcuts_lambda,
            optcuts_distortion_bound=args.optcuts_distortion_bound,
            optcuts_initial_cut_option=args.optcuts_initial_cut_option,
            optcuts_use_bijectivity=getattr(
                args,
                "optcuts_use_bijectivity",
                DEFAULT_RUN_CONFIG.optcuts.optcuts_use_bijectivity,
            ),
            use_input_uv=args.optcuts_initialization == "provided",
            residue_fragmentation_weight=args.residue_fragmentation_weight,
            timeout_sec=args.optcuts_timeout,
            expected_binary_sha256=getattr(
                args,
                "expected_optcuts_sha256",
                DEFAULT_RUN_CONFIG.optcuts.expected_binary_sha256,
            ),
        ).for_headless(),
        visualization=replace(
            DEFAULT_RUN_CONFIG.visualization,
            show_plot=args.show,
            min_points=args.min_points,
            residue_scope=getattr(
                args, "residue_scope", "patch" if (
                    getattr(args, "map_style", "markers") == "footprints" or getattr(args, "view", "atlas") == "surface"
                ) else DEFAULT_RUN_CONFIG.visualization.residue_scope,
            ),
            color_by_interaction_type=getattr(
                args,
                "color_by_interaction_type",
                DEFAULT_RUN_CONFIG.visualization.color_by_interaction_type,
            ),
            use_geometric_interaction_fallback=args.geometric_interaction_fallback or args.interaction_source == "geometric",
            **_visualization_overrides(args),
        ),
    )
    try:
        run_interface_mapping(config, log=log)
    except TopoPPIError as exc:
        log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
