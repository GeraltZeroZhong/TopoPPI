"""Command-line interface for TopoPPI."""

from __future__ import annotations

import argparse
import json
import logging
import math
import sys
from dataclasses import asdict, replace
from difflib import get_close_matches
from pathlib import Path
from typing import Optional, Sequence

from topoppi import __version__
from topoppi.config import DEFAULT_RUN_CONFIG, TopoPPIRunConfig
from topoppi.errors import ExecutionCancelled, PluginError, TopoPPIError
from topoppi.logging_utils import setup_logging


def run_interface_mapping(*args, **kwargs):
    """Preserve the mapping entry point while loading its renderer on demand."""
    from topoppi.pipeline import run_interface_mapping as run

    return run(*args, **kwargs)


def _platform_help():
    from topoppi.plugins import get_registry

    commands = get_registry().commands()
    text = "Core commands: interface, geometry, representation, plugins, render. "
    if commands:
        text += "Installed plugin commands: " + ", ".join(commands) + ". "
    return text


def build_parser() -> argparse.ArgumentParser:
    defaults = DEFAULT_RUN_CONFIG
    parser = argparse.ArgumentParser(
        prog="topoppi",
        description="Create an annotated 2D atlas or 3D interface view from a protein complex.",
        epilog=(
            _platform_help()
            + "\n"
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
        metavar="BACKEND",
        default=defaults.interaction_source,
        help="Chemical backend name (default: prolif), or geometric for heavy-atom contacts",
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
        "--interactions", "--prolif",
        dest="prolif",
        default=argparse.SUPPRESS,
        help="Existing interaction JSON; omit to generate one beside the output image (--prolif is an alias)",
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
        help="Use distance-based interaction assignments when chemical interaction generation fails",
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
        if args.export_atlas and Path(args.export_atlas).resolve() == Path(args.output).resolve():
            raise ValueError("Choose separate image and atlas output paths.")
        annotation_file = getattr(args, "annotation_file", None)
        if args.export_atlas and annotation_file and (
            Path(args.export_atlas).expanduser().resolve() == Path(annotation_file).expanduser().resolve()
        ):
            raise ValueError("Choose separate annotation input and atlas output paths.")
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
            save_atlas(args.export_atlas, document.patches, document.visualizer, run_metadata=document.metadata)
        import matplotlib.pyplot as plt

        plt.close(figure)
        log.info("Saved %s view to %s using the stored atlas geometry.", style.get("view", "atlas"), args.output)
    except (TopoPPIError, OSError, ValueError, KeyError) as exc:
        log.error("%s", exc)
        return 1
    return 0


def _plugin_management(argv) -> int:
    from topoppi.plugins import get_registry

    parser = argparse.ArgumentParser(prog="topoppi plugins", description="Inspect optional installed plugins.")
    commands = parser.add_subparsers(dest="command", required=True)
    listing = commands.add_parser("list", help="List providers without importing plugin code")
    listing.add_argument("--json", action="store_true", help="Print machine-readable metadata")
    detail = commands.add_parser("info", help="Read one plugin's registration metadata")
    detail.add_argument("name")
    detail.add_argument("--json", action="store_true", help="Print machine-readable metadata")
    args = parser.parse_args(argv)
    registry = get_registry()
    if args.command == "list":
        values = registry.list()
        if args.json:
            print(json.dumps([asdict(value) for value in values], indent=2))
        elif not values:
            print("No optional TopoPPI plugins are installed.")
        else:
            for value in values:
                print(f"{value.name}\t{value.version}\t{value.description}")
    else:
        value = registry.info(args.name)
        if args.json:
            print(json.dumps(asdict(value), indent=2))
        else:
            print(f"{value.name} {value.version}\n{value.description}")
            print("Capabilities: " + (", ".join(value.capabilities) or "none"))
    return 0


def _core_analysis(command, argv) -> int:
    """Run structure/contact/geometry operations without loading the renderer."""
    from topoppi import compute_geometry, get_interface, load_structure, save_interface
    from topoppi.config import SurfaceConfig, TopologyConfig
    from topoppi.execution import ExecutionContext

    parser = argparse.ArgumentParser(prog=f"topoppi {command}", description=(
        "Extract protein partner contacts." if command == "interface" else
        "Compute and save native protein-interface geometry."
    ))
    parser.add_argument("structure", help="Input PDB or mmCIF structure")
    parser.add_argument("-A", "--chains-a", "--chain-a", default="A", help="First partner's comma-separated chain IDs")
    parser.add_argument("-B", "--chains-b", "--chain-b", default="B", help="Second partner's comma-separated chain IDs")
    parser.add_argument("--contact-distance", type=float, default=6., help="Heavy-atom contact cutoff in angstroms")
    parser.add_argument("-o", "--output", required=command == "geometry", help=(
        "Contact JSON destination; defaults to stdout" if command == "interface" else "Native interface NPZ destination"
    ))
    if command == "geometry":
        parser.add_argument("--res", type=float, default=SurfaceConfig().grid_resolution, help="Surface grid spacing in angstroms")
        parser.add_argument("--cutoff", type=float, default=TopologyConfig().distance_cutoff, help="Interface face cutoff in angstroms")
    args = parser.parse_args(argv)
    if args.output and Path(args.output).expanduser().resolve() == Path(args.structure).expanduser().resolve():
        raise ValueError("Choose an output path separate from the input structure.")
    context = ExecutionContext()
    with context.operation(command):
        structure = load_structure(args.structure)
        interface = get_interface(structure, args.chains_a, args.chains_b, contact_distance=args.contact_distance)
        if command == "interface":
            record = {"interface": interface.metadata(), "contacts": interface.contacts()}
            text = json.dumps(record, indent=2) + "\n"
            if args.output:
                output = Path(args.output)
                output.parent.mkdir(parents=True, exist_ok=True)
                output.write_text(text, encoding="utf-8")
            else:
                print(text, end="")
        else:
            geometry = compute_geometry(interface, surface=SurfaceConfig(grid_resolution=args.res),
                                        topology=TopologyConfig(distance_cutoff=args.cutoff))
            output = save_interface(geometry, args.output)
            print(json.dumps({"output": str(output), "interface": geometry.metadata()}, indent=2))
    return 0


def _platform_dispatch(argv):
    """Return a handled command's exit code, or None for legacy mapping args."""
    if not argv or argv[0].startswith("-"):
        return None
    command = argv[0]
    try:
        if command == "plugins":
            return _plugin_management(argv[1:])
        if command in {"interface", "geometry"}:
            return _core_analysis(command, argv[1:])
        if command == "representation":
            from topoppi_representation.cli import main as representation_main

            return representation_main(argv[1:])
        if command == "render":
            return _render_atlas(argv[1:])
        from topoppi.plugins import get_registry

        registry = get_registry()
        # Entry-point names are top-level plugin namespaces. Selecting one
        # loads only that registration module, never other installed plugins.
        if registry.available(command):
            previous_argv = sys.argv
            try:
                sys.argv = [f"topoppi {command}", *argv[1:]]
                result = registry.run_command(command, argv[1:])
            finally:
                sys.argv = previous_argv
            return int(result or 0)
        if "." not in command and "/" not in command and "\\" not in command and not Path(command).exists():
            candidates = {"plugins", "interface", "geometry", "representation", "render", *registry.commands()}
            suggestions = get_close_matches(command, sorted(candidates), n=1, cutoff=0.8)
            if suggestions:
                raise ValueError(f"Unknown command '{command}'. Did you mean 'topoppi {suggestions[0]}'? "
                                 "Use 'topoppi --help' to list commands.")
            registry.info(command)  # A missing plugin gets the actionable boundary error.
    except ExecutionCancelled as exc:
        print(str(exc), file=sys.stderr)
        return 130
    except (PluginError, TopoPPIError, OSError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return None


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    handled = _platform_dispatch(argv)
    if handled is not None:
        return handled
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
