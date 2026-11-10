"""Export physicochemical descriptors from selected protein partners."""

from __future__ import annotations

import argparse
import contextlib
import csv
import json
import math
import sys
from pathlib import Path


def _descriptors(args):
    output = args.output.expanduser().resolve()
    if any((output / name).exists() for name in ("descriptors.csv", "metadata.json", "residues.csv", "contacts.csv")):
        raise ValueError(f"Descriptor results already exist in {output}. Choose a new output directory.")

    from topoppi import __version__, get_interface, load_structure

    from .descriptors.physical import SCHEMA, chemical_area_descriptors, residue_physics, summarize

    interface = get_interface(load_structure(args.structure), args.chains_a, args.chains_b)
    residues, contacts = residue_physics(interface, n_points=args.n_points, probe_radius=args.probe_radius)
    values = {**summarize(residues, contacts), **chemical_area_descriptors(residues)}
    output.mkdir(parents=True, exist_ok=True)
    files = {"descriptors": "descriptors.csv", "metadata": "metadata.json"}
    with (output / files["descriptors"]).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=values)
        writer.writeheader()
        writer.writerow(values)
    for requested, name, frame in ((args.residues, "residues", residues), (args.contacts, "contacts", contacts)):
        if requested:
            files[name] = f"{name}.csv"
            frame.to_csv(output / files[name], index=False)
    metadata = {
        "format": "topoppi-physical-descriptors",
        "version": 1,
        "topoppi_version": __version__,
        "schema": SCHEMA,
        "structure": str(Path(args.structure).expanduser().resolve()),
        "partner_groups": [list(group) for group in sorted(partner.chains for partner in interface.partners)],
        "numbering": interface.numbering,
        "parameters": {"sasa_points": args.n_points, "probe_radius_angstrom": args.probe_radius,
                       "contact_distance_angstrom": 5.5},
        "n_descriptors": len(values),
        "n_residues": len(residues),
        "n_contact_pairs": len(contacts),
        "files": files,
    }
    (output / files["metadata"]).write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return {**metadata, "descriptors": values, "files": {name: str(output / file) for name, file in files.items()}}


def main(argv=None):
    parser = argparse.ArgumentParser(prog="topoppi representation", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    descriptors = commands.add_parser("descriptors", help="Export 477 physicochemical descriptors from a structure")
    descriptors.add_argument("--structure", required=True, help="Input PDB or mmCIF structure")
    descriptors.add_argument("-A", "--chains-a", required=True, help="First partner's comma-separated chain IDs")
    descriptors.add_argument("-B", "--chains-b", required=True, help="Second partner's comma-separated chain IDs")
    descriptors.add_argument("-o", "--output", "--output-dir", type=Path, required=True, help="New results directory")
    descriptors.add_argument("--residues", action="store_true", help="Also save residue SASA, chemistry and contacts")
    descriptors.add_argument("--contacts", action="store_true", help="Also save cross-partner residue contact pairs")
    descriptors.add_argument("--n-points", type=int, default=960, help="SASA sample points per atom (default: 960)")
    descriptors.add_argument("--probe-radius", type=float, default=1.4, help="SASA probe radius in angstroms (default: 1.4)")
    descriptors.add_argument("--json", action="store_true", help="Print descriptors and saved-file locations as JSON")
    args = parser.parse_args(argv)
    if args.n_points <= 0:
        parser.error("--n-points must be a positive integer")
    if not math.isfinite(args.probe_radius) or args.probe_radius <= 0:
        parser.error("--probe-radius must be a finite positive number")
    from topoppi.errors import TopoPPIError

    try:
        with contextlib.redirect_stdout(sys.stderr):
            result = _descriptors(args)
    except ModuleNotFoundError as exc:
        print(f"Descriptor dependency '{exc.name}' is unavailable. Install it with: "
              "python -m pip install 'topoppi[representation]'", file=sys.stderr)
        return 1
    except (OSError, ValueError, TopoPPIError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    if args.json:
        print(json.dumps(result, indent=2, allow_nan=False))
    else:
        print(f"Saved {result['n_descriptors']} descriptors: {result['files']['descriptors']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
