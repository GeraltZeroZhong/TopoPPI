"""Chemistry backend discovery and explicit interaction detection."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def main(argv=None):
    parser = argparse.ArgumentParser(prog="topoppi chemistry", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("backends", help="List installed perception backends")
    compute = commands.add_parser("detect", help="Detect interactions between selected chains")
    compute.add_argument("structure")
    compute.add_argument("--chain-a", required=True)
    compute.add_argument("--chain-b", required=True)
    compute.add_argument("--backend", default="prolif")
    compute.add_argument("--output", type=Path, help="Write canonical interaction JSON; default: standard output")
    args = parser.parse_args(argv)
    from topoppi_chemistry import api

    if args.command == "backends":
        for name in api.backends():
            print(name)
        return 0
    from topoppi import get_interface, load_structure
    from topoppi.errors import TopoPPIError
    from topoppi.json_utils import dump_json_atomic

    try:
        if args.output and args.output.expanduser().resolve() == Path(args.structure).expanduser().resolve():
            raise ValueError("Choose an output path separate from the input structure.")
        result = api.detect(get_interface(load_structure(args.structure), args.chain_a, args.chain_b), backend=args.backend)
        if args.output:
            dump_json_atomic(result.to_document(), args.output)
        else:
            print(json.dumps(result.to_document(), indent=2))
    except (OSError, ValueError, RuntimeError, TopoPPIError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
