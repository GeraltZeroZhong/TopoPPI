"""Export the bundled structure example or an independently generated peptide pair."""

from __future__ import annotations

import argparse
import sys
from importlib import resources
from pathlib import Path


def copy_example(destination: str | Path, *, synthetic: bool = False) -> Path:
    """Copy 1BVK (chains B/C), or the synthetic peptide pair (chains A/B)."""
    folder = Path(destination).expanduser().resolve()
    filename = "synthetic_complex.pdb" if synthetic else "1bvk.pdb"
    source = resources.files("topoppi.assets").joinpath("examples", filename)
    data = source.read_bytes()
    folder.mkdir(parents=True, exist_ok=True)
    target = folder / filename
    if target.exists():
        if target.read_bytes() != data:
            raise FileExistsError(f"{target} already contains a different structure. Choose another output folder.")
    else:
        target.write_bytes(data)
    return target


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", help="Folder for the example structure and your analysis outputs")
    parser.add_argument("--synthetic", action="store_true", help="Use an independently generated peptide pair (chains A/B)")
    args = parser.parse_args()
    try:
        path = copy_example(args.output_dir, synthetic=args.synthetic)
    except OSError as exc:
        print(f"topoppi.examples: {exc}", file=sys.stderr)
        return 1
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
