"""Rank interface residues with the official model or a supplied fitted model."""
from __future__ import annotations

import argparse
import contextlib
import json
import sys
from pathlib import Path


def main(argv=None):
    prog = "topoppi hotspot" if Path(sys.argv[0]).name == "topoppi hotspot" else "topoppi-hotspot"
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    inference = sub.add_parser("rank-complex", help="Rank one complete WT partner using a frozen model")
    inference.add_argument("--structure", required=True, help="WT structure with author chain/residue numbering")
    inference.add_argument("--chains-a", required=True, help="First complete partner's comma-separated chain IDs")
    inference.add_argument("--chains-b", required=True, help="Second complete partner's comma-separated chain IDs")
    inference.add_argument("--partner", required=True, help="Exactly one complete input partner group to rank")
    inference.add_argument("--model-dir", help="Optional custom model folder; defaults to the bundled official model")
    inference.add_argument("--output-dir", "--output", required=True, help="New results directory for rankings and model metadata")
    inference.add_argument("--json", action="store_true", help="Print a JSON result summary; all progress goes to stderr")
    inference.add_argument("--optcuts-bin", help="Residue-aware OptCuts executable for the default UV model")
    inference.add_argument("--embedding-registry", help="Existing WT chain/embedding registry, when needed")
    inference.add_argument("--esmc-checkpoint", help="Explicit local ESMC checkpoint for required WT embeddings")
    inference.add_argument("--device", choices=("cpu", "cuda"), default="cpu", help="Device for fresh ESMC inference")
    args = parser.parse_args(argv)
    from topoppi.errors import TopoPPIError

    output = Path(args.output_dir).expanduser().resolve()
    try:
        if any((output / name).exists() for name in (
            "rankings.csv", "candidate_coverage.csv", "features.csv", "representation.bin", "interface.npz", "metadata.json",
        )):
            raise ValueError(f"Ranking results already exist in {output}. Choose a new output directory.")
        with contextlib.redirect_stdout(sys.stderr):
            from .inference import rank_complex

            result, metadata = rank_complex(args.structure, args.chains_a, args.chains_b, args.partner,
                                           args.model_dir, args.output_dir, embedding_registry=args.embedding_registry,
                                           esmc_checkpoint=args.esmc_checkpoint, device=args.device,
                                           optcuts_bin=args.optcuts_bin)
    except KeyboardInterrupt:
        print("Ranking cancelled.", file=sys.stderr)
        return 130
    except (TopoPPIError, OSError, ValueError, RuntimeError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    if args.json:
        summary = {
            "analysis": str(output / "metadata.json"), "prediction_status": "ok",
            "model": metadata.get("model_scope"), "model_id": metadata.get("model_id"),
            "model_version": metadata.get("model_version"),
            "partners": metadata.get("partner_groups"), "selected_partner": metadata.get("selected_partner"),
            "score_type": metadata.get("score_type"), "n_ranked_candidates": len(result),
            "target": metadata.get("target"), "units": metadata.get("units"),
            "temperature_K": metadata.get("temperature_K"),
            "rankings": result.astype(object).where(result.notna(), None).to_dict("records"),
            "files": {name: str(output / filename) for name, filename in (
                ("rankings", "rankings.csv"), ("candidate_coverage", "candidate_coverage.csv"),
                ("metadata", "metadata.json"))},
        }
        summary["files"].update({name: str(output / filename)
                                 for name, filename in metadata.get("files", {}).items()})
        print(json.dumps(summary, indent=2, allow_nan=False))
    else:
        print(f"Ranked {len(result)} standard non-Ala interface residues; results: {output / 'rankings.csv'}")
    return 0
