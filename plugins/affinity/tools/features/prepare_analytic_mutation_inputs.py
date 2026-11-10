#!/usr/bin/env python3
"""Replace native curvature pooling in unlabeled mutation inputs only."""
import argparse
from pathlib import Path

from topoppi_affinity.analytic_preprocessing import prepare_analytic_mutation_inputs

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", required=True, type=Path)
parser.add_argument("--graph-manifest", required=True, type=Path)
parser.add_argument("--features", required=True, type=Path)
parser.add_argument("--contact-features", required=True, type=Path)
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()
prepare_analytic_mutation_inputs(args.input, args.graph_manifest, args.features, args.contact_features, args.output)
