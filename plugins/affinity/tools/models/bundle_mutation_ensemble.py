#!/usr/bin/env python3
"""Package two previously fitted mutation models into a fixed equal ensemble."""
import argparse
from pathlib import Path

from topoppi_affinity.mutation_ensemble import bundle_mutation_ensemble, bundle_scalar_prior_ensemble

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--base-artifact', type=Path, required=True)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--residual-artifact', type=Path)
group.add_argument('--prior-method', help='Frozen scalar prior, e.g. BA-Cycle-partner-MC8')
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
report = (bundle_mutation_ensemble(args.base_artifact, args.residual_artifact, args.output) if args.residual_artifact else
          bundle_scalar_prior_ensemble(args.base_artifact, args.prior_method, args.output))
print(report['combination'])
