#!/usr/bin/env python3
"""Predict mutation ddg from unlabeled native features and a frozen scalar prior."""
import argparse
from pathlib import Path

import joblib

from topoppi_affinity.mutation_ensemble import predict_mutation_ensemble, predict_scalar_prior_ensemble

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--artifact', type=Path, required=True)
parser.add_argument('--base-features', type=Path, required=True)
parser.add_argument('--residual-features', type=Path, help='Required for a residual-component artifact')
parser.add_argument('--prior', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
bundle = joblib.load(args.artifact)
if 'residual' in bundle:
    if args.residual_features is None:
        parser.error('--residual-features is required for this residual-component artifact')
    result = predict_mutation_ensemble(bundle, args.base_features, args.residual_features, args.prior)
else:
    result = predict_scalar_prior_ensemble(bundle, args.base_features, args.prior)
args.output.parent.mkdir(parents=True, exist_ok=True)
result.to_csv(args.output, index=False)
print(f'{len(result)} predictions; {result.components_used.value_counts().to_dict()}')
