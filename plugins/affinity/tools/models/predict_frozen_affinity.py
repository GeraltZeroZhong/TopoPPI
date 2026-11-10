#!/usr/bin/env python3
"""Predict pKd/Kd/reference-temperature ΔG or mutation ΔΔG from unlabeled features."""
import argparse
from pathlib import Path

from topoppi_affinity.artifacts import predict_frozen

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--model', type=Path, required=True)
parser.add_argument('--features', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--reference-temperature-K', type=float, default=298.15)
args = parser.parse_args()
prediction = predict_frozen(args.model, args.features, temperature_K=args.reference_temperature_K)
args.output.parent.mkdir(parents=True, exist_ok=True)
prediction.to_csv(args.output, index=False)
print(f'Saved {len(prediction)} predictions to {args.output}')
