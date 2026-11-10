#!/usr/bin/env python3
"""Fit a named-feature regressor with a caller-selected configuration."""
import argparse
import json
from pathlib import Path

from topoppi_affinity.artifacts import fit_frozen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--training-table', type=Path, required=True)
    parser.add_argument('--features', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--target', choices=['pkd', 'ddg'], required=True)
    parser.add_argument('--representation', required=True)
    parser.add_argument('--family', choices=['ridge', 'rbf', 'trees'], required=True)
    parser.add_argument('--parameters', type=json.loads, required=True)
    parser.add_argument('--antisymmetric', action='store_true',
                        help='Enforce f(-x)=-f(x) for signed mutation features')
    args = parser.parse_args()
    report = fit_frozen(args.training_table, args.features, args.output, target=args.target,
        representation=args.representation, family=args.family, parameters=args.parameters,
        antisymmetric=args.antisymmetric)
    print(json.dumps({k: v for k, v in report.items() if not k.startswith('train_')}, indent=2))


if __name__ == '__main__':
    main()
