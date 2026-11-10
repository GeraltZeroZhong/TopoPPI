#!/usr/bin/env python3
"""Fit a scalar-prior residual model using grouped training-only selection."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from topoppi_affinity.residual_artifacts import fit_residual_model


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--training-table', type=Path, required=True)
    parser.add_argument('--features', type=Path, required=True)
    parser.add_argument('--prior', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--task', choices=('mutation_within', 'mutation_family'), required=True)
    parser.add_argument('--representation', choices=('surface', 'contact'), default='surface')
    parser.add_argument('--seed', type=int, default=11)
    args = parser.parse_args()
    report = fit_residual_model(args.training_table, args.features, args.prior, args.output,
        task=args.task, representation=args.representation, seed=args.seed)
    print(json.dumps({'artifact': str(args.output), 'n_train': report['n_train'],
                      'parameters': report['parameters'], 'origin_slope': report['origin_slope']}))


if __name__ == '__main__':
    main()
