#!/usr/bin/env python3
"""Predict mutation ddG from unlabeled prior CSV and frozen-feature NPZ inputs."""
from __future__ import annotations

import argparse
from pathlib import Path

from topoppi_affinity.residual_artifacts import predict_residual_affinity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifact', type=Path, required=True)
    parser.add_argument('--prior', type=Path, required=True)
    parser.add_argument('--features', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = predict_residual_affinity(args.artifact, args.prior, args.features, args.output)
    print(f'{result.status.eq("ok").sum()}/{len(result)} predictions written to {args.output}')


if __name__ == '__main__':
    main()
