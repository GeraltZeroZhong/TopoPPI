#!/usr/bin/env python3
"""Predict absolute affinity from unlabeled sequence/interface and conventional-descriptor NPZ inputs."""
from __future__ import annotations

import argparse
from pathlib import Path

from topoppi_affinity.legacy_kernel import predict_legacy_kernel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifact', required=True, type=Path)
    parser.add_argument('--plm-features', required=True, type=Path)
    parser.add_argument('--descriptor-features', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--temperature-k', type=float, default=298.15,
                        help='Reference temperature for thermodynamic conversion; no temperature extrapolation')
    args = parser.parse_args()
    frame = predict_legacy_kernel(args.artifact, args.plm_features, args.descriptor_features, temperature_K=args.temperature_k)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output, index=False)
    print(f'Saved {len(frame)} absolute-affinity predictions to {args.output}')


if __name__ == '__main__':
    main()
