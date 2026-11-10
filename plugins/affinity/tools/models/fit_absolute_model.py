#!/usr/bin/env python3
"""Fit a two-kernel affinity model from an explicit training table and features."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from topoppi_affinity.frozen_probes import load_features
from topoppi_affinity.legacy_kernel import fit_legacy_bundle


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--training-table', type=Path, required=True,
                        help='CSV with sample_id, complex_id, family_group and pkd')
    parser.add_argument('--sequence-features', type=Path, required=True)
    parser.add_argument('--descriptor-features', type=Path, required=True,
                        help='NPZ with sample_ids, X, feature_names and feature_groups')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--alpha', type=float, default=1.)
    parser.add_argument('--gamma-factor', type=float, default=.1)
    parser.add_argument('--descriptor-weight', type=float, default=.25)
    parser.add_argument('--representation', choices=('surface', 'contact'), default='surface')
    args = parser.parse_args()
    training = pd.read_csv(args.training_table, keep_default_na=False, dtype={'sample_id': str})
    if training.sample_id.duplicated().any():
        raise ValueError('Training sample IDs must be unique')
    sequence, sequence_names = load_features(args.sequence_features, training)
    descriptors, descriptor_names = load_features(args.descriptor_features, training)
    with np.load(args.descriptor_features, allow_pickle=False) as data:
        descriptor_groups = data['feature_groups']
    bundle = fit_legacy_bundle(sequence, sequence_names, descriptors, descriptor_names,
        descriptor_groups, training, alpha=args.alpha, gamma_factor=args.gamma_factor,
        legacy_weight=args.descriptor_weight, use_surface=args.representation == 'surface')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(bundle, args.output)
    args.output.with_suffix('.json').write_text(json.dumps(bundle['report'], indent=2) + '\n')
    print(json.dumps({'artifact': str(args.output), 'n_train': bundle['report']['n_train'],
                      'parameters': bundle['report']['parameters']}))


if __name__ == '__main__':
    main()
