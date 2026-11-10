#!/usr/bin/env python3
"""Replace pooled curvature fields using explicitly supplied analytic graphs."""
from __future__ import annotations

import argparse
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from topoppi_affinity.analytic_curvature import update_pooled_curvature
from topoppi_affinity.frozen_probes import load_features
from topoppi_affinity.graph_features import load_graph


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True,
                        help='CSV with sample_id, source, complex_id and mutations')
    parser.add_argument('--graphs', type=Path, required=True,
                        help='CSV with source, complex_id and graph_path')
    parser.add_argument('--features', type=Path, required=True,
                        help='NPZ containing the original named pooled features')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--target', choices=['pkd', 'ddg'], required=True)
    args = parser.parse_args()
    frame = pd.read_csv(args.manifest, keep_default_na=False,
                        dtype={'sample_id': str, 'source': str, 'complex_id': str})
    registry = pd.read_csv(args.graphs, keep_default_na=False,
                           dtype={'source': str, 'complex_id': str})
    paths = registry.set_index(['source', 'complex_id']).graph_path.to_dict()
    matrix, names = load_features(args.features, frame)

    @lru_cache(maxsize=16)
    def graph(key):
        return load_graph(paths[key])

    corrected, ids = [], []
    for index, row in enumerate(frame.itertuples()):
        corrected.append(update_pooled_curvature(matrix[index], names,
            graph((row.source, row.complex_id)), target=args.target,
            mutations=getattr(row, 'mutations', '')))
        ids.append(row.sample_id)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('wb') as stream:
        np.savez_compressed(stream, X=np.asarray(corrected), sample_ids=np.asarray(ids, dtype=str),
                            feature_names=names)


if __name__ == '__main__':
    main()
