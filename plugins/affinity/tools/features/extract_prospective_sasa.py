#!/usr/bin/env python3
"""Compute the existing SASA-weighted ESMC readout for an unlabeled manifest."""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from topoppi_affinity.attribution_features import buried_sasa_weights, save_features, weighted_pool
from topoppi_affinity.graph_dataset import GraphDataset
from topoppi_affinity.graph_features import load_graph


def extract(manifest, graphs, sequences, output):
    columns = ['sample_id', 'source', 'complex_id', 'structure_path', 'chains_a', 'chains_b', 'mutations']
    frame = pd.read_csv(manifest, usecols=columns, keep_default_na=False)
    registry = pd.read_csv(graphs, keep_default_na=False)
    dataset = GraphDataset(frame.assign(pkd=0.), registry, sequences, max_nodes=None, target='pkd')
    vectors, coverage = {}, []
    torch.set_num_threads(1)
    for index, row in enumerate(frame.to_dict('records')):
        started = time.monotonic()
        try:
            data = dataset[index]
            graph = registry.loc[registry.source.eq(row['source']) & registry.complex_id.eq(row['complex_id'])].iloc[0]
            weights = buried_sasa_weights(row, load_graph(graph.graph_path))
            vectors[row['sample_id']] = weighted_pool(data.esm.numpy(), data.side.numpy(), weights)
            result = dict(status='ok', error='')
        except (KeyError, ValueError, OSError, IndexError) as error:
            result = dict(status='error', error=f'{type(error).__name__}: {error}')
        coverage.append(dict(sample_id=row['sample_id'], seconds=time.monotonic()-started, **result))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(coverage).to_csv(output.with_suffix('.coverage.csv'), index=False)
    if not vectors:
        raise ValueError('No SASA readout succeeded; see coverage CSV')
    width = len(next(iter(vectors.values())))
    save_features(output, frame.sample_id, [vectors.get(sample, np.full(width, np.nan)) for sample in frame.sample_id],
                  [f'weight.sasa.{i}' for i in range(width)])
    return dict(n_declared=len(frame), success=len(vectors), labels_read=False,
                protocol='Isolated-partner minus complex SASA; probe1.4A/100points; residue-weighted ESMC')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('manifest', 'graphs', 'sequences', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    print(json.dumps(extract(**vars(parser.parse_args())), indent=2))
