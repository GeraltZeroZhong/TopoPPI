#!/usr/bin/env python3
"""Pool existing ESMC/TopoPPI inputs for an explicit, possibly unlabeled manifest."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from threadpoolctl import threadpool_limits

from topoppi_affinity.embedding_baselines import pooled_features, zero_surface_slots
from topoppi_affinity.graph_dataset import GraphDataset

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--graphs', type=Path, required=True, help='Explicit graph manifest CSV')
parser.add_argument('--sequences', type=Path, required=True, help='Sequence embedding registry directory')
parser.add_argument('--manifest', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--target', choices=['pkd', 'ddg'], required=True)
args = parser.parse_args()
frame = pd.read_csv(args.manifest, keep_default_na=False)
frame[args.target] = 0.  # GraphDataset's y field is unused by label-free pooling.
args.output.mkdir(parents=True, exist_ok=True)
torch.set_num_threads(1)
dataset = GraphDataset(frame, args.graphs, args.sequences, max_nodes=None, target=args.target)
values, names = {'plm_only': [], 'plm_surface': []}, {}
with threadpool_limits(limits=1):
    for index in range(len(dataset)):
        data = dataset[index]
        for name in values:
            vector, names[name] = pooled_features(data, target=args.target, use_surface=name == 'plm_surface')
            values[name].append(vector)
        if (index + 1) % 100 == 0:
            print(f'Pooled {index+1}/{len(dataset)}', flush=True)
for name, vectors in values.items():
    np.savez_compressed(args.output / f'{args.target}_{name}.npz', X=np.asarray(vectors),
                        sample_ids=frame.sample_id.to_numpy(dtype=str), feature_names=names[name])
np.savez_compressed(args.output / f'{args.target}_plm_contact.npz',
                    X=zero_surface_slots(np.asarray(values['plm_surface']), names['plm_surface']),
                    sample_ids=frame.sample_id.to_numpy(dtype=str), feature_names=names['plm_surface'])
(args.output / f'{args.target}_pooling.json').write_text(json.dumps({
    'manifest': str(args.manifest), 'n_samples': len(frame), 'labels_used': False,
    'encoder': 'frozen ESMC 600M', 'max_nodes': None,
    'same_representation_as': 'embedding_baselines.pooled_features'}, indent=2) + '\n')
