"""Label-independent embedding alignment and feature normalization."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


def all_residue_embeddings(scan_dir, frame):
    """Map Hotspot-owned WT embeddings, including native Ala context nodes."""
    root = Path(scan_dir)
    chains = pd.read_csv(root / 'chain_sequences.csv', keep_default_na=False)
    lookup = {}
    for row in chains.itertuples(index=False):
        with np.load(root / 'sequence_embeddings' / f'{row.sequence_id}.npz', allow_pickle=False) as data:
            sequence, values = str(data['sequence']), data['embedding']
        keys = json.loads(row.residue_keys)
        if len(keys) != len(sequence) or len(values) != len(keys):
            raise ValueError('WT embedding registry is not aligned with structure residues')
        for key, wt, value in zip(keys, sequence, values, strict=True):
            lookup[row.complex_id, key] = (wt, value)
    result = []
    for row in frame.itertuples(index=False):
        wt, value = lookup[row.complex_id, row.residue_key]
        if row.wt != wt:
            raise ValueError('WT identity differs from the frozen ESMC input')
        result.append(value)
    return np.asarray(result, dtype=np.float64)


def relative_new_columns(frame, columns):
    """Two partner-only channels for each new geometric scalar."""
    grouped = frame.groupby(['complex_id', 'side'], sort=False)
    extra = {}
    for column in columns:
        values = grouped[column]
        mean, std = values.transform('mean'), values.transform('std', ddof=0)
        extra[f'new_relative.z.{column}'] = (frame[column] - mean) / std.where(std > 0, 1.)
        extra[f'new_relative.pct.{column}'] = (values.rank(method='average') - .5) / values.transform('size')
    return pd.concat([frame, pd.DataFrame(extra, index=frame.index)], axis=1), list(extra)


def center_within_queries(values, queries):
    result = np.asarray(values, dtype=float).copy()
    for query in np.unique(queries):
        mask = queries == query
        result[mask] -= result[mask].mean(axis=0)
    return result
