"""Measured-state difference augmentation within identical assay conditions.

This uses the standard thermodynamic state-difference identity, related to
mutation-path augmentation. Both endpoints must be measured training variants;
no individual mutation label or intermediate-state affinity is invented.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd

from topoppi_affinity.embedding_baselines import pooled_features
from topoppi_affinity.joint_training import WT_MEASUREMENT_KEYS


def training_variant_pairs(frame, train, *, seed=11):
    """Bound augmentation to <= original train size, round-robin by assay group."""
    keys = ['source', 'complex_id', *WT_MEASUREMENT_KEYS]
    selected = frame.loc[train].sort_values('sample_id').drop_duplicates(keys + ['mutations'])
    random = np.random.default_rng(seed)
    queues = []
    for _, group in selected.groupby(keys, sort=True, dropna=False):
        variants = list(group.itertuples())
        options = list(combinations(range(len(variants)), 2))
        random.shuffle(options)
        pairs = []
        for i, j in options:
            a, b = variants[i], variants[j]
            pairs.append({'left_id': a.sample_id, 'right_id': b.sample_id, 'complex_id': a.complex_id,
                          'family_group': a.family_group, 'ddg': float(b.ddg) - float(a.ddg)})
        if pairs:
            queues.append(pairs)
    result, position = [], 0
    while len(result) < int(np.asarray(train).sum()):
        added = False
        for queue in queues:
            if position < len(queue):
                result.append(queue[position])
                added = True
                if len(result) >= int(np.asarray(train).sum()):
                    break
        if not added:
            break
        position += 1
    return pd.DataFrame(result, columns=['left_id', 'right_id', 'complex_id', 'family_group', 'ddg'])


def measured_pair_features(left, right):
    """Both observed variants live on the same WT graph and residue ordering."""
    pair = left.clone()
    pair.x, pair.x_mut = left.x_mut, right.x_mut
    pair.esm, pair.esm_mut = left.esm_mut, right.esm_mut
    pair.mutation_mask = pair.x[:, :21].argmax(dim=1) != pair.x_mut[:, :21].argmax(dim=1)
    return pooled_features(pair, target='ddg', use_surface=True)
