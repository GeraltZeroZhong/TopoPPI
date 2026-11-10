"""Named feature alignment and parameter candidates for frozen representations."""
from __future__ import annotations

import numpy as np


def load_features(path, frame):
    """Read an NPZ and align its unique sample IDs to the requested table."""
    with np.load(path, allow_pickle=False) as data:
        ids = data['sample_ids'].astype(str)
        if len(set(ids)) != len(ids):
            raise ValueError('Feature NPZ contains duplicate sample IDs')
        lookup = {sample: i for i, sample in enumerate(ids)}
        return data['X'][[lookup[str(sample)] for sample in frame.sample_id]], data['feature_names']


def candidates(family, dimension):
    if family == 'trees':
        return [{'min_samples_leaf': value} for value in (2, 5)]
    if family == 'ridge':
        return [{'alpha': value} for value in (100., 1000., 10000., 100000.)]
    return [{'alpha': value, 'gamma': factor / dimension, 'gamma_factor': factor}
            for value in (.1, 1., 10.) for factor in (.1, 1.)]
