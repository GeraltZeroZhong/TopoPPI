"""Localize the contribution of native footprints with matched permutations."""

import numpy as np


def permute_footprints(data, seed):
    """Preserve each partner/contact stratum's entire joint surface distribution."""
    result = data.clone()
    generator = np.random.default_rng(seed)
    side, contact = data.side.numpy(), data.contact_mask.numpy()
    for partner in (0, 1):
        for interface in (False, True):
            indices = np.flatnonzero((side == partner) & (contact == interface))
            permutation = generator.permutation(indices)
            result.x[indices, 30:] = data.x[permutation, 30:]
            result.x_mut[indices, 30:] = data.x_mut[permutation, 30:]
    return result
