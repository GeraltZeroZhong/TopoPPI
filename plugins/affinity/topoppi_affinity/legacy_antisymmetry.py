"""Reverse augmentation and odd projection of the existing legacy descriptors.

This is conditional antisymmetry in a fixed reference structural context. It
does not claim relaxed mutant structures or a globally transitive state energy.
"""
from dataclasses import dataclass

import numpy as np


def reverse_legacy(matrix, names):
    names = np.asarray(names).astype(str)
    lookup = {name: i for i, name in enumerate(names)}
    result = np.array(matrix, copy=True)
    for name, index in lookup.items():
        if name.startswith('mutation.wt_'):
            other = lookup[name.replace('mutation.wt_', 'mutation.mt_', 1)]
            result[..., index], result[..., other] = matrix[..., other], matrix[..., index]
        elif name.startswith(('mutation.delta_', 'mutation.contact_delta_', 'mutation.surface_delta_')):
            result[..., index] *= -1
    return result


@dataclass
class ConditionalOddTree:
    estimator: object
    feature_names: np.ndarray

    def predict(self, matrix):
        return .5*(self.estimator.predict(matrix) - self.estimator.predict(reverse_legacy(matrix, self.feature_names)))
