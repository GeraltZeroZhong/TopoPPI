"""State-pair kernel regression for absolute and mutation observations."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import solve
from scipy.sparse import csr_matrix
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.preprocessing import StandardScaler

from topoppi_affinity.models import complex_weights
from topoppi_affinity.thermodynamics import R_KCAL


@dataclass
class StateKernelRegressor:
    scaler: object
    training_states: np.ndarray
    coefficients: np.ndarray
    mean: float
    gamma: float

    def predict(self, states):
        return rbf_kernel(self.scaler.transform(states), self.training_states, gamma=self.gamma) @ self.coefficients + self.mean

    def difference(self, wild_type, mutant, temperature_K):
        return -R_KCAL * np.log(10.) * np.asarray(temperature_K) * (self.predict(mutant) - self.predict(wild_type))


def fit_state_kernel(states, operator, labels, weights, *, alpha, gamma_factor):
    """Rows of L select absolute states or -RT ln(10) times MT-WT."""
    scaler = StandardScaler().fit(states)
    x = scaler.transform(states)
    gamma = gamma_factor / x.shape[1]
    kernel = rbf_kernel(x, gamma=gamma)
    operator = csr_matrix(operator)
    row_sum = np.asarray(operator.sum(axis=1)).reshape(-1)
    absolute = np.isclose(row_sum, 1.)
    mean = float(np.average(np.asarray(labels)[absolute], weights=np.asarray(weights)[absolute])) if absolute.any() else 0.
    covariance = np.asarray((operator @ kernel) @ operator.T)
    covariance.flat[::len(covariance) + 1] += alpha / np.asarray(weights)
    dual = solve(covariance, np.asarray(labels) - row_sum * mean, assume_a='pos')
    return StateKernelRegressor(scaler, x, np.asarray(operator.T @ dual).reshape(-1), mean, gamma)


def wt_key(row):
    return f'{row.source}|{row.complex_id}|WT'


def mt_key(row):
    return f'{row.source}|{row.complex_id}|{row.mutations}'


def make_observations(absolute, mutation, states, *, balance_views=True):
    """Return only training states and the sparse linear observation operator."""
    keys = list(dict.fromkeys([wt_key(row) for row in absolute.itertuples()] +
                             [key for row in mutation.itertuples() for key in (wt_key(row), mt_key(row))]))
    lookup = {key: i for i, key in enumerate(keys)}
    entries, columns, rows, labels, weights = [], [], [], [], []
    primary_size = max(len(absolute), len(mutation), 1)
    for frame, target in [(absolute, 'pkd'), (mutation, 'ddg')]:
        if not len(frame):
            continue
        group_weight = complex_weights(frame.complex_id)
        if balance_views:
            group_weight = group_weight * primary_size / len(frame)
        for row, weight in zip(frame.itertuples(), group_weight, strict=True):
            observation = len(labels)
            labels.append(getattr(row, target))
            weights.append(weight)
            if target == 'pkd':
                rows.append(observation)
                columns.append(lookup[wt_key(row)])
                entries.append(1.)
            else:
                factor = R_KCAL * np.log(10.) * row.temperature_K
                rows.extend([observation, observation])
                columns.extend([lookup[wt_key(row)], lookup[mt_key(row)]])
                entries.extend([factor, -factor])
    operator = csr_matrix((entries, (rows, columns)), shape=(len(labels), len(keys)))
    return np.asarray([states[key] for key in keys]), operator, np.asarray(labels), np.asarray(weights)
