"""Group-DRO exponentiated family risks with exact weighted Ridge updates.

The adversarial update follows Sagawa et al. (ICLR 2020); replacing stochastic
neural updates with a convex weighted linear minimization is an adaptation.
Only development rows enter these runs. Family risk gives each complex equal
weight within its family, independent of its number of measured mutations.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

from topoppi_affinity.evaluation import regression_metrics
from topoppi_affinity.models import complex_weights


def family_risk_matrix(complex_ids, family_ids):
    """Rows average observations within complex, then complexes within family."""
    frame = pd.DataFrame({'complex': complex_ids, 'family': family_ids})
    if frame.groupby('complex').family.nunique().max() != 1:
        raise ValueError('Each complex must belong to one family')
    families = np.unique(frame.family.astype(str))
    matrix = np.zeros((len(families), len(frame)), dtype=float)
    for index, family in enumerate(families):
        mask = frame.family.astype(str).eq(family).to_numpy()
        weights = complex_weights(frame.loc[mask, 'complex'])
        matrix[index, mask] = weights / weights.sum()
    return families, matrix


def adversarial_update(probabilities, risks, eta):
    """Author q_g <- q_g exp(eta R_g), normalized, using training risk only."""
    logits = np.log(probabilities) + eta * np.asarray(risks)
    weights = np.exp(logits - logits.max())
    return weights / weights.sum()


@dataclass
class WeightedLinearRegressor:
    estimator: object
    scaler: object
    target_center: float
    target_scale: float
    antisymmetric: bool

    def predict(self, matrix):
        transformed = self.scaler.transform(np.asarray(matrix, float))
        result = self.estimator.predict(transformed)
        if self.antisymmetric:
            result = (result - self.estimator.predict(self.scaler.transform(-np.asarray(matrix, float)))) / 2.
        return result * self.target_scale + self.target_center


class FixedTrainingCoordinates:
    """A shared train-only coordinate system for all weight objectives."""
    def __init__(self, x, y, complex_ids, antisymmetric):
        x, y = np.asarray(x, float), np.asarray(y, float)
        weights = complex_weights(complex_ids)
        self.antisymmetric = antisymmetric
        self.center = 0. if antisymmetric else float(np.average(y, weights=weights))
        self.scale = max(float(np.sqrt(np.average((y-self.center)**2, weights=weights))), 1e-12)
        self.y = (y-self.center) / self.scale
        self.scaler = StandardScaler()
        if antisymmetric:
            self.scaler.fit(np.concatenate([x, -x]), sample_weight=np.concatenate([weights, weights]) / 2.)
        else:
            self.scaler.fit(x, sample_weight=weights)
        self.x = self.scaler.transform(x)

    def fit(self, row_weights, alpha):
        # Odd linear fitting equals explicit +/- augmentation with half weights:
        # centered inputs, zero intercept, identical normal equations.
        estimator = Ridge(alpha=alpha, fit_intercept=not self.antisymmetric)
        estimator.fit(self.x, self.y, sample_weight=row_weights)
        return WeightedLinearRegressor(estimator, self.scaler, self.center, self.scale, self.antisymmetric)


def training_path(coordinates, complex_ids, family_ids, *, alpha, eta=.01, steps=20):
    """Yield a path depending solely on training labels/features/groups."""
    families, risk_map = family_risk_matrix(complex_ids, family_ids)
    probabilities = np.full(len(families), 1. / len(families))
    for iteration in range(steps + 1):
        weights = len(complex_ids) * (probabilities @ risk_map)
        model = coordinates.fit(weights, alpha)
        error = coordinates.y - model.estimator.predict(coordinates.x)
        risks = risk_map @ (error ** 2)
        yield iteration, model, families, probabilities.copy(), risks
        probabilities = adversarial_update(probabilities, risks, eta)


def family_metrics(observed, predicted, complex_ids, family_ids):
    _, matrix = family_risk_matrix(complex_ids, family_ids)
    error = np.asarray(predicted) - np.asarray(observed)
    family_mae, family_mse = matrix @ abs(error), matrix @ error**2
    return {**regression_metrics(observed, predicted, complex_ids),
            'family_macro_mae': float(family_mae.mean()), 'family_macro_mse': float(family_mse.mean()),
            'worst_family_mae': float(family_mae.max()), 'worst_family_mse': float(family_mse.max()),
            'n_families': len(matrix)}
