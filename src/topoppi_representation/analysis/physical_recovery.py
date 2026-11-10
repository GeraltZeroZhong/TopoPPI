"""Optional linear and RBF physical-target recovery using training-only scaling."""

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.kernel_ridge import KernelRidge
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def family_weights(groups):
    groups = np.asarray(groups)
    _, index, counts = np.unique(groups, return_inverse=True, return_counts=True)
    weight = 1.0 / counts[index]
    return weight * len(weight) / weight.sum()


@dataclass
class PhysicalTargetReader:
    preprocessing: object
    o_models: list
    x_model: object
    x_training_mean: np.ndarray
    specification: dict

    def predict(self, features):
        values = self.preprocessing.transform(features)
        result = np.full((len(values), 13), np.nan)
        for side, model in enumerate(self.o_models):
            if model is not None:
                result[:, side] = np.clip(model.predict(values), 0.0, 1.0)
        if self.x_model is not None:
            distribution = np.maximum(self.x_model.predict(values), 0.0)
            total = distribution.sum(axis=1)
            positive = total > 0
            distribution[positive] /= total[positive, None]
            distribution[~positive] = self.x_training_mean
            result[:, 2:] = distribution
        return result


def fit_physical_reader(features, samples, specification):
    """Fit one bounded-output readout on fitting groups only."""
    pre = make_pipeline(
        SimpleImputer(strategy="median", add_indicator=True, keep_empty_features=True), StandardScaler()
    )
    z = pre.fit_transform(features)
    targets = np.array([s.target for s in samples], float)
    masks = np.array([s.target_mask for s in samples], bool)
    groups = np.array([s.family_id for s in samples])

    def fit(take, y):
        if not np.any(take):
            return None
        model = (
            Ridge(alpha=specification["alpha"])
            if specification["kind"] == "ridge"
            else KernelRidge(alpha=specification["alpha"], kernel="rbf", gamma=specification["gamma"])
        )
        model.fit(z[take], y[take], sample_weight=family_weights(groups[take]))
        return model

    o = [fit(masks[:, side], targets[:, side]) for side in (0, 1)]
    x_valid = masks[:, 2:].all(axis=1)
    x = fit(x_valid, targets[:, 2:])
    fallback = np.full(11, np.nan)
    if np.any(x_valid):
        fallback = np.average(targets[x_valid, 2:], axis=0, weights=family_weights(groups[x_valid]))
        fallback = np.maximum(fallback, 0.0)
        fallback /= fallback.sum()
    return PhysicalTargetReader(pre, o, x, fallback, specification)


def validation_loss(samples, predictions):
    """Family-equal O MSE / half X-L2, with valid objectives equal per case."""
    rows = []
    for sample, predicted in zip(samples, predictions, strict=False):
        parts = []
        o = sample.target_mask[:2] & np.isfinite(predicted[:2])
        if o.any():
            parts.append(float(np.square(predicted[:2] - sample.target[:2])[o].mean()))
        if sample.target_mask[2:].all() and np.isfinite(predicted[2:]).all():
            parts.append(float(np.square(predicted[2:] - sample.target[2:]).sum() / 2))
        if parts:
            rows.append({"family": sample.family_id, "loss": np.mean(parts)})
    if not rows:
        return float("inf")
    return float(pd.DataFrame(rows).groupby("family").loss.mean().mean())
