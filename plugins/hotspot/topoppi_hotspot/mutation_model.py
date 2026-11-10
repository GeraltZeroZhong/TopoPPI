"""Portable mutation-effect regression for alanine-scan hotspot ranking."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


class MutationSurfaceRegressor:
    """Predict ΔΔG from named physical and decoded-surface descriptors.

    The fitted arrays contain median imputation, missing-value indicators,
    standardization and an RBF kernel-ridge model. Prediction runs on CPU.
    """

    def __init__(self, metadata_path):
        self.path = Path(metadata_path)
        self.metadata = json.loads(self.path.read_text(encoding="utf-8"))
        self.feature_names = tuple(self.metadata["feature_names"])
        self.target = self.metadata["target"]
        with np.load(self.path.parent / self.metadata["weights"], allow_pickle=False) as data:
            self.arrays = {name: data[name] for name in data.files}

    def predict(self, X, feature_names):
        """Return mutation ΔΔG in kcal/mol; positive values mean weaker binding."""
        lookup = {str(name): i for i, name in enumerate(feature_names)}
        missing = set(self.feature_names) - lookup.keys()
        if missing:
            raise ValueError(f"Missing model features: {sorted(missing)[:5]}")
        values = np.asarray(X, dtype=float)[:, [lookup[name] for name in self.feature_names]]
        if np.isinf(values).any():
            raise ValueError("Mutation model inputs contain infinite values")
        absent = np.isnan(values)
        filled = np.where(absent, self.arrays["median"], values)
        indicators = absent[:, self.arrays["indicator_indices"].astype(int)]
        transformed = np.column_stack([filled, indicators])
        transformed = (transformed - self.arrays["mean"]) / self.arrays["scale"]
        support = self.arrays["training_vectors"]
        support_norm = np.einsum("ij,ij->i", support, support)
        gamma = float(self.metadata["training"]["parameters"]["gamma"])
        predictions = np.empty(len(transformed), dtype=float)
        for start in range(0, len(transformed), 256):
            block = transformed[start:start + 256]
            distance = -2.0 * block @ support.T
            distance += np.einsum("ij,ij->i", block, block)[:, None]
            distance += support_norm[None, :]
            np.maximum(distance, 0., out=distance)
            np.multiply(distance, -gamma, out=distance)
            np.exp(distance, out=distance)
            predictions[start:start + len(block)] = distance @ self.arrays["dual_coefficient"]
        return predictions
