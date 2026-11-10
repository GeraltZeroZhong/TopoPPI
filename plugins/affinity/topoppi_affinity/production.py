"""Portable default affinity models with named physical inputs and fitted scaling."""
from __future__ import annotations

import json
from importlib.resources import files
from pathlib import Path

import numpy as np


class PhysicalRegressor:
    """Fitted preprocessing followed by linear or RBF kernel ridge regression."""

    def __init__(self, metadata_path):
        self.path = Path(metadata_path)
        self.metadata = json.loads(self.path.read_text(encoding="utf-8"))
        self.target = self.metadata["target"]
        self.family = self.metadata["training"]["family"]
        if self.family not in {"ridge", "rbf"}:
            raise ValueError(f"Unsupported physical regressor: {self.family}")
        self.feature_names = tuple(self.metadata["feature_names"])
        with np.load(self.path.parent / self.metadata["weights"], allow_pickle=False) as data:
            self.arrays = {name: data[name] for name in data.files}

    def predict(self, X, feature_names):
        names = list(map(str, feature_names))
        lookup = {name: i for i, name in enumerate(names)}
        missing = set(self.feature_names) - lookup.keys()
        if missing:
            raise ValueError(f"Missing model features: {sorted(missing)[:5]}")
        values = np.asarray(X, dtype=float)[:, [lookup[name] for name in self.feature_names]]
        if np.isinf(values).any():
            raise ValueError("Physical model inputs contain infinite values")
        absent = np.isnan(values)
        filled = np.where(absent, self.arrays["median"], values)
        indicators = absent[:, self.arrays["indicator_indices"].astype(int)]
        transformed = np.column_stack([filled, indicators])
        transformed = (transformed - self.arrays["mean"]) / self.arrays["scale"]
        if self.family == "ridge":
            predictions = transformed @ self.arrays["coefficient"] + self.arrays["intercept"]
        else:
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
        if self.target == "ddg" and "mutation.n_mutations" in lookup:
            predictions[np.asarray(X)[:, lookup["mutation.n_mutations"]] == 0] = 0.
        return predictions


def bundled_model(task):
    """Locate the model supplied with this installed plugin."""
    if task not in {"absolute", "mutation"}:
        raise ValueError("Default model task must be absolute or mutation")
    return Path(str(files("topoppi_affinity").joinpath("data", "production", task, "model.json")))


def load_model(path):
    path = Path(path)
    if path.suffix == ".json":
        return PhysicalRegressor(path)
    import joblib
    return joblib.load(path)


def mutation_physical_features(descriptor_bundle, mutations, temperature, *, atomic_context):
    """Return the shared 430 physical descriptors and temperature for a variant."""
    from topoppi_affinity.features import mutation_features

    values = dict(descriptor_bundle["absolute_features"])
    values.update(mutation_features(descriptor_bundle, mutations))
    atomic, quality = atomic_context.describe(mutations)
    if quality["status"] != "ok":
        raise ValueError(f"Mutation SASA context: {quality['status']}")
    values.update(atomic)
    values["temperature_ratio"] = float(temperature) / 298.15
    return values


def physical_feature_rows(prepared, descriptor_bundle, queries, model):
    """Calculate the model's exact physical scheme from the supplied complex."""
    if model.target == "pkd":
        from topoppi_representation import physical_descriptors
        values = physical_descriptors(prepared, n_points=960)
        rows = [values for _ in range(len(queries))]
    else:
        from topoppi_affinity.mutation_physics import AtomicQueryContext, compute_residue_context
        context = AtomicQueryContext(compute_residue_context(prepared))
        rows = [mutation_physical_features(descriptor_bundle, row.mutations, row.temperature_K,
                                          atomic_context=context) for row in queries.itertuples()]
    names = list(model.feature_names)
    return np.asarray([[row[name] for name in names] for row in rows], dtype=float), names
