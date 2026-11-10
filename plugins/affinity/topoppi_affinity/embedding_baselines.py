"""Frozen representations, normalized kernels and reusable regression estimators."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.kernel_ridge import KernelRidge
from sklearn.linear_model import Ridge
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.preprocessing import KernelCenterer, StandardScaler

from topoppi_affinity.graph_features import FEATURE_NAMES
from topoppi_affinity.models import complex_weights

REPRESENTATIONS = ("plm_only", "plm_surface", "legacy_v1", "plm_contact")


def _mean(values, mask):
    return values[mask].mean(axis=0) if np.any(mask) else np.zeros(values.shape[1], dtype=np.float64)


def _symmetric_sides(left, right):
    return np.concatenate([left + right, np.abs(left - right)])


def pooled_features(data, *, target, use_surface):
    """Return a fixed vector and exact names; all reductions are label-free."""
    wt, mt = np.asarray(data.esm).astype(np.float64), np.asarray(data.esm_mut).astype(np.float64)
    x, x_mut = np.asarray(data.x).astype(np.float64), np.asarray(data.x_mut).astype(np.float64)
    side, contact = np.asarray(data.side), np.asarray(data.contact_mask)
    mutation = np.asarray(data.mutation_mask)
    parts, names = [], []

    def add(name, vector):
        flat = np.asarray(vector).reshape(-1)
        parts.append(flat)
        names.extend(f"{name}.{index}" for index in range(len(flat)))

    if target == "pkd":
        for scope, mask in (("global", np.ones(len(wt), dtype=bool)), ("contact", contact)):
            add(f"plm.{scope}.sum_absdiff", _symmetric_sides(_mean(wt, mask & (side == 0)), _mean(wt, mask & (side == 1))))
        if use_surface:
            surface = x[:, 30:]
            for statistic in ("mean", "sum", "std"):
                values = [getattr(surface[side == partner], statistic)(axis=0) for partner in (0, 1)]
                add(f"surface.{statistic}.sum_absdiff", _symmetric_sides(*values))
            weights = x[:, FEATURE_NAMES.index("surface.footprint_area")]
            pooled = [(wt[side == partner] * weights[side == partner, None]).sum(axis=0)
                      / max(weights[side == partner].sum(), 1.) for partner in (0, 1)]
            add("surface.area_weighted_plm.sum_absdiff", _symmetric_sides(*pooled))
    else:
        delta, context = mt - wt, (mt + wt) / 2.
        chemical_delta = x_mut[:, :26] - x[:, :26]
        add("plm.mutation_delta_sum", delta[mutation].sum(axis=0))
        add("plm.mutation_delta_context_sum", (delta[mutation] * context[mutation]).sum(axis=0))
        add("plm.contact_delta_mean", _mean(delta, contact))
        add("chemical.delta_sum", chemical_delta[mutation].sum(axis=0))
        if use_surface:
            # WT geometry is fixed; this symmetric context also preserves the
            # odd representation if a future pair supplies different geometry.
            surface = (x[:, 30:] + x_mut[:, 30:]) / 2.
            add("surface.chemical_delta_context", np.einsum("ni,nj->ij", chemical_delta[mutation], surface[mutation]))
            weights = surface[:, 1]
            add("surface.mutation_area_delta_sum", (delta[mutation] * weights[mutation, None]).sum(axis=0))
            add("surface.interface_area_delta_mean", (delta * weights[:, None]).sum(axis=0) / max(weights.sum(), 1.))
    return np.concatenate(parts).astype(np.float32), np.asarray(names)


def zero_surface_slots(matrix, feature_names):
    """Same-dimensionality control: preserve every non-surface feature exactly."""
    result = np.asarray(matrix).copy()
    result[..., np.char.startswith(np.asarray(feature_names).astype(str), "surface.")] = 0.
    return result


@dataclass
class FrozenRegressor:
    estimator: object
    scaler: object
    antisymmetric: bool
    target_center: float = 0.

    def predict(self, x):
        x = np.asarray(x, dtype=np.float64)
        values = self.scaler.transform(x) if self.scaler is not None else x
        prediction = self.estimator.predict(values)
        if self.antisymmetric:
            reverse = self.scaler.transform(-x) if self.scaler is not None else -x
            prediction = .5 * (prediction - self.estimator.predict(reverse))
        return np.asarray(prediction + self.target_center, dtype=float)


def fit_frozen_regressor(x, y, groups, *, family, params, antisymmetric=False, seed=11):
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    weights = complex_weights(groups)
    target_center = 0. if antisymmetric else float(np.average(y, weights=weights))
    y = y - target_center
    if antisymmetric:
        x, y = np.concatenate([x, -x]), np.concatenate([y, -y])
        weights = np.concatenate([weights, weights]) / 2.
    if family == "trees":
        scaler = None
        estimator = ExtraTreesRegressor(n_estimators=192, min_samples_leaf=params["min_samples_leaf"],
                                        max_features=.75, random_state=seed, n_jobs=2)
    elif family == "catboost":
        from catboost import CatBoostRegressor
        scaler = None
        estimator = CatBoostRegressor(iterations=params["iterations"], depth=params["depth"],
                                       learning_rate=.03, loss_function="RMSE", thread_count=params.get("thread_count", 2),
                                       random_seed=seed, verbose=False, allow_writing_files=False,
                                       task_type="CPU")
    else:
        scaler = StandardScaler().fit(x, sample_weight=weights)
        x = scaler.transform(x)
        estimator = (Ridge(alpha=params["alpha"]) if family in {"ridge", "ridge_extended"} else
                     KernelRidge(alpha=params["alpha"], kernel="rbf", gamma=params["gamma"]))
    estimator.fit(x, y, sample_weight=weights)
    return FrozenRegressor(estimator, scaler, antisymmetric, target_center)


def candidate_parameters(family, dimension):
    if family == "ridge":
        return [{"alpha": alpha} for alpha in (1., 10., 100.)]
    if family == "ridge_extended":
        return [{"alpha": alpha} for alpha in (1000., 10000.)]
    if family == "catboost":
        return [{"depth": depth, "iterations": iterations} for depth in (4, 6) for iterations in (400, 800)]
    if family == "rbf":
        return [{"alpha": alpha, "gamma": factor / dimension, "gamma_factor": factor}
                for alpha in (.1, 1., 10.) for factor in (.1, 1., 10.)]
    return [{"min_samples_leaf": leaf} for leaf in (2, 5)]


@dataclass
class NormalizedKernelBlock:
    """Training-fitted feature scaling, kernel centering, and unit trace/n."""
    columns: np.ndarray
    scaler: object
    training: np.ndarray
    centerer: object
    normalization: float
    gamma: float

    def transform(self, matrix):
        values = self.scaler.transform(np.asarray(matrix[:, self.columns], dtype=np.float64))
        gram = rbf_kernel(values, self.training, gamma=self.gamma)
        return self.centerer.transform(gram) / self.normalization


def fit_kernel_block(matrix, columns, sample_weight, gamma_factor):
    values = np.asarray(matrix[:, columns], dtype=np.float64)
    scaler = StandardScaler().fit(values, sample_weight=sample_weight)
    values = scaler.transform(values)
    gamma = gamma_factor / values.shape[1]
    gram = rbf_kernel(values, gamma=gamma)
    centerer = KernelCenterer().fit(gram)
    centered = centerer.transform(gram)
    normalization = max(float(np.trace(centered) / len(centered)), 1e-12)
    block = NormalizedKernelBlock(np.asarray(columns), scaler, values, centerer, normalization, gamma)
    return block, centered / normalization


@dataclass
class MultiKernelRegressor:
    """Existing convex kernel fusion plus precomputed sklearn KernelRidge.

    Mixture weights are chosen from a bounded validation grid. This does not
    claim to reproduce SimpleMKL's SVM gradient optimization algorithm.
    """
    estimator: object
    blocks: tuple
    surface_weight: float
    antisymmetric: bool
    target_center: float

    def gram(self, matrix):
        values = np.asarray(matrix, dtype=np.float64)
        result = np.zeros((len(values), len(self.blocks[0].training)))
        for weight, block in zip((1. - self.surface_weight, self.surface_weight), self.blocks, strict=True):
            if weight > 0:
                result += weight * block.transform(values)
        return result

    def predict(self, matrix):
        prediction = self.estimator.predict(self.gram(matrix))
        if self.antisymmetric:
            prediction = .5 * (prediction - self.estimator.predict(self.gram(-np.asarray(matrix))))
        return prediction + self.target_center
