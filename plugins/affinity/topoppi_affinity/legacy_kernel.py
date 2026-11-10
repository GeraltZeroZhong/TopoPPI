"""Two-block normalized kernel regression with named sequence and structural features."""
from __future__ import annotations

import joblib
import numpy as np
import pandas as pd
from sklearn.kernel_ridge import KernelRidge
from threadpoolctl import threadpool_limits

from topoppi_affinity.embedding_baselines import MultiKernelRegressor, fit_kernel_block
from topoppi_affinity.models import complex_weights


def concatenate_blocks(plm, plm_names, legacy, legacy_names, legacy_groups, *, use_surface):
    """Keep both kernel dimensions fixed while removing surface from both blocks."""
    plm, legacy = np.array(plm, copy=True), np.array(legacy, copy=True)
    if not use_surface:
        plm[:, np.char.startswith(np.asarray(plm_names).astype(str), 'surface.')] = 0.
        legacy[:, np.asarray(legacy_groups).astype(str) == 'surface'] = 0.
    matrix = np.concatenate([plm, legacy], axis=1)
    names = np.concatenate([np.char.add('esm_block.', np.asarray(plm_names).astype(str)),
                            np.char.add('legacy_block.', np.asarray(legacy_names).astype(str))])
    second = np.arange(matrix.shape[1]) >= plm.shape[1]
    return matrix, names, second


def fit_legacy_bundle(plm, plm_names, legacy, legacy_names, legacy_groups, training, *,
                      alpha=1., gamma_factor=.1, legacy_weight=.25, use_surface=True):
    """Refit explicitly fixed parameters; no split selection occurs here."""
    if not (0. <= legacy_weight <= 1.) or alpha <= 0 or gamma_factor <= 0:
        raise ValueError("Convex kernel weight must be in [0, 1]; alpha and gamma_factor must be positive")
    matrix, names, second = concatenate_blocks(plm, plm_names, legacy, legacy_names, legacy_groups,
                                               use_surface=use_surface)
    weights = complex_weights(training.complex_id)
    labels = training.pkd.to_numpy(float)
    center = float(np.average(labels, weights=weights))
    with threadpool_limits(limits=2):
        left, left_gram = fit_kernel_block(matrix, ~second, weights, gamma_factor)
        right, right_gram = fit_kernel_block(matrix, second, weights, gamma_factor)
        gram = (1.-legacy_weight)*left_gram + legacy_weight*right_gram
        estimator = KernelRidge(alpha=alpha, kernel='precomputed').fit(gram, labels-center, sample_weight=weights)
    model = MultiKernelRegressor(estimator, (left, right), legacy_weight, False, center)
    report = {'target': 'pkd', 'method': 'Sequence and structural normalized two-kernel regression',
              'parameters': {'alpha': alpha, 'gamma_factor': gamma_factor, 'legacy_weight': legacy_weight},
              'use_surface': use_surface, 'n_train': len(training), 'n_features': matrix.shape[1],
              'n_complexes': training.complex_id.nunique(), 'n_families': training.family_group.nunique(),
              'train_sample_ids': training.sample_id.tolist(), 'train_family_groups': sorted(training.family_group.unique()),
              'kernel_normalization': 'Two training-only scalers, kernel centerers and trace/n normalizers saved in model',
              'affinity_units': 'pKd=-log10(Kd/M); deltaG=-R*T*ln(10)*pKd in kcal/mol at explicit reference temperature',
              'selection': 'Caller-specified fixed parameters',
              'surface_scope': 'Native three-dimensional surfaces and residue footprints; no UV optimization',
              'contact_control': 'When use_surface=false, zero ESM surface.* and saved legacy surface groups'}
    return {'model': model, 'feature_names': names, 'plm_feature_names': np.asarray(plm_names),
            'legacy_feature_names': np.asarray(legacy_names), 'legacy_feature_groups': np.asarray(legacy_groups),
            'report': report}


def predict_legacy_kernel(artifact, plm_features, legacy_features, *, temperature_K=298.15):
    """Two label-free NPZ inputs: sample_ids, X, feature_names, in any row order."""
    from topoppi_affinity.thermodynamics import pkd_to_dg, pkd_to_kd

    bundle = joblib.load(artifact)
    with np.load(plm_features, allow_pickle=False) as source:
        ids, plm, plm_names = source['sample_ids'].astype(str), source['X'], source['feature_names']
    with np.load(legacy_features, allow_pickle=False) as source:
        legacy_ids, legacy, legacy_names = source['sample_ids'].astype(str), source['X'], source['feature_names']
    if len(set(ids)) != len(ids) or len(set(legacy_ids)) != len(legacy_ids):
        raise ValueError('Prediction feature caches must have unique sample IDs')
    if set(ids) != set(legacy_ids):
        raise ValueError('Sequence and conventional-descriptor prediction caches must contain identical sample IDs')
    if not np.array_equal(plm_names, bundle['plm_feature_names']) or not np.array_equal(legacy_names, bundle['legacy_feature_names']):
        raise ValueError('Input feature columns differ from the fitted sequence or conventional-descriptor schema')
    lookup = {sample: index for index, sample in enumerate(legacy_ids)}
    legacy = legacy[[lookup[sample] for sample in ids]]
    matrix, names, _ = concatenate_blocks(plm, plm_names, legacy, legacy_names, bundle['legacy_feature_groups'],
                                          use_surface=bundle['report']['use_surface'])
    if not np.array_equal(names, bundle['feature_names']):
        raise ValueError('Combined feature schema differs from fitted model')
    with threadpool_limits(limits=2):
        predicted = bundle['model'].predict(matrix)
    return pd.DataFrame({'sample_id': ids, 'predicted': predicted, 'predicted_pkd': predicted,
        'predicted_kd_M': pkd_to_kd(predicted),
        'predicted_dg_at_reference_kcal_mol': pkd_to_dg(predicted, temperature_K),
        'reference_temperature_K': temperature_K})
