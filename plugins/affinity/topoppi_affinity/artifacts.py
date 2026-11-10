"""Fit named-feature regressors and apply their saved preprocessing."""
from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from topoppi_affinity.embedding_baselines import fit_frozen_regressor
from topoppi_affinity.frozen_probes import load_features
from topoppi_affinity.thermodynamics import REFERENCE_TEMPERATURE_K, pkd_to_dg, pkd_to_kd


def fit_frozen(training, features, output, *, target, representation, family,
               parameters, antisymmetric=False):
    """Fit a fixed configuration on an explicit table and feature NPZ.

    The training table supplies sample_id, complex_id, family_group and target.
    Feature rows are aligned by sample_id; preprocessing is fitted on these rows.
    """
    frame = training.copy() if isinstance(training, pd.DataFrame) else pd.read_csv(
        training, keep_default_na=False, dtype={'sample_id': str})
    if target not in ('pkd', 'ddg'):
        raise ValueError('Target must be pkd or ddg')
    if frame.sample_id.duplicated().any():
        raise ValueError('Training sample IDs must be unique')
    matrix, names = load_features(features, frame)
    with threadpool_limits(limits=2):
        model = fit_frozen_regressor(matrix, frame[target], frame.complex_id,
            family=family, params=parameters, antisymmetric=antisymmetric)
    report = {'target': target, 'representation': representation, 'family': family,
              'parameters': parameters, 'antisymmetric': antisymmetric,
              'n_train': len(frame), 'n_complexes': frame.complex_id.nunique(),
              'n_families': frame.family_group.nunique(), 'n_features': matrix.shape[1],
              'train_sample_ids': frame.sample_id.tolist(),
              'train_family_groups': sorted(frame.family_group.unique()),
              'selection': 'Caller-specified fixed configuration',
              'affinity_units': 'pKd=-log10(Kd/M)' if target == 'pkd' else 'kcal/mol; mutant minus WT'}
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump({'model': model, 'feature_names': names, 'report': report}, output)
    output.with_suffix('.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def predict_frozen(artifact, features, *, temperature_K=REFERENCE_TEMPERATURE_K):
    """Feature NPZ only needs sample_ids, X and feature_names, with no target."""
    bundle = joblib.load(artifact)
    with np.load(features, allow_pickle=False) as source:
        if not np.array_equal(source['feature_names'], bundle['feature_names']):
            raise ValueError('Input feature columns differ from the fitted model')
        with threadpool_limits(limits=2):
            prediction = bundle['model'].predict(source['X'])
        frame = pd.DataFrame({'sample_id': source['sample_ids'].astype(str), 'predicted': prediction})
    target = bundle['report']['target']
    if target == 'pkd':
        frame['predicted_pkd'] = prediction
        frame['predicted_kd_M'] = pkd_to_kd(prediction)
        frame['predicted_dg_at_reference_kcal_mol'] = pkd_to_dg(prediction, temperature_K)
        frame['reference_temperature_K'] = temperature_K
    else:
        frame['predicted_ddg_kcal_mol'] = prediction
    return frame
