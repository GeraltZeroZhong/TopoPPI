"""Refit and unlabeled inference for frozen-prior plus Topo kernel residuals."""
from __future__ import annotations

import json
import time
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from topoppi_affinity.embedding_baselines import zero_surface_slots
from topoppi_affinity.residual_baselines import align_prior, fit_kernel_residual

PRIOR_COLUMNS = ['sample_id', 'predicted', 'target', 'method', 'status', 'error']


def _read_features(path, sample_ids):
    with np.load(path, allow_pickle=False) as source:
        ids = source['sample_ids'].astype(str)
        names, matrix = source['feature_names'].astype(str), source['X']
    if len(set(ids)) != len(ids):
        raise ValueError('Duplicate sample IDs in feature NPZ')
    if matrix.shape != (len(ids), len(names)):
        raise ValueError('Feature matrix shape does not match sample IDs and feature names')
    lookup = {sample: index for index, sample in enumerate(ids)}
    missing = set(sample_ids) - set(lookup)
    if missing:
        raise ValueError(f'Feature NPZ lacks requested sample IDs: {sorted(missing)[:5]}')
    return matrix[[lookup[str(sample)] for sample in sample_ids]], names, ids


def predict_residual(artifact, features, prior):
    """Return unlabeled predictions; artifact may be a loaded bundle or joblib path.

    Features is an NPZ path; prior is a CSV path or a DataFrame with PRIOR_COLUMNS.
    Output preserves prior row order and explicit failures.
    """
    bundle = artifact if isinstance(artifact, dict) else joblib.load(artifact)
    metadata, model = bundle['report'], bundle['model']
    prior = prior[PRIOR_COLUMNS].copy() if isinstance(prior, pd.DataFrame) else pd.read_csv(
        prior, usecols=PRIOR_COLUMNS, keep_default_na=False, dtype={'sample_id': str})
    prior['error'] = prior.error.fillna('')
    if prior.sample_id.duplicated().any():
        raise ValueError('Duplicate sample IDs in prior CSV')
    if not prior.target.eq('ddg').all():
        raise ValueError('Prior target must be ddg in kcal/mol')
    if not prior.method.eq(metadata['prior_method']).all():
        raise ValueError(f"Artifact requires prior method {metadata['prior_method']}")
    matrix, names, ids = _read_features(features, prior.sample_id)
    if set(ids) != set(prior.sample_id):
        raise ValueError('Prior CSV and feature NPZ must contain identical sample-ID sets')
    if not np.array_equal(names, bundle['feature_names']):
        raise ValueError('Feature columns differ from the fitted residual artifact')
    if metadata['representation'] == 'contact':
        matrix = zero_surface_slots(matrix, names)
    values = pd.to_numeric(prior.predicted, errors='coerce').to_numpy(dtype=float)
    prior_ok = prior.status.eq('ok').to_numpy() & np.isfinite(values)
    feature_ok = np.isfinite(matrix).all(axis=1)
    usable = prior_ok & feature_ok
    result = prior[['sample_id']].copy()
    result['target'], result['method'] = 'ddg', metadata['method']
    result['prior'] = values
    result['scaled_prior'], result['residual'], result['predicted'] = np.nan, np.nan, np.nan
    result['status'] = np.where(usable, 'ok', 'error')
    result['error'] = prior.error
    result.loc[~prior_ok & result.error.eq(''), 'error'] = 'Failed or nonfinite prior prediction'
    result.loc[prior_ok & ~feature_ok, 'error'] = 'Nonfinite residual features'
    if usable.any():
        with threadpool_limits(limits=2):
            predicted = model.predict(values[usable], matrix[usable])
        scaled = model.slope * values[usable]
        result.loc[usable, ['scaled_prior', 'residual', 'predicted']] = np.column_stack((scaled, predicted-scaled, predicted))
        result.loc[usable, 'error'] = ''
    return result


def predict_residual_affinity(artifact, prior_path, features_path, output):
    """CLI-facing inference writer; neither input needs observed affinity labels."""
    artifact, output = Path(artifact), Path(output)
    bundle = joblib.load(artifact)
    metadata = bundle['report']
    result = predict_residual(bundle, features_path, prior_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)
    output.with_suffix('.json').write_text(json.dumps({
        'artifact': str(artifact.resolve()), 'method': metadata['method'], 'target': 'ddg',
        'units': metadata['units'], 'prior_method': metadata['prior_method'], 'n_attempted': len(result),
        'n_success': int(result.status.eq('ok').sum()), 'observed_labels_loaded': False,
        'components': 'predicted = scaled_prior + residual; scaled_prior = fitted_origin_slope * prior',
    }, indent=2) + '\n')
    return result


def fit_residual_model(training, features, prior_path, output, *, task='mutation_within', representation='surface', seed=11):
    """Select RBF parameters within supplied training groups, then fit successful rows."""
    if task not in ('mutation_within', 'mutation_family') or representation not in ('surface', 'contact'):
        raise ValueError('Supported refits are mutation_within/family and surface/contact')
    prior_path, output = Path(prior_path), Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    full_frame = training.copy() if isinstance(training, pd.DataFrame) else pd.read_csv(
        training, keep_default_na=False, dtype={'sample_id': str})
    predictions = pd.read_csv(prior_path, usecols=PRIOR_COLUMNS, keep_default_na=False, dtype={'sample_id': str})
    frame, prior, _, coverage = align_prior(full_frame, predictions)
    methods = predictions.loc[predictions.sample_id.isin(frame.sample_id), 'method'].unique().tolist()
    if len(methods) != 1:
        raise ValueError('A residual artifact must use exactly one prior method')
    matrix, names, _ = _read_features(features, frame.sample_id)
    if representation == 'contact':
        matrix = zero_surface_slots(matrix, names)
    if not np.isfinite(matrix).all():
        raise ValueError('Refit features contain nonfinite entries')
    start = time.monotonic()
    with threadpool_limits(limits=2):
        model, params, candidates = fit_kernel_residual(matrix, prior, frame, task=task, seed=seed)
    metadata = {
        'method': f'{methods[0]}-Topo-{representation}-residual', 'target': 'ddg',
        'units': 'kcal/mol; positive weakens binding', 'task': task, 'representation': representation,
        'parameters': params, 'origin_slope': model.slope, 'seed': seed,
        'n_train': len(frame), 'n_failed_prior': len(full_frame)-len(frame), 'n_features': len(names),
        'train_sample_ids': frame.sample_id.tolist(), 'train_complex_ids': frame.complex_id.tolist(),
        'train_family_groups': frame.family_group.tolist(), 'prior_method': methods[0],
        'prior_source': str(prior_path.resolve()),
        'prior_provenance': json.loads(prior_path.with_suffix('.json').read_text())
            if prior_path.with_suffix('.json').exists() else {'status': 'metadata_not_provided'},
        'selection': 'Two inner training folds grouped by variant or family according to task; minimum macro MAE',
        'weighting': 'Equal total weight per complex; odd residual reverse augmentation with half weights',
        'contact_policy': 'All surface.* feature slots zeroed in both refit and inference for contact artifacts',
        'input_contract': 'Unlabeled prior CSV in kcal/mol and NPZ with sample_ids/X/feature_names; exact sample-ID alignment',
        'symmetry': 'Odd residual and no-intercept calibration; raw prior symmetry is not newly established',
        'elapsed_seconds': time.monotonic()-start,
    }
    joblib.dump({'model': model, 'feature_names': names, 'report': metadata}, output)
    output.with_suffix('.json').write_text(json.dumps(metadata, indent=2) + '\n')
    coverage.to_csv(output.with_suffix('.coverage.csv'), index=False)
    pd.DataFrame(candidates).to_csv(output.with_suffix('.inner_candidates.csv'), index=False)
    return metadata
