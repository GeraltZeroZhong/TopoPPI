"""Fixed equal-weight fusion of fitted mutation predictors.

The components retain their own preprocessing and training metadata. No stacking
weights or confirmation labels are fitted by this module.
"""
from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits


def bundle_mutation_ensemble(base_artifact, residual_artifact, output):
    base, residual = joblib.load(base_artifact), joblib.load(residual_artifact)
    if base['report']['target'] != 'ddg' or residual['report']['target'] != 'ddg':
        raise ValueError('Both ensemble components must predict mutation ddg')
    report = {
        'target': 'ddg', 'units': 'kcal/mol; mutant minus WT; positive weakens binding',
        'combination': '0.5 * base + 0.5 * prior-and-kernel-residual prediction',
        'missing_prior': 'Use the base prediction alone and report components_used=base_only',
        'selection': 'Fixed equal weights; no stacking fit or confirmation labels',
        'base': base['report'], 'residual': residual['report'],
        'symmetry': 'No additional odd or thermodynamic-cycle guarantee from averaging',
        'source_artifacts': [str(base_artifact), str(residual_artifact)],
    }
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump({'base': base, 'residual': residual, 'report': report}, output)
    output.with_suffix('.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def predict_mutation_ensemble(artifact, base_features, residual_features, prior):
    """Predict unlabeled inputs; preserve the prior input's sample order."""
    from topoppi_affinity.residual_artifacts import predict_residual

    bundle = joblib.load(artifact) if isinstance(artifact, (str, Path)) else artifact
    residual = predict_residual(bundle['residual'], residual_features, prior)
    with np.load(base_features, allow_pickle=False) as features:
        if not np.array_equal(features['feature_names'], bundle['base']['feature_names']):
            raise ValueError('Base input feature columns differ from its fitted model')
        ids = features['sample_ids'].astype(str)
        if len(np.unique(ids)) != len(ids) or set(ids) != set(residual.sample_id):
            raise ValueError('Base and residual inputs require the same unique sample IDs')
        with threadpool_limits(limits=2):
            prediction = bundle['base']['model'].predict(features['X'])
    result = residual.copy()
    result['method'] = 'fixed_equal_mutation_ensemble'
    result['base_prediction'] = result.sample_id.map(pd.Series(prediction, index=ids))
    result['residual_prediction'] = result['predicted']
    valid = result.status.eq('ok') & np.isfinite(result.predicted)
    result['predicted'] = result.base_prediction
    result.loc[valid, 'predicted'] = .5 * (
        result.loc[valid, 'base_prediction'] + result.loc[valid, 'residual_prediction'])
    result['components_used'] = np.where(valid, 'equal_ensemble', 'base_only')
    result['predicted_ddg_kcal_mol'] = result.predicted
    result['prior_status'] = result.status
    result['status'] = np.where(valid, 'ok', 'base_only')
    return result


def bundle_scalar_prior_ensemble(base_artifact, prior_method, output):
    """Package the selected family-generalization recipe with a frozen BA score."""
    base = joblib.load(base_artifact)
    if base['report']['target'] != 'ddg':
        raise ValueError('The base model must predict mutation ddg')
    report = {
        'target': 'ddg', 'units': 'kcal/mol; mutant minus WT; positive weakens binding',
        'combination': '0.5 * base + 0.5 * frozen RT-scaled BA-Cycle score',
        'prior_method': prior_method, 'base': base['report'],
        'selection': 'Fixed equal weights selected on development; no confirmation labels',
        'missing_prior': 'Use base alone and explicitly report base_only',
        'symmetry': 'No exact antisymmetry or thermodynamic cycle guarantee for the WT-context prior',
        'source_artifacts': [str(base_artifact)],
    }
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump({'base': base, 'prior_method': prior_method, 'report': report}, output)
    output.with_suffix('.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def predict_scalar_prior_ensemble(artifact, base_features, prior):
    """Predict from one learned base and one supplied frozen scalar-prior CSV."""
    bundle = joblib.load(artifact) if isinstance(artifact, (str, Path)) else artifact
    prior = prior.copy() if isinstance(prior, pd.DataFrame) else pd.read_csv(
        prior, keep_default_na=False, dtype={'sample_id': str})
    if not prior.target.eq('ddg').all() or not prior.method.eq(bundle['prior_method']).all():
        raise ValueError('The prior method and ddg target must match the fitted recipe')
    with np.load(base_features, allow_pickle=False) as features:
        ids = features['sample_ids'].astype(str)
        if prior.sample_id.duplicated().any() or len(np.unique(ids)) != len(ids) or set(ids) != set(prior.sample_id):
            raise ValueError('Base features and scalar prior require the same unique sample IDs')
        if not np.array_equal(features['feature_names'], bundle['base']['feature_names']):
            raise ValueError('Base input feature columns differ from its fitted model')
        with threadpool_limits(limits=2):
            prediction = bundle['base']['model'].predict(features['X'])
    result = prior[['sample_id', 'status', 'error']].copy()
    result['base_prediction'] = result.sample_id.map(pd.Series(prediction, index=ids))
    result['prior_prediction'] = pd.to_numeric(prior.predicted, errors='coerce')
    usable = result.status.eq('ok') & np.isfinite(result.prior_prediction)
    result['predicted'] = result.base_prediction
    result.loc[usable, 'predicted'] = .5*(result.loc[usable, 'base_prediction']+result.loc[usable, 'prior_prediction'])
    result['components_used'] = np.where(usable, 'equal_ensemble', 'base_only')
    result['prior_status'] = result.status
    result['status'] = np.where(usable, 'ok', 'base_only')
    result['target'], result['method'] = 'ddg', 'fixed_base_BA-Cycle_ensemble'
    result['predicted_ddg_kcal_mol'] = result.predicted
    return result
