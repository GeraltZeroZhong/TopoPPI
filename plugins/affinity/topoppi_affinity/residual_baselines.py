"""Scalar-prior calibration and grouped kernel residual fitting."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

from topoppi_affinity.embedding_baselines import (
    candidate_parameters,
    fit_frozen_regressor,
)
from topoppi_affinity.evaluation import regression_metrics
from topoppi_affinity.models import complex_weights

METHODS = ('prior', 'origin_calibrated', 'plm_contact_residual', 'plm_surface_residual')


def align_prior(frame, predictions):
    """Retain successful finite predictions once, preserving manifest row order."""
    if frame.sample_id.duplicated().any() or predictions.sample_id.duplicated().any():
        raise ValueError('Duplicate sample IDs in manifest or prior predictions')
    if not (predictions.target == 'ddg').all():
        raise ValueError('Residual prior must predict ddg in kcal/mol')
    prediction = predictions.set_index('sample_id').reindex(frame.sample_id)
    values = pd.to_numeric(prediction.predicted, errors='coerce').to_numpy(dtype=float)
    usable = prediction.status.eq('ok').to_numpy() & np.isfinite(values)
    coverage = frame[['sample_id', 'complex_id']].copy()
    coverage['included'] = usable
    coverage['prior_status'] = prediction.status.fillna('missing').to_numpy()
    coverage['error'] = prediction.error.fillna('').to_numpy()
    coverage.loc[(coverage.prior_status == 'ok') & ~usable, 'error'] = 'Nonfinite prior prediction'
    return frame.loc[usable].copy().reset_index(drop=True), values[usable], usable, coverage


def fit_origin_slope(prior, labels, groups):
    """Complex-weighted least squares with zero intercept, using training only."""
    prior, labels = np.asarray(prior, float), np.asarray(labels, float)
    weights = complex_weights(groups)
    denominator = np.sum(weights * prior * prior)
    if denominator <= 0 or not np.isfinite(denominator):
        raise ValueError('An origin slope requires finite, nonzero training priors')
    return float(np.sum(weights * prior * labels) / denominator)


@dataclass
class ResidualRegressor:
    slope: float
    residual_model: object = None

    def predict(self, prior, features=None):
        predicted = self.slope * np.asarray(prior, float)
        if self.residual_model is not None:
            predicted = predicted + self.residual_model.predict(features)
        return predicted


def inner_partitions(frame, task, seed=11):
    """Two inner folds preserve the outer task's mutation or family boundary."""
    if 'within' in task:
        units = frame[['complex_id', 'mutations']].drop_duplicates().reset_index(drop=True)
        units['unit'] = units.complex_id + '|' + units.mutations
        groups = units.unit.to_numpy()
        row_units = frame.complex_id + '|' + frame.mutations
        splitter = GroupKFold(n_splits=2, shuffle=True, random_state=seed)
    else:
        column = 'family_group' if 'family' in task else 'initial_interaction_group'
        units = frame[['complex_id', column]].drop_duplicates('complex_id').reset_index(drop=True)
        groups = units[column].to_numpy()
        units = units.assign(unit=units.complex_id)
        row_units = frame.complex_id
        splitter = GroupKFold(n_splits=2)
    for train, validation in splitter.split(units, groups=groups):
        yield (row_units.isin(set(units.iloc[train].unit)).to_numpy(),
               row_units.isin(set(units.iloc[validation].unit)).to_numpy())


def fit_kernel_residual(matrix, prior, frame, *, task, seed=11):
    """Choose RBF parameters on inner OOF predictions and fit all supplied rows."""
    labels, groups = frame.ddg.to_numpy(dtype=float), frame.complex_id.to_numpy()
    splits = list(inner_partitions(frame, task, seed))
    candidates, best_params, best_score = [], None, float('inf')
    for index, params in enumerate(candidate_parameters('rbf', matrix.shape[1])):
        predicted = np.full(len(frame), np.nan)
        for train, validation in splits:
            slope = fit_origin_slope(prior[train], labels[train], groups[train])
            model = fit_frozen_regressor(
                matrix[train], labels[train] - slope * prior[train], groups[train],
                family='rbf', params=params, antisymmetric=True, seed=seed,
            )
            predicted[validation] = slope * prior[validation] + model.predict(matrix[validation])
        metrics = regression_metrics(labels, predicted, groups)
        candidates.append({'candidate': index, **params, **metrics})
        if metrics['macro_mae'] < best_score:
            best_params, best_score = params, metrics['macro_mae']
    slope = fit_origin_slope(prior, labels, groups)
    model = fit_frozen_regressor(matrix, labels - slope * prior, groups, family='rbf',
                                 params=best_params, antisymmetric=True, seed=seed)
    return ResidualRegressor(slope, model), best_params, candidates


def paired_oof_comparisons(output, *, n_bootstrap=2000, seed=42):
    """Paired surface increments; family or complex clusters retain dependence."""
    output = Path(output)
    rows = []
    for surface_path in sorted((output / 'oof_predictions').glob('*_plm_surface_residual.csv')):
        task = surface_path.name.removesuffix('_plm_surface_residual.csv')
        surface = pd.read_csv(surface_path)
        for comparator in ('origin_calibrated', 'plm_contact_residual'):
            path = output / 'oof_predictions' / f'{task}_{comparator}.csv'
            if not path.exists():
                continue
            control = pd.read_csv(path)
            joined = surface.merge(control[['sample_id', 'predicted']], on='sample_id',
                                   suffixes=('_surface', '_control'), validate='one_to_one')
            if not len(joined) == len(surface) == len(control):
                raise ValueError('Paired residual comparison requires identical OOF sample IDs')
            joined['difference'] = abs(joined.predicted_surface - joined.observed) - abs(joined.predicted_control - joined.observed)
            cluster = 'family_group' if task.endswith('family') else 'complex_id'
            micro = joined.groupby(cluster).difference.agg(['sum', 'count'])
            per_complex = joined.groupby(['family_group', 'complex_id']).difference.mean().reset_index()
            macro = per_complex.groupby(cluster).difference.agg(['sum', 'count']).reindex(micro.index)
            generator = np.random.default_rng(seed)
            selected = generator.integers(len(micro), size=(n_bootstrap, len(micro)))
            draws = [values['sum'].to_numpy()[selected].sum(1) / values['count'].to_numpy()[selected].sum(1)
                     for values in (micro, macro)]
            rows.append({
                'task': task, 'contrast': 'plm_surface_residual-minus-' + comparator, 'n': len(joined),
                'cluster': cluster, 'n_clusters': len(micro), 'mae_difference': float(joined.difference.mean()),
                'mae_difference_ci_low': float(np.quantile(draws[0], .025)),
                'mae_difference_ci_high': float(np.quantile(draws[0], .975)),
                'macro_mae_difference': float(per_complex.difference.mean()),
                'macro_difference_ci_low': float(np.quantile(draws[1], .025)),
                'macro_difference_ci_high': float(np.quantile(draws[1], .975)),
            })
    pd.DataFrame(rows).to_csv(output / 'paired_oof_comparisons.csv', index=False)
    return rows
