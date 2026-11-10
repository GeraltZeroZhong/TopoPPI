"""Coverage-aware reporting for the fixed affinity benchmark; no fitting.

Default reporting requires every declared fold, sample and method. Partial
reports identify their incomplete coverage separately from full out-of-fold predictions.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .benchmark_data import TASKS
from .evaluation import regression_metrics

EXTERNAL = {
    'prodigy': ('absolute', 'PRODIGY', 'core'),
    'ba': ('mutation', 'ba_raw', 'core'),
    'stab': ('mutation', 'stab_raw', 'core'),
    'mirage': ('absolute', 'MIRAGE', 'exposure_reference'),
    'rde': ('mutation', 'RDE-Net', 'exposure_reference'),
    'stab_binding': ('mutation', 'StaB-ddG-binding', 'exposure_reference'),
}
PRIMARY = {
    'absolute_interaction': 'topoppi_surface', 'absolute_family': 'topoppi_surface',
    'mutation_variant': 'known_surface_ensemble',
    'mutation_interaction': 'generalization_surface_ensemble',
    'mutation_family': 'generalization_surface_ensemble',
}


def _bool(values):
    return values.astype(str).str.lower().eq('true')


def _read(path):
    return pd.read_csv(path, keep_default_na=False)


def _validated_predictions(predictions, ids, methods, *, context, allow_partial):
    required = {'sample_id', 'method', 'predicted', 'status'}
    if not required.issubset(predictions):
        raise ValueError(f'{context}: missing prediction columns {sorted(required - set(predictions))}')
    if predictions.duplicated(['sample_id', 'method']).any():
        raise ValueError(f'{context}: each sample/method must occur exactly once, including failure rows.')
    if set(predictions.sample_id) - set(ids) or set(predictions.method) - set(methods):
        raise ValueError(f'{context}: prediction contains undeclared samples or methods.')
    grid = pd.MultiIndex.from_product([list(ids), list(methods)], names=['sample_id', 'method'])
    selected = predictions.set_index(['sample_id', 'method'])
    missing = grid.difference(selected.index)
    if len(missing) and not allow_partial:
        raise ValueError(f'{context}: {len(missing)} declared sample/method rows are missing.')
    selected = selected.reindex(grid).reset_index()
    selected['status'] = selected.status.fillna('not_run')
    for column in ['error', 'components_used', 'training_overlap']:
        if column not in selected:
            selected[column] = ''
        selected[column] = selected[column].fillna('')
    selected['predicted'] = pd.to_numeric(selected.predicted, errors='coerce')
    nonfinite = selected.status.eq('ok') & ~np.isfinite(selected.predicted)
    selected.loc[nonfinite, 'status'] = 'nonfinite_prediction'
    selected.loc[nonfinite, 'error'] = 'Author/model marked ok but prediction is nonfinite'
    return selected, len(missing)


def collect_predictions(root: Path, *, allow_partial=False):
    """Validate declared inventory and attach labels only in this scoring step."""
    root = Path(root)
    plan = json.loads((root / 'execution_plan.json').read_text())
    methods_by_task = plan['methods_by_task']
    expected_tasks = tuple(plan['tasks'])
    if set(expected_tasks) != set(TASKS):
        raise ValueError('Execution plan must declare the five supported benchmark tasks.')
    manifests = {kind: _read(root / 'manifests' / f'{kind}.csv') for kind in ['absolute', 'mutation']}
    for kind, frame in manifests.items():
        label = 'pkd' if kind == 'absolute' else 'ddg'
        if frame.sample_id.duplicated().any() or not np.isfinite(pd.to_numeric(frame[label], errors='coerce')).all():
            raise ValueError(f'{kind}: duplicate sample ID or invalid original target.')
    external, issues = {}, []
    requested_external = list(plan.get('core_external', ['prodigy', 'ba', 'stab']))
    requested_external += [name for name in plan.get('exposed_references', [])
                           if (root / 'priors' / f'{name}.csv').exists()]
    for key in requested_external:
        kind, method, scope = EXTERNAL[key]
        path = root / 'priors' / f'{key}.csv'
        if path.exists():
            raw = _read(path)
            if 'target' in raw and not raw.target.eq('pkd' if kind == 'absolute' else 'ddg').all():
                raise ValueError(f'{path}: target/unit mismatch.')
            raw['method'] = method
        elif not allow_partial:
            raise ValueError(f'Missing required full-pool external predictions: {path}')
        else:
            issues.append(f'Missing external prediction file: {path}')
            raw = pd.DataFrame(columns=['sample_id', 'method', 'predicted', 'status'])
        raw, n_missing = _validated_predictions(raw, manifests[kind].sample_id, [method],
                                                context=str(path), allow_partial=allow_partial)
        if n_missing:
            issues.append(f'{key}: {n_missing} declared external rows missing')
        raw['comparison_scope'] = scope
        external[key] = raw
    all_predictions, inventory = [], []
    for task in expected_tasks:
        kind, _ = TASKS[task]
        frame = manifests[kind]
        split = _read(root / 'splits' / f'{task}.csv')
        if split.sample_id.duplicated().any() or set(split.sample_id) != set(frame.sample_id):
            raise ValueError(f'{task}: split must contain exactly the declared manifest sample IDs.')
        if set(split.fold) != set(range(plan['n_folds'])):
            raise ValueError(f'{task}: split fold IDs disagree with execution plan.')
        metadata = frame[['sample_id', 'source', 'complex_id', 'pdb_id', 'n_mutations', 'family_group',
                          'interaction_group', 'variant_group', 'pkd' if kind == 'absolute' else 'ddg']].copy()
        metadata = metadata.rename(columns={'pkd' if kind == 'absolute' else 'ddg': 'observed'})
        extra = ['fold', 'complex_seen_in_training', 'evaluation_eligible', 'fresh_only',
                 'fresh_row', 'fresh_label', 'fresh_variant', 'fresh_interaction', 'fresh_family',
                 'exposure_roles', 'exposure_manifests']
        metadata = metadata.merge(split[['sample_id'] + extra], on='sample_id', validate='one_to_one')
        metadata['complex_key'] = metadata.source + ':' + metadata.complex_id
        for column in ['complex_seen_in_training', 'evaluation_eligible', 'fresh_only', 'fresh_row',
                       'fresh_label', 'fresh_variant', 'fresh_interaction', 'fresh_family']:
            metadata[column] = _bool(metadata[column])
        metadata['prodigy_nonexact'] = False
        if kind == 'absolute' and 'prodigy' in external:
            prodigy = external['prodigy'].set_index('sample_id')
            if 'exact_training_pdb_overlap' in prodigy:
                known_flag = prodigy.exact_training_pdb_overlap.astype(str).str.lower().isin(['true', 'false'])
                nonexact = known_flag & ~_bool(prodigy.exact_training_pdb_overlap)
            else:
                nonexact = prodigy.training_overlap.isin(['possible', 'excluded'])
            metadata['prodigy_nonexact'] = metadata.sample_id.map(nonexact).fillna(False)
        parts = []
        for fold in range(plan['n_folds']):
            ids = split.loc[split.fold.eq(fold), 'sample_id']
            folder = root / 'runs' / task / f'fold_{fold}'
            path = folder / 'predictions.csv'
            complete = (folder / 'complete.json').exists()
            if not complete and not allow_partial:
                raise ValueError(f'Fold not marked complete: {task}/{fold}')
            if not complete:
                issues.append(f'Incomplete fold: {task}/{fold}')
            if path.exists():
                raw = _read(path)
                if 'task' in raw and not raw.task.eq(task).all():
                    raise ValueError(f'{path}: task mismatch.')
                if 'fold' in raw and not pd.to_numeric(raw.fold).eq(fold).all():
                    raise ValueError(f'{path}: fold mismatch.')
                if 'target' in raw and not raw.target.eq('pkd' if kind == 'absolute' else 'ddg').all():
                    raise ValueError(f'{path}: target mismatch.')
            else:
                raw = pd.DataFrame(columns=['sample_id', 'method', 'predicted', 'status'])
            raw, n_missing = _validated_predictions(raw, ids, methods_by_task[task],
                                                    context=str(path), allow_partial=allow_partial)
            if n_missing:
                issues.append(f'{task}/{fold}: {n_missing} declared prediction rows missing')
            inventory.append({'task': task, 'fold': fold, 'complete_marker': complete,
                              'declared_samples': len(ids), 'declared_methods': len(methods_by_task[task]),
                              'missing_prediction_rows': n_missing})
            parts.append(raw[['sample_id', 'method', 'predicted', 'status', 'error', 'components_used',
                              'training_overlap']])
        task_predictions = pd.concat(parts, ignore_index=True)
        task_predictions['comparison_scope'] = 'core'
        for key, raw in external.items():
            if EXTERNAL[key][0] != kind:
                continue
            method = EXTERNAL[key][1]
            if method in methods_by_task[task]:
                # BA/StaB raw rows already emitted by the fold runner. Verify
                # their usable values rather than double-count the same prior.
                emitted = task_predictions[task_predictions.method.eq(method)].set_index('sample_id')
                prior = raw.set_index('sample_id').reindex(emitted.index)
                if not allow_partial and not emitted.status.eq('ok').equals(prior.status.eq('ok')):
                    raise ValueError(f'{task}: emitted {method} usable coverage differs from its frozen prior.')
                both = emitted.status.eq('ok') & prior.status.eq('ok')
                if not np.allclose(emitted.loc[both, 'predicted'], prior.loc[both, 'predicted'], rtol=1e-10, atol=1e-10):
                    raise ValueError(f'{task}: emitted {method} differs from its frozen prior.')
                mask = task_predictions.method.eq(method)
                task_predictions.loc[mask, 'training_overlap'] = task_predictions.loc[mask, 'sample_id'].map(
                    prior.training_overlap)
            else:
                task_predictions = pd.concat([task_predictions, raw[[
                    'sample_id', 'method', 'predicted', 'status', 'error', 'components_used',
                    'training_overlap', 'comparison_scope']]], ignore_index=True)
        merged = task_predictions.merge(metadata, on='sample_id', validate='many_to_one')
        merged['task'] = task
        merged['target'] = 'pkd' if kind == 'absolute' else 'ddg'
        merged['usable'] = merged.status.eq('ok') & np.isfinite(merged.predicted)
        all_predictions.append(merged)
    return pd.concat(all_predictions, ignore_index=True), pd.DataFrame(inventory), issues, plan


def strata(frame, task):
    """Freshness layers use the task-eligible cohort; known strata are explicit."""
    all_rows = pd.Series(True, index=frame.index)
    masks = {'all': all_rows, 'evaluation_eligible': frame.evaluation_eligible}
    if task.startswith('mutation'):
        masks.update(single=frame.n_mutations.eq(1), multiple=frame.n_mutations.gt(1))
        if task == 'mutation_variant':
            known = frame.complex_seen_in_training
            masks.update(known_complex=known, known_single=known & frame.n_mutations.eq(1),
                         known_multiple=known & frame.n_mutations.gt(1), unseen_complex=~known)
    for name in ['fresh_row', 'fresh_variant', 'fresh_interaction', 'fresh_family', 'fresh_only']:
        if name != 'fresh_variant' or task.startswith('mutation'):
            masks[name] = frame.evaluation_eligible & frame[name]
    if task.startswith('absolute'):
        masks['prodigy_nonexact'] = frame.prodigy_nonexact
    return masks


def metrics(frame):
    if not len(frame):
        return {'n': 0, 'n_groups': 0, 'n_families': 0, 'mae': np.nan, 'rmse': np.nan,
                'macro_mae': np.nan, 'family_balanced_mae': np.nan}
    result = regression_metrics(frame.observed.to_numpy(float), frame.predicted.to_numpy(float), frame.complex_key)
    errors = frame[['complex_key', 'family_group']].copy()
    errors['ae'] = np.abs(frame.predicted.to_numpy(float) - frame.observed.to_numpy(float))
    per_complex = errors.groupby(['family_group', 'complex_key']).ae.mean()
    result['family_balanced_mae'] = float(per_complex.groupby('family_group').mean().mean())
    result['n_families'] = int(frame.family_group.nunique())
    return result


def _bootstrap_ratio(sums, counts, n_bootstrap, seed):
    sums, counts = np.asarray(sums, dtype=float), np.asarray(counts, dtype=float)
    if len(sums) < 2 or not n_bootstrap:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    draws = []
    for start in range(0, n_bootstrap, 100):
        selected = rng.integers(len(sums), size=(min(100, n_bootstrap - start), len(sums)))
        draws.extend((sums[selected].sum(axis=1) / counts[selected].sum(axis=1)).tolist())
    return tuple(np.quantile(draws, [.025, .975]))


def paired_comparison(left, right, *, task, n_bootstrap=2000, seed=42):
    """Compare on the exact intersection, retaining equal-complex/family views."""
    common = left[left.usable].merge(right.loc[right.usable, ['sample_id', 'predicted']],
                                    on='sample_id', suffixes=('_left', '_right'), validate='one_to_one')
    if not len(common):
        return {'n_matched': 0, 'n_complexes': 0, 'n_families': 0}
    common['ae_left'] = np.abs(common.predicted_left - common.observed)
    common['ae_right'] = np.abs(common.predicted_right - common.observed)
    common['delta'] = common.ae_left - common.ae_right
    unit = 'family_group' if task.endswith('_family') else (
        'interaction_group' if task.endswith('_interaction') else 'complex_key')
    rows = common.groupby(unit).delta.agg(['sum', 'count'])
    complex_columns = list(dict.fromkeys([unit, 'complex_key']))
    complexes = common.groupby(complex_columns, as_index=False).delta.mean()
    macro = complexes.groupby(unit).delta.agg(['sum', 'count'])
    families = common.groupby(['family_group', 'complex_key']).delta.mean().groupby('family_group').mean()
    result = {'n_matched': len(common), 'n_complexes': common.complex_key.nunique(),
              'n_families': common.family_group.nunique(), 'bootstrap_unit': unit,
              'bootstrap_n_groups': len(rows), 'family_bootstrap_n_groups': len(families),
              'mae_difference': float(common.delta.mean()),
              'macro_mae_difference': float(complexes.delta.mean()),
              'family_balanced_mae_difference': float(families.mean()),
              'left_mae': float(np.abs(common.predicted_left - common.observed).mean()),
              'right_mae': float(np.abs(common.predicted_right - common.observed).mean())}
    for side in ['left', 'right']:
        per_complex = common.groupby(['family_group', 'complex_key'])[f'ae_{side}'].mean()
        result[f'{side}_macro_mae'] = float(per_complex.mean())
        result[f'{side}_family_balanced_mae'] = float(per_complex.groupby('family_group').mean().mean())
    for name, sums, counts in [('mae', rows['sum'], rows['count']),
                               ('macro_mae', macro['sum'], macro['count']),
                               ('family_balanced_mae', families, np.ones(len(families)))]:
        low, high = _bootstrap_ratio(sums, counts, n_bootstrap, seed)
        result[f'{name}_difference_ci_low'], result[f'{name}_difference_ci_high'] = low, high
    return result



def primary_tables(all_metrics, all_pairs, plan):
    """Keep prespecified primary strata and orient every core pair consistently."""
    primary = plan.get('primary_methods', PRIMARY)
    strata_by_task = {
        task: (['known_complex', 'known_single', 'known_multiple'] if task == 'mutation_variant' else
               ['all', 'single', 'multiple'] if task.startswith('mutation') else ['all'])
        for task in TASKS
    }
    metrics_mask = pd.Series(False, index=all_metrics.index)
    rows = []
    for task, method in primary.items():
        selected = (all_metrics.task.eq(task) & all_metrics.method.eq(method)
                    & all_metrics.stratum.isin(strata_by_task[task]) & all_metrics.comparison_scope.eq('core'))
        metrics_mask |= selected
        pairs = all_pairs[all_pairs.task.eq(task) & all_pairs.stratum.isin(strata_by_task[task])
                          & all_pairs.comparison_scope.eq('core')
                          & (all_pairs.left_method.eq(method) | all_pairs.right_method.eq(method))]
        for _, original in pairs.iterrows():
            row = original.to_dict()
            if row['right_method'] == method:
                row['left_method'], row['right_method'] = row['right_method'], row['left_method']
                for metric in ['mae', 'macro_mae', 'family_balanced_mae']:
                    left, right = f'left_{metric}', f'right_{metric}'
                    row[left], row[right] = original.get(right, np.nan), original.get(left, np.nan)
                    difference = f'{metric}_difference'
                    row[difference] = -original.get(difference, np.nan)
                    low, high = difference + '_ci_low', difference + '_ci_high'
                    row[low], row[high] = -original.get(high, np.nan), -original.get(low, np.nan)
            row['primary_method'], row['comparator'] = row['left_method'], row['right_method']
            rows.append(row)
    return all_metrics[metrics_mask].copy(), pd.DataFrame(rows)

def summarize_benchmark(root: Path, *, output: Path | None = None, allow_partial=False,
                         n_bootstrap=2000, seed=42):
    root = Path(root)
    output = Path(output) if output is not None else root / 'report'
    matched, inventory, issues, plan = collect_predictions(root, allow_partial=allow_partial)
    metric_rows, pair_rows, coverage = [], [], []
    for (task, fold, method), part in matched.groupby(['task', 'fold', 'method'], sort=True):
        coverage.append({'task': task, 'fold': fold, 'method': method, 'comparison_scope': part.comparison_scope.iloc[0],
                         'n_declared': len(part), 'n_attempted': int(part.status.ne('not_run').sum()),
                         'n_usable': int(part.usable.sum()), 'n_failed': int((~part.usable & part.status.ne('not_run')).sum()),
                         'n_task_eligible': int(part.evaluation_eligible.sum()),
                         'n_eligible_usable': int((part.usable & part.evaluation_eligible).sum()),
                         'n_base_only': int((part.usable & part.components_used.eq('base_only')).sum()),
                         'status_counts': json.dumps(part.status.value_counts().to_dict(), sort_keys=True)})
    for task, task_frame in matched.groupby('task', sort=True):
        for stratum, mask in strata(task_frame, task).items():
            subset = task_frame[mask]
            methods = sorted(task_frame.method.unique())
            by_method = {method: subset[subset.method.eq(method)] for method in methods}
            for method, part in by_method.items():
                scope = task_frame.loc[task_frame.method.eq(method), 'comparison_scope'].iloc[0]
                metric_rows.append({'task': task, 'stratum': stratum, 'method': method,
                                    'comparison_scope': scope, 'n_declared': len(part),
                                    'n_usable': int(part.usable.sum()), **metrics(part[part.usable])})
            core = [method for method in methods if task_frame.loc[task_frame.method.eq(method),
                                                                   'comparison_scope'].iloc[0] == 'core']
            pairs = [(left, right, 'core') for left, right in itertools.combinations(core, 2)]
            primary = plan.get('primary_methods', PRIMARY)[task]
            pairs += [(method, primary, 'exposure_reference') for method in methods if method not in core]
            for left, right, scope in pairs:
                result = paired_comparison(by_method[left], by_method[right], task=task,
                                           n_bootstrap=n_bootstrap, seed=seed)
                pair_rows.append({'task': task, 'stratum': stratum, 'left_method': left, 'right_method': right,
                                  'comparison_scope': scope, **result})
    output.mkdir(parents=True, exist_ok=True)
    all_metrics, all_pairs = pd.DataFrame(metric_rows), pd.DataFrame(pair_rows)
    primary_metrics, primary_pairs = primary_tables(all_metrics, all_pairs, plan)
    primary_metrics.to_csv(output / 'primary_results.csv', index=False)
    primary_pairs.to_csv(output / 'primary_comparisons.csv', index=False)
    for scope, prefix in [('core', ''), ('exposure_reference', 'reference_')]:
        all_metrics[all_metrics.comparison_scope.eq(scope)].to_csv(output / f'{prefix}metrics.csv', index=False)
        all_pairs[all_pairs.comparison_scope.eq(scope)].to_csv(output / f'{prefix}paired_comparisons.csv', index=False)
    pd.DataFrame(coverage).to_csv(output / 'coverage.csv', index=False)
    inventory.to_csv(output / 'fold_inventory.csv', index=False)
    matched.to_csv(output / 'matched_predictions.csv', index=False)
    matched[~matched.usable].to_csv(output / 'failed_predictions.csv', index=False)
    metadata = json.loads((root / 'protocol.json').read_text())
    report = {
        'scope': 'Partial benchmark out-of-fold coverage' if allow_partial else 'Completed expanded fixed-recipe benchmark OOF',
        'allow_partial': allow_partial, 'completeness_issues': issues,
        'completed_folds': int(inventory.complete_marker.sum()), 'declared_folds': len(inventory),
        'n_bootstrap': n_bootstrap, 'seed': seed, 'paired_difference': 'left error minus right error; negative favors left',
        'bootstrap': 'Task grouping for row/complex-macro MAE; family resampling for family-balanced MAE. Fewer than two clusters: CI unavailable.',
        'family_balanced': 'Average absolute error within complex, then average complexes within family, then equal family average.',
        'primary_comparison_direction': 'primary minus comparator; negative favors primary; core comparators only',
        'primary_matched_points': 'Both methods MAE, complex macro MAE and family-balanced MAE use the same paired successful sample IDs.',
        'primary_stratum': {task: 'known_complex' if task == 'mutation_variant' else 'all' for task in TASKS},
        'freshness': metadata['freshness'], 'fresh_family_interpretation': metadata['fresh_family_interpretation'],
        'prodigy_nonexact': 'The same nonexact-PDB subset is used for every absolute method; homologous and pretraining exposure can remain.',
        'exposed_references': 'MIRAGE, binding-supervised RDE-Net and final StaB-ddG are reported in exposure-reference tables.',
        'coverage': 'Failed rows retain their status; comparisons use exact common successful sample IDs and preserve missing predictions.',
        'selection_limit': 'Development, auxiliary and previously reserved labels informed the model recipes. Family independence is assessed separately from row novelty; comparisons apply to this benchmark and its reported exposure conditions.',
        'units': {'absolute': 'pKd; Kd and reference-temperature deltaG are derived outputs', 'mutation': 'kcal/mol, positive weaker binding'},
        'training_performed_by_report': False,
    }
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    text = [f'# {report["scope"]}', '', report['fresh_family_interpretation'], '',
            'All results use the declared held-out fold. Failure coverage is retained and paired differences use identical successful IDs.',
            'Mutation variant primary evaluation is the strict known-complex layer; all-row and unseen-complex strata are separate.',
            'The same PRODIGY nonexact-PDB strata are applied to every absolute method; homologous exposure can remain.',
            'Binding-supervised and potentially overlapping reference models are reported separately with their exposure conditions.',
            'Development, auxiliary and previously reserved labels informed recipe selection. Bootstrap intervals describe resampling uncertainty conditional on those recipes.', '',
            f'Fold completion: {report["completed_folds"]}/{report["declared_folds"]}.',
            'Prespecified primary strata: primary_results.csv; primary-minus-comparator differences: primary_comparisons.csv.',
            'All-method metrics: metrics.csv; all pairs: paired_comparisons.csv; reference_* files contain exposed methods.',
            'Coverage/status: coverage.csv, failed_predictions.csv; row alignment: matched_predictions.csv.', '']
    (output / 'README.md').write_text('\n'.join(text))
    return report
