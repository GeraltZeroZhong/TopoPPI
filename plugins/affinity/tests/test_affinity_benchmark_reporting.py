from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from topoppi_affinity.benchmark_data import TASKS
from topoppi_affinity.benchmark_reporting import (
    PRIMARY,
    collect_predictions,
    metrics,
    paired_comparison,
    summarize_benchmark,
)


def make_benchmark(root):
    for directory in ['manifests', 'splits', 'priors']:
        (root / directory).mkdir(parents=True)
    frames = {}
    for kind, n in [('absolute', 5), ('mutation', 11)]:
        rows = []
        for index in range(n):
            system = index if kind == 'absolute' else index // 2
            rows.append(dict(sample_id=f'{kind}{index}', source=kind, complex_id=f'c{system}',
                             pdb_id=f'p{system}', n_mutations=0 if kind == 'absolute' else 1 + index % 2,
                             family_group=f'f{system}', interaction_group=f'i{system}', variant_group=f'v{index}',
                             pkd=index * .5 - 1, ddg=index * .5 - 1))
        frames[kind] = pd.DataFrame(rows)
        frames[kind].to_csv(root / 'manifests' / f'{kind}.csv', index=False)
    methods = {}
    priors = {}
    for key, kind, method in [('prodigy', 'absolute', 'PRODIGY'), ('ba', 'mutation', 'BA-Cycle-partner-MC8'),
                              ('stab', 'mutation', 'StaB-ddG-stability-only'), ('rde', 'mutation', 'RDE-Net')]:
        f = frames[kind]
        raw = pd.DataFrame(dict(sample_id=f.sample_id, method=method, predicted=np.arange(len(f)) * .2,
                                target='pkd' if kind == 'absolute' else 'ddg', status='ok', error=''))
        if key == 'stab':
            raw.loc[len(raw) - 1, ['predicted', 'status']] = [np.nan, 'error']
        if key == 'prodigy':
            raw['exact_training_pdb_overlap'] = [True, False, False, False, False]
        raw.to_csv(root / 'priors' / f'{key}.csv', index=False)
        priors[key] = raw.set_index('sample_id')
    for task, (kind, _) in TASKS.items():
        frame = frames[kind]
        declared = [PRIMARY[task], 'contact_control']
        if kind == 'mutation':
            declared.append('ba_raw')
            if task == 'mutation_variant':
                declared.append('stab_raw')
        methods[task] = declared
        split = frame[['sample_id']].copy()
        split['fold'] = (np.arange(len(frame)) if task == 'mutation_variant' else
                         frame.complex_id.str[1:].astype(int).to_numpy()) % 5
        seen = frame.assign(fold=split.fold).groupby('complex_id').fold.transform('nunique').gt(1)
        split['complex_seen_in_training'] = seen if task == 'mutation_variant' else False
        split['evaluation_eligible'] = split.complex_seen_in_training if task == 'mutation_variant' else True
        for name in ['fresh_row', 'fresh_label', 'fresh_variant', 'fresh_interaction']:
            split[name] = np.arange(len(frame)) % 2 == 0
        split['fresh_family'] = kind == 'absolute'
        split['fresh_only'] = split.fresh_row & split.evaluation_eligible
        split['exposure_roles'], split['exposure_manifests'] = '', ''
        split.to_csv(root / 'splits' / f'{task}.csv', index=False)
        for fold in range(5):
            folder = root / 'runs' / task / f'fold_{fold}'
            folder.mkdir(parents=True)
            rows = []
            for index in split.index[split.fold.eq(fold)]:
                sample = frame.iloc[index]
                for method in declared:
                    predicted, status = sample.pkd + (.5 if method == PRIMARY[task] else 1.), 'ok'
                    if method in ('ba_raw', 'stab_raw'):
                        p = priors['ba' if method == 'ba_raw' else 'stab'].loc[sample.sample_id]
                        predicted, status = p.predicted, p.status
                    if method == 'contact_control' and index == len(frame) - 1:
                        predicted, status = np.nan, 'missing_features'
                    rows.append(dict(sample_id=sample.sample_id, method=method, predicted=predicted,
                                     status=status, components_used='base_only' if index == len(frame) - 1 else '',
                                     task=task, fold=fold))
            pd.DataFrame(rows).to_csv(folder / 'predictions.csv', index=False)
            (folder / 'complete.json').write_text('{}')
    (root / 'execution_plan.json').write_text(json.dumps(dict(n_folds=5, tasks=list(TASKS),
        methods_by_task=methods, primary_methods=PRIMARY, core_external=['prodigy', 'ba', 'stab'],
        exposed_references=['rde'])))
    (root / 'protocol.json').write_text(json.dumps(dict(freshness={'fresh_row': 'exact sample ID not previously used'},
        fresh_family_interpretation='Mutation fresh-family count is zero; expanded CV, not fresh confirmation.')))
    return frames


def test_complete_report_coverage_same_ids_and_exposure_references(tmp_path):
    frames = make_benchmark(tmp_path)
    report = summarize_benchmark(tmp_path, n_bootstrap=30)
    assert report['completed_folds'] == 25
    assert not report['completeness_issues']
    output = tmp_path / 'report'
    main = pd.read_csv(output / 'metrics.csv')
    references = pd.read_csv(output / 'reference_metrics.csv')
    assert 'RDE-Net' not in set(main.method)
    assert set(references.method) == {'RDE-Net'}
    pairs = pd.read_csv(output / 'paired_comparisons.csv')
    pair = pairs[(pairs.task == 'absolute_family') & (pairs.stratum == 'all') &
                 (pairs.left_method == 'contact_control') & (pairs.right_method == 'topoppi_surface')].iloc[0]
    assert pair.n_matched == 4
    assert pair.mae_difference == pytest.approx(.5)
    assert pair.left_macro_mae == pytest.approx(1.)
    assert pair.right_macro_mae == pytest.approx(.5)
    assert pair.left_family_balanced_mae == pytest.approx(1.)
    assert pair.right_family_balanced_mae == pytest.approx(.5)
    primary = pd.read_csv(output / 'primary_results.csv')
    assert set(primary.loc[primary.task.eq('mutation_variant'), 'stratum']) == {
        'known_complex', 'known_single', 'known_multiple'}
    primary_pairs = pd.read_csv(output / 'primary_comparisons.csv')
    oriented = primary_pairs[(primary_pairs.task == 'absolute_family') &
                             (primary_pairs.comparator == 'contact_control')].iloc[0]
    assert oriented.left_method == 'topoppi_surface'
    assert oriented.primary_method == 'topoppi_surface'
    assert oriented.mae_difference == pytest.approx(-.5)
    assert oriented.mae_difference_ci_high == pytest.approx(-.5)
    assert oriented.left_macro_mae == pytest.approx(.5)
    assert oriented.right_macro_mae == pytest.approx(1.)
    assert oriented.left_family_balanced_mae == pytest.approx(.5)
    assert oriented.right_family_balanced_mae == pytest.approx(1.)
    assert set(primary_pairs.comparison_scope) == {'core'}
    nonexact = main[(main.task == 'absolute_family') & (main.stratum == 'prodigy_nonexact')]
    assert set(nonexact.n_declared) == {4}
    assert main.loc[(main.task == 'mutation_family') & (main.stratum == 'fresh_family'), 'n'].eq(0).all()
    matched = pd.read_csv(output / 'matched_predictions.csv')
    assert not matched.duplicated(['task', 'sample_id', 'method']).any()
    assert len(matched[(matched.task == 'mutation_variant') & (matched.method == 'ba_raw')]) == len(frames['mutation'])
    known = main[(main.task == 'mutation_variant') & (main.stratum == 'known_complex')]
    assert set(known.n_declared) == {10}
    coverage = pd.read_csv(output / 'coverage.csv')
    assert coverage.n_failed.sum() > 0
    assert coverage.n_base_only.sum() > 0


def test_incomplete_run_requires_explicit_partial_evaluation_and_duplicate_rows_always_fail(tmp_path):
    make_benchmark(tmp_path)
    folder = tmp_path / 'runs/absolute_family/fold_0'
    (folder / 'complete.json').unlink()
    with pytest.raises(ValueError, match='not marked complete'):
        collect_predictions(tmp_path)
    partial = summarize_benchmark(tmp_path, allow_partial=True, n_bootstrap=0)
    assert partial['scope'] == 'Partial benchmark out-of-fold coverage'
    assert partial['completeness_issues']
    assert (tmp_path / 'report/report.json').exists()
    path = folder / 'predictions.csv'
    data = pd.read_csv(path)
    pd.concat([data, data.iloc[[0]]]).to_csv(path, index=False)
    with pytest.raises(ValueError, match='exactly once'):
        collect_predictions(tmp_path, allow_partial=True)


def test_old_test_union_external_file_cannot_pass_full_pool_validation(tmp_path):
    make_benchmark(tmp_path)
    path = tmp_path / 'priors/prodigy.csv'
    pd.read_csv(path).iloc[:2].to_csv(path, index=False)
    with pytest.raises(ValueError, match='declared sample/method rows are missing'):
        collect_predictions(tmp_path)


def test_family_weighting_and_cluster_ci_are_not_row_weighted():
    frame = pd.DataFrame(dict(sample_id=['a1', 'a2', 'b', 'c'], complex_key=['a', 'a', 'b', 'c'],
                              family_group=['f1', 'f1', 'f1', 'f2'], interaction_group=['i1', 'i1', 'i2', 'i3'],
                              observed=[0., 0., 0., 0.], predicted=[0., 2., 5., 9.], usable=True))
    result = metrics(frame)
    assert result['mae'] == pytest.approx(4.)
    assert result['macro_mae'] == pytest.approx(5.)
    assert result['family_balanced_mae'] == pytest.approx(6.)
    reference = frame.copy()
    reference['predicted'] = 0.
    paired = paired_comparison(frame, reference, task='mutation_family', n_bootstrap=100, seed=1)
    assert paired['family_balanced_mae_difference'] == pytest.approx(6.)
    assert paired['bootstrap_unit'] == 'family_group'
    assert paired['bootstrap_n_groups'] == 2
    assert paired['family_balanced_mae_difference_ci_low'] == pytest.approx(3.)
    assert paired['family_balanced_mae_difference_ci_high'] == pytest.approx(9.)
    variant = paired_comparison(frame, reference, task='mutation_variant', n_bootstrap=30)
    assert variant['bootstrap_n_groups'] == 3
    singleton = paired_comparison(frame.iloc[:2], reference.iloc[:2], task='mutation_family', n_bootstrap=30)
    assert np.isnan(singleton['mae_difference_ci_low'])
