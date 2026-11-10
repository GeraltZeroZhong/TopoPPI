"""Candidate double-mutant cycles and explicitly reviewed epistasis evaluation.

Metadata matches are candidates, never automatically verified experiments.
Delta-delta G is MT minus WT; epsilon = ddG(AB) - ddG(A) - ddG(B).
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from topoppi_affinity.evaluation import regression_metrics

CONDITION_COLUMNS = [
    'complex_id', 'sequence_a', 'sequence_b', 'reference', 'temperature_K',
    'temperature_status', 'method', 'kd_wt_M',
]
LIMITATIONS = [
    'Matching metadata does not verify pH, buffer, construct background or measurement comparability.',
    'Blank SKEMPI Notes do not establish matched experimental conditions.',
    'Repeated measurements and shared singletons do not form independent cycles.',
    'Use the same fitted model and fold for A, B and AB; separately fitted OOF rows cannot be combined.',
    'All current mutation families have project historical exposure; these are not fresh confirmation data.',
]


def _json(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True)


def candidate_cycles(manifest: pd.DataFrame, raw: pd.DataFrame) -> pd.DataFrame:
    """Match explicit singletons to doubles with identical available conditions.

Each output row is a distinct condition/mutation-set cycle, with all replicate
IDs retained. Raw SKEMPI row numbers use the existing one-based sample ID scheme.
    """
    singletons, doubles = defaultdict(list), defaultdict(list)
    aliases = defaultdict(list)
    records = manifest.to_dict('records')
    for row in records:
        key = tuple(row[column] for column in CONDITION_COLUMNS)
        mutations = tuple(sorted(row['mutations'].split(',')))
        if row['n_mutations'] == 1:
            singletons[key + (mutations[0],)].append(row)
        elif row['n_mutations'] == 2:
            doubles[key + (mutations,)].append(row)
        variant = row.get('variant_group') or (row['complex_id'], row['sequence_a'], row['sequence_b'], mutations)
        aliases[variant].append(row['sample_id'])
    notes = {f'skempi2:{index + 1:05d}': str(row.get('Notes', '')) for index, row in enumerate(raw.to_dict('records'))}
    result = []
    for key, ab in sorted(doubles.items(), key=lambda item: str(item[0])):
        prefix, mutations = key[:-1], key[-1]
        a, b = singletons.get(prefix + (mutations[0],)), singletons.get(prefix + (mutations[1],))
        if not a or not b:
            continue
        row = {column: ab[0][column] for column in CONDITION_COLUMNS}
        row.update(cycle_id=f'cycle_{len(result):05d}', mutation_a=mutations[0], mutation_b=mutations[1],
                   family_group=ab[0]['family_group'], interaction_group=ab[0].get('interaction_group', ''),
                   review_status='candidate', review_note='')
        all_notes = {}
        for state, members in [('a', a), ('b', b), ('ab', ab)]:
            ids = sorted(member['sample_id'] for member in members)
            row[f'{state}_ids'] = _json(ids)
            row[f'{state}_ddg_by_id'] = _json({member['sample_id']: float(member['ddg']) for member in members})
            equivalents = set()
            for member in members:
                variant = member.get('variant_group') or (
                    member['complex_id'], member['sequence_a'], member['sequence_b'],
                    tuple(sorted(member['mutations'].split(','))))
                equivalents.update(aliases[variant])
            row[f'{state}_equivalent_ids'] = _json(sorted(equivalents))
            for sample_id in ids:
                all_notes[sample_id] = notes[sample_id]
        row['ambiguous_replicates'] = any(len(members) > 1 for members in (a, b, ab))
        row['notes_by_id'] = _json(all_notes)
        row['notes_status'] = 'present_requires_review' if any(all_notes.values()) else 'blank_not_verified'
        row['variant_fold_compatible'] = any(
            x['fold_variant'] == y['fold_variant'] == z['fold_variant'] for x in a for y in b for z in ab)
        result.append(row)
    columns = CONDITION_COLUMNS + ['cycle_id', 'mutation_a', 'mutation_b', 'family_group', 'interaction_group',
        'review_status', 'review_note'] + [f'{state}_{suffix}' for state in ('a', 'b', 'ab')
        for suffix in ('ids', 'ddg_by_id', 'equivalent_ids')] + [
        'ambiguous_replicates', 'notes_by_id', 'notes_status', 'variant_fold_compatible']
    return pd.DataFrame(result, columns=columns)


def apply_reviews(candidates: pd.DataFrame, reviews: pd.DataFrame) -> pd.DataFrame:
    """Only explicit include decisions with rationale become verified cycles.

For replicated states the reviewer must select one measurement ID. This avoids
silently averaging potentially different conditions. Unreviewed cycles remain
candidates; exclusion is recorded, not silently discarded.
    """
    reviews = reviews.fillna('')
    if reviews.cycle_id.duplicated().any():
        raise ValueError('Each cycle must have at most one review decision')
    unknown = set(reviews.cycle_id) - set(candidates.cycle_id)
    if unknown:
        raise ValueError(f'Review contains unknown cycle IDs: {sorted(unknown)}')
    reviewed = candidates.copy()
    for state in ('a', 'b', 'ab'):
        reviewed[f'{state}_sample_id'] = ''
        reviewed[f'{state}_observed_ddg'] = np.nan
    lookup = reviews.set_index('cycle_id').to_dict('index')
    for index, cycle in reviewed.iterrows():
        review = lookup.get(cycle.cycle_id, {})
        decision = str(review.get('decision', '')).strip()
        if not decision:
            continue
        if decision not in ('include', 'exclude'):
            raise ValueError(f'{cycle.cycle_id}: decision must be include or exclude')
        note = str(review.get('review_note', '')).strip()
        if not note:
            raise ValueError(f'{cycle.cycle_id}: explicit review_note is required')
        reviewed.loc[index, 'review_note'] = note
        reviewed.loc[index, 'review_status'] = 'excluded' if decision == 'exclude' else 'verified'
        if decision == 'exclude':
            continue
        for state in ('a', 'b', 'ab'):
            ids = json.loads(cycle[f'{state}_ids'])
            selected = str(review.get(f'{state}_sample_id', '')).strip()
            selected = selected or (ids[0] if len(ids) == 1 else '')
            if selected not in ids:
                raise ValueError(f'{cycle.cycle_id}: explicitly select {state}_sample_id from {ids}')
            reviewed.loc[index, f'{state}_sample_id'] = selected
            reviewed.loc[index, f'{state}_observed_ddg'] = json.loads(cycle[f'{state}_ddg_by_id'])[selected]
    return reviewed


def prepare_epistasis(manifest, raw, output, *, review_csv=None):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(manifest, keep_default_na=False)
    candidates = candidate_cycles(frame, pd.read_csv(raw, sep=';', keep_default_na=False))
    candidates.to_csv(output / 'candidate_cycles.csv', index=False)
    template = candidates[['cycle_id']].copy()
    for column in ('decision', 'review_note', 'a_sample_id', 'b_sample_id', 'ab_sample_id'):
        template[column] = ''
    # Keep a filled reviewer file separate from this regenerated blank template.
    template.to_csv(output / 'review_template.csv', index=False)
    reviewed = apply_reviews(candidates, pd.read_csv(review_csv, keep_default_na=False) if review_csv else template)
    reviewed.to_csv(output / 'reviewed_cycles.csv', index=False)
    reviewed[reviewed.review_status.eq('verified')].to_csv(output / 'verified_cycles.csv', index=False)
    report = {'manifest': str(manifest), 'raw': str(raw), 'condition_columns': CONDITION_COLUMNS,
              'candidate_cycles': len(candidates),
              'candidate_double_measurements': sum(len(json.loads(ids)) for ids in candidates.ab_ids),
              'complexes': int(candidates.complex_id.nunique()), 'families': int(candidates.family_group.nunique()),
              'review_status_counts': reviewed.review_status.value_counts().to_dict(),
              'labels_scored': False, 'limitations': LIMITATIONS}
    (output / 'preparation.json').write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    return report


def evaluate_epistasis(cycles, predictions, training_ids, *, allow_trained_singletons=False,
                      cluster_column='family_group', n_bootstrap=2000, seed=42):
    """Evaluate reviewed cycles under a common fitted model/fold.

``predictions`` has sample_id, model_id, fold, predicted. ``training_ids`` maps
each (model_id, fold) tuple to the full set of training sample IDs. Model IDs
must identify a fitted model, not merely an architecture. All aliases of a
state are checked. AB exposure is always rejected; A/B exposure is optional
and explicitly separated as known_singletons, never pooled with blind cycles.
    """
    if cycles.empty or not cycles.review_status.eq('verified').all():
        raise ValueError('Evaluation requires nonempty explicitly verified cycles')
    if predictions.duplicated(['model_id', 'fold', 'sample_id']).any():
        raise ValueError('Duplicate prediction for a sample/model/fold')
    rows = []
    for cycle in cycles.to_dict('records'):
        found = False
        for (model_id, fold), group in predictions.groupby(['model_id', 'fold'], sort=False):
            selected = [cycle[f'{state}_sample_id'] for state in ('a', 'b', 'ab')]
            if not set(selected).issubset(set(group.sample_id)):
                continue
            found = True
            if (model_id, fold) not in training_ids:
                raise ValueError(f'Missing training IDs for {(model_id, fold)}')
            training = set(training_ids[(model_id, fold)])
            overlaps = {state: training.intersection(json.loads(cycle[f'{state}_equivalent_ids']))
                        for state in ('a', 'b', 'ab')}
            if overlaps['ab']:
                raise ValueError(f'{cycle["cycle_id"]}: AB or equivalent variant is in training')
            trained_singletons = bool(overlaps['a'] or overlaps['b'])
            if trained_singletons and not allow_trained_singletons:
                raise ValueError(f'{cycle["cycle_id"]}: singleton training exposure; not a blind cycle')
            values = group.set_index('sample_id').loc[selected, 'predicted'].to_numpy(float)
            observed = np.asarray([cycle[f'{state}_observed_ddg'] for state in ('a', 'b', 'ab')], dtype=float)
            if not np.isfinite(values).all() or not np.isfinite(observed).all():
                raise ValueError('Nonfinite cycle observations or predictions')
            row = {column: cycle[column] for column in ('cycle_id', 'complex_id', 'family_group')}
            row.update(model_id=model_id, fold=fold,
                evaluation_scope='known_singletons' if trained_singletons else 'blind_cycle',
                epsilon_observed=float(observed[2] - observed[0] - observed[1]),
                epsilon_predicted=float(values[2] - values[0] - values[1]),
                zero_epsilon_predicted=0., ab_observed=float(observed[2]), ab_predicted=float(values[2]),
                additive_predicted=float(values[0] + values[1]),
                training_overlap_ids=_json({state: sorted(ids) for state, ids in overlaps.items()}))
            rows.append(row)
        if not found:
            raise ValueError(f'{cycle["cycle_id"]}: A/B/AB lack predictions from the same model and fold')
    result = pd.DataFrame(rows)
    summary = []
    # Separate fitted models/folds and exposure scopes. Pooling ensembles is a
    # distinct declared analysis, never an implicit averaging of these rows.
    for (model_id, fold, scope), part in result.groupby(['model_id', 'fold', 'evaluation_scope']):
        truth, pred = part.epsilon_observed.to_numpy(), part.epsilon_predicted.to_numpy()
        delta = np.abs(pred - truth) - np.abs(truth)
        clusters = part[cluster_column].astype(str).to_numpy()
        unique = np.unique(clusters)
        indices = [np.flatnonzero(clusters == cluster) for cluster in unique]
        interval, macro_interval = None, None
        if n_bootstrap and len(indices) > 1:
            rng = np.random.default_rng(seed)
            values = []
            macro_values = []
            for _ in range(n_bootstrap):
                drawn = rng.integers(len(indices), size=len(indices))
                values.append(float(delta[np.concatenate([indices[index] for index in drawn])].mean()))
                macro_values.append(float(np.mean([delta[indices[index]].mean() for index in drawn])))
            interval = np.quantile(values, [.025, .975]).tolist()
            macro_interval = np.quantile(macro_values, [.025, .975]).tolist()
        summary.append({'model_id': model_id, 'fold': int(fold), 'evaluation_scope': scope,
            'cluster_column': cluster_column, 'prediction': regression_metrics(truth, pred, clusters),
            'zero_epsilon': regression_metrics(truth, np.zeros(len(part)), clusters),
            'paired_mae_minus_zero': float(delta.mean()), 'paired_mae_minus_zero_95_ci': interval,
            'paired_cluster_macro_mae_minus_zero': float(np.mean([delta[index].mean() for index in indices])),
            'paired_cluster_macro_mae_minus_zero_95_ci': macro_interval,
            'n_bootstrap': n_bootstrap, 'seed': seed})
    return result, summary
