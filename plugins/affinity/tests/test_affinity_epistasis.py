"""Epistasis requires matched experimental cycles and a common held-out model."""
import json

import numpy as np
import pandas as pd
import pytest

from topoppi_affinity.epistasis import apply_reviews, candidate_cycles, evaluate_epistasis


def fixture_cycles(*, notes=None):
    rows = []
    for index, (mutations, ddg) in enumerate([('AA1V', 1.), ('GA2D', 2.), ('AA1V,GA2D', 5.), ('AA1V', 1.2)], 1):
        rows.append(dict(sample_id=f'skempi2:{index:05d}', complex_id='test:A:B', sequence_a='AG', sequence_b='GG',
            reference='paper', temperature_K=298., temperature_status='reported', method='SPR', kd_wt_M=1e-8,
            mutations=mutations, n_mutations=len(mutations.split(',')), ddg=ddg, family_group='family1',
            interaction_group='interaction1', fold_variant=0, variant_group=f'variant_{mutations}'))
    raw = pd.DataFrame({'Notes': notes if notes is not None else ['', '', '', 'different pH']})
    return candidate_cycles(pd.DataFrame(rows), raw)


def verified_cycle():
    cycles = fixture_cycles()
    return apply_reviews(cycles, pd.DataFrame([dict(cycle_id='cycle_00000', decision='include',
        review_note='Synthetic fixture conditions matched; select first A measurement', a_sample_id='skempi2:00001')]))


def predictions():
    return pd.DataFrame(dict(sample_id=['skempi2:00001', 'skempi2:00002', 'skempi2:00003'],
                             predicted=[.8, 1.9, 4.8], model_id=['fit1'] * 3, fold=[0] * 3))


def test_candidates_retain_replicates_and_notes_without_verification():
    cycles = fixture_cycles()
    assert len(cycles) == 1
    assert cycles.iloc[0].review_status == 'candidate'
    assert json.loads(cycles.iloc[0].a_ids) == ['skempi2:00001', 'skempi2:00004']
    assert cycles.iloc[0].ambiguous_replicates
    assert cycles.iloc[0].notes_status == 'present_requires_review'
    # Removing notes must not turn the same measurements into verified cycles.
    blank = fixture_cycles(notes=[''] * 4)
    assert blank.iloc[0].notes_status == 'blank_not_verified'
    with pytest.raises(ValueError, match='verified'):
        evaluate_epistasis(blank, predictions(), {('fit1', 0): []})


def test_manual_include_requires_selection_and_exclusions_are_retained():
    review = pd.DataFrame([dict(cycle_id='cycle_00000', decision='include', review_note='Reviewed')])
    with pytest.raises(ValueError, match='select a_sample_id'):
        apply_reviews(fixture_cycles(), review)
    review['decision'] = 'exclude'
    assert apply_reviews(fixture_cycles(), review).iloc[0].review_status == 'excluded'
    assert verified_cycle().iloc[0].a_observed_ddg == 1.


def test_common_model_fold_and_epsilon_zero_baseline():
    rows, summary = evaluate_epistasis(verified_cycle(), predictions(), {('fit1', 0): []}, n_bootstrap=20)
    assert rows.iloc[0].epsilon_observed == 2.
    assert rows.iloc[0].epsilon_predicted == pytest.approx(2.1)
    assert rows.iloc[0].additive_predicted == pytest.approx(2.7)
    assert summary[0]['prediction']['mae'] == pytest.approx(.1)
    assert summary[0]['zero_epsilon']['mae'] == 2.
    assert summary[0]['paired_mae_minus_zero'] == pytest.approx(-1.9)
    assert summary[0]['paired_mae_minus_zero_95_ci'] is None
    split = predictions()
    split.loc[0, 'fold'] = 1
    with pytest.raises(ValueError, match='same model and fold'):
        evaluate_epistasis(verified_cycle(), split, {('fit1', 0): [], ('fit1', 1): []})


def test_training_aliases_and_known_singletons_are_explicit():
    with pytest.raises(ValueError, match='singleton training exposure'):
        evaluate_epistasis(verified_cycle(), predictions(), {('fit1', 0): ['skempi2:00004']})
    rows, summary = evaluate_epistasis(verified_cycle(), predictions(), {('fit1', 0): ['skempi2:00004']},
                                      allow_trained_singletons=True)
    assert rows.iloc[0].evaluation_scope == summary[0]['evaluation_scope'] == 'known_singletons'
    with pytest.raises(ValueError, match='AB or equivalent'):
        evaluate_epistasis(verified_cycle(), predictions(), {('fit1', 0): ['skempi2:00003']},
                          allow_trained_singletons=True)
    with pytest.raises(ValueError, match='Missing training IDs'):
        evaluate_epistasis(verified_cycle(), predictions(), {})


def test_cluster_bootstrap_preserves_common_cycle_units():
    cycles = pd.concat([verified_cycle()] * 2, ignore_index=True)
    cycles.loc[1, ['cycle_id', 'family_group', 'complex_id']] = ['cycle_1', 'family2', 'other:A:B']
    other_predictions = predictions()
    other_predictions['sample_id'] = ['other:A', 'other:B', 'other:AB']
    for state, sample in zip(('a', 'b', 'ab'), other_predictions.sample_id, strict=True):
        cycles.loc[1, f'{state}_sample_id'] = sample
        cycles.loc[1, f'{state}_equivalent_ids'] = json.dumps([sample])
    # Two distinct synthetic cycles under the same fitted model, with identical effects.
    rows, summary = evaluate_epistasis(cycles, pd.concat([predictions(), other_predictions]),
                                      {('fit1', 0): []}, n_bootstrap=100)
    assert len(rows) == 2
    assert np.allclose(summary[0]['paired_mae_minus_zero_95_ci'], [-1.9, -1.9])
