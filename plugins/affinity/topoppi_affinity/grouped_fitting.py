"""Grouped inner selection for absolute-affinity descriptor attribution."""
from __future__ import annotations

import json
import time
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from threadpoolctl import threadpool_limits

from topoppi_affinity.embedding_baselines import fit_frozen_regressor
from topoppi_affinity.thermodynamics import pkd_to_dg, pkd_to_kd

ALPHAS = (10., 100., 1000., 10000., 100000.)
FRACTIONS = (.25, .5, .75, 1.)








def learning_subset(train, fraction, seed):
    """Nested prefixes of whole family components; actual row counts are reported."""
    groups = np.array(sorted(train.family_group.unique()))
    np.random.default_rng(seed).shuffle(groups)
    counts = train.groupby('family_group').size().reindex(groups).to_numpy()
    n = int(np.searchsorted(counts.cumsum(), np.ceil(len(train)*fraction))+1)
    n = min(len(groups), max(3, n))
    return train[train.family_group.isin(groups[:n])].copy()


def fit_nested(train, heldout, features, output, *, alphas=ALPHAS, n_splits=3):
    """The heldout argument may contain arbitrary labels; none are consumed."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    heldout = heldout[['sample_id', 'complex_id']].copy()
    if set(train.sample_id) & set(heldout.sample_id):
        raise ValueError('Training and heldout rows overlap')
    x, xv = features.aligned(train.sample_id), features.aligned(heldout.sample_id)
    good_train, good_test = np.isfinite(x).all(axis=1), np.isfinite(xv).all(axis=1)
    used = train[good_train].reset_index(drop=True)
    x = x[good_train]
    if used.family_group.nunique() < n_splits:
        raise ValueError('Insufficient training families for declared inner folds')
    labels = used.pkd.to_numpy(float)
    groups = used.complex_key if 'complex_key' in used else used.complex_id
    folds = list(GroupKFold(n_splits=n_splits).split(x, groups=used.family_group))
    assignments = np.full(len(used), -1)
    inner_predictions = {float(alpha): np.full(len(used), np.nan) for alpha in alphas}
    with threadpool_limits(limits=1):
        for fold, (fitting, valid) in enumerate(folds):
            assignments[valid] = fold
            for alpha in alphas:
                model = fit_frozen_regressor(x[fitting], labels[fitting], np.asarray(groups)[fitting],
                                            family='ridge', params={'alpha': float(alpha)})
                inner_predictions[float(alpha)][valid] = model.predict(x[valid])
        scores = []
        for alpha, prediction in inner_predictions.items():
            errors = pd.Series(np.abs(prediction-labels))
            scores.append(dict(alpha=alpha, inner_macro_mae=float(errors.groupby(np.asarray(groups)).mean().mean()),
                               inner_mae=float(errors.mean())))
        best = min(scores, key=lambda row: (row['inner_macro_mae'], row['alpha']))['alpha']
        model = fit_frozen_regressor(x, labels, groups, family='ridge', params={'alpha': best})
        prediction = np.full(len(heldout), np.nan)
        if good_test.any():
            prediction[good_test] = model.predict(xv[good_test])
    result = heldout.copy()
    result['predicted'] = prediction
    result['status'] = np.where(good_test & np.isfinite(prediction), 'ok', 'error')
    result['error'] = np.where(result.status.eq('ok'), '', 'Missing/nonfinite feature or prediction')
    result['predicted_kd_M'] = np.nan
    result['predicted_dg_298_15_kcal_mol'] = np.nan
    successful = result.status.eq('ok')
    result.loc[successful, 'predicted_kd_M'] = pkd_to_kd(prediction[successful])
    result.loc[successful, 'predicted_dg_298_15_kcal_mol'] = pkd_to_dg(prediction[successful], 298.15)
    result.to_csv(output/'predictions.csv', index=False)
    used[['sample_id', 'family_group']].assign(inner_fold=assignments).to_csv(output/'inner_assignment.csv', index=False)
    used[['sample_id']].to_csv(output/'training_ids.csv', index=False)
    train.loc[~good_train, ['sample_id']].to_csv(output/'excluded_training_ids.csv', index=False)
    heldout[['sample_id']].to_csv(output/'test_ids.csv', index=False)
    pd.DataFrame(scores).to_csv(output/'selection.csv', index=False)
    joblib.dump(dict(model=model, feature_names=features.feature_names, selected_alpha=best), output/'model.joblib')
    report = dict(n_train_declared=len(train), n_train_usable=len(used), n_train_families=used.family_group.nunique(),
                  n_test_declared=len(heldout), n_test_success=int(result.status.eq('ok').sum()),
                  selected_alpha=best, inner_folds=n_splits, alphas=list(alphas),
                  fit_count=len(alphas)*n_splits+1, selection='pooled inner OOF complex-macro MAE',
                  heldout_targets_read=False, seconds=time.monotonic()-started)
    (output/'fit.json').write_text(json.dumps(report, indent=2)+'\n')
    return result, report










