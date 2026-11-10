"""An unlabeled deployment input must reproduce saved predictions and units."""
import json

import joblib
import numpy as np
import pandas as pd
import pytest

from topoppi_affinity.artifacts import fit_frozen, predict_frozen
from topoppi_affinity.thermodynamics import R_KCAL


@pytest.mark.parametrize('target,prefix', [('pkd', 'absolute'), ('ddg', 'mutation')])
def test_refit_and_unlabeled_prediction_roundtrip(tmp_path, target, prefix):
    features = tmp_path / 'features'
    features.mkdir()
    frame = pd.DataFrame({'sample_id': ['a', 'b', 'c', 'd'],
                          'complex_id': ['A', 'A', 'B', 'C'],
                          'family_group': ['f1', 'f1', 'f2', 'f3'],
                          target: [6., 8., 7., 9.] if target == 'pkd' else [-1., 1., 0., 2.]})
    matrix = np.array([[1., 0.], [2., 1.], [1., 1.], [3., 2.]])
    names = np.array(['plm.0', 'surface.0'])
    source = features / f'{target}_plm_surface.npz'
    np.savez(source, X=matrix, sample_ids=frame.sample_id.to_numpy(dtype=str), feature_names=names)
    artifact = tmp_path / 'model.joblib'
    fit_frozen(frame, source, artifact, target=target, representation='plm_surface',
               family='ridge', parameters={'alpha': 1.}, antisymmetric=target == 'ddg')
    output = predict_frozen(artifact, source, temperature_K=310.)
    bundle = joblib.load(artifact)
    np.testing.assert_allclose(output.predicted, bundle['model'].predict(matrix))
    assert json.loads(artifact.with_suffix('.json').read_text())['train_sample_ids'] == frame.sample_id.tolist()
    if target == 'pkd':
        np.testing.assert_allclose(output.predicted_kd_M, 10.**-output.predicted)
        np.testing.assert_allclose(output.predicted_dg_at_reference_kcal_mol,
                                   -R_KCAL * 310. * np.log(10.) * output.predicted)
    else:
        np.testing.assert_allclose(output.predicted_ddg_kcal_mol, -bundle['model'].predict(-matrix))
