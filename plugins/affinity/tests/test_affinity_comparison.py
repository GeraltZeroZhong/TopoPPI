import json

import numpy as np
import pandas as pd
import pytest

from topoppi_affinity.comparison import common_cohort_metrics, paired_mae_difference, read_external


def test_paired_difference_sign_and_group_bootstrap():
    report = paired_mae_difference([1, 1, 2], [1, 1, 2], [3, 3, 4], ['a', 'a', 'b'], repetitions=30)
    assert report['external_minus_topoppi_mae'] == 2
    assert report['paired_group_bootstrap_95_ci'] == [2, 2]
    assert report['n_groups'] == 2


def test_failed_predictions_remain_missing_and_provenance_is_required(tmp_path):
    path = tmp_path / 'method.csv'
    pd.DataFrame({'sample_id': ['one', 'two'], 'predicted': [1., np.nan], 'target': ['ddg']*2,
                  'method': ['published']*2, 'status': ['ok', 'error'], 'error': ['', 'missing residue']}).to_csv(path, index=False)
    path.with_suffix('.json').write_text(json.dumps({'protocol': 'native_pretrained', 'training_overlap': 'known',
                                                   'source_url': 'https://example.org/original-code'}))
    frame, metadata = read_external(path)
    assert frame.valid_prediction.tolist() == [True, False]
    assert np.isnan(frame.predicted.iloc[1])
    assert metadata['training_overlap'] == 'known'
    pd.concat([frame, frame.iloc[:1]]).to_csv(path, index=False)
    with pytest.raises(ValueError, match='unique sample IDs'):
        read_external(path)


def test_common_cohort_excludes_unsupported_rows_for_every_method():
    frame = pd.DataFrame([
        ["single", "a", "one", "c1", 1., 1., 2.],
        ["single", "a", "two", "c2", 10., 0., 0.],
        ["single", "b", "one", "c1", 1., 1., 3.],
    ], columns=["experiment", "method", "sample_id", "complex_id", "observed", "topoppi_prediction", "predicted"])
    result = common_cohort_metrics(frame).set_index("method")
    assert result.n_common.tolist() == [1, 1, 1]
    assert result.n_external_methods.tolist() == [2, 2, 2]
    assert result.loc["TopoPPI", "mae"] == 0
    assert result.loc["a", "mae"] == 1
    assert result.loc["b", "mae"] == 2
