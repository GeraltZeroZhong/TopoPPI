"""Derived labels must never use held-out variants or incompatible conditions."""
import numpy as np
import pandas as pd

from topoppi_affinity.variant_augmentation import training_variant_pairs


def test_only_training_endpoints_with_identical_condition_enter_difference_labels():
    data = pd.DataFrame({'sample_id': list('abcdef'), 'source': 'skempi', 'complex_id': 'c', 'family_group': 'f',
                         'sequence_a': 'AAA', 'sequence_b': 'CCC', 'reference': 'study', 'method': 'assay',
                         'temperature_K': [298.] * 5 + [310.], 'kd_wt_M': 1e-8,
                         'mutations': ['A1V', 'A1L', 'A1F', 'A2V', 'A2L', 'A2F'], 'ddg': np.arange(6.)})
    train = np.array([True, True, True, False, False, True])
    before = training_variant_pairs(data, train)
    assert set(before.left_id) | set(before.right_id) == set('abc')
    lookup = data.set_index('sample_id').ddg
    np.testing.assert_allclose(before.ddg, before.right_id.map(lookup) - before.left_id.map(lookup))
    data.loc[~train, 'ddg'] += 1000.
    pd.testing.assert_frame_equal(before, training_variant_pairs(data, train))
