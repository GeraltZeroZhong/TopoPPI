
import numpy as np
import pytest

from topoppi_affinity.ba_cycle_prior import BACycleEncoder, cycle_score, mpnn_inputs


def test_cycle_sign_sum_and_full_partner_cancellation():
    # Complex NLL changes are [2, 5]; isolated partner changes are [1, 3].
    probs = np.zeros((2, 3, 2, 21))
    probs[:, 0, :, 1] = [-2., -5.]
    probs[:, 1, :, 1] = [-1., -100.]
    probs[:, 2, :, 1] = [-100., -3.]
    value = cycle_score(probs, [-1, 0, 1], [0, 0], [1, 1], [0, 1])
    np.testing.assert_array_equal(value, [3., 3.])
    np.testing.assert_array_equal(cycle_score(probs, [-1, 0, 1], [0, 0], [0, 0], [0, 1]), 0.)


def small_graph():
    coords = np.arange(36, dtype=np.float32).reshape(3, 4, 3)
    return {'aa_index': np.array([0, 1, 2]), 'backbone': coords, 'backbone_mask': np.ones((3, 4), bool),
            'chain': np.array(['B', 'A', 'C']), 'residue_number': np.array([1, 1, 1]),
            'insertion': np.array(['', '', '']), 'side': np.array([0, 0, 1]),
            'residue_key': np.array(['B:CYS:1', 'A:ALA:1', 'C:ASP:1'])}


def test_partner_mask_preserves_internal_chain_and_missing_backbone_fails():
    graph = small_graph()
    inputs = mpnn_inputs(graph)
    assert inputs['order'].tolist() == [1, 0, 2]
    np.testing.assert_array_equal(inputs['side'], [0, 0, 1])
    np.testing.assert_array_equal(inputs['chain_encoding_all'], [1, 2, 3])
    np.testing.assert_array_equal(inputs['residue_idx'], [1, 102, 203])
    graph['backbone_mask'][0, 3] = False
    encoder = BACycleEncoder.__new__(BACycleEncoder)
    with pytest.raises(ValueError, match='Mutation sites lack'):
        encoder.site_log_probs(graph, np.array([0]), mc_samples=2)


