
import joblib
import numpy as np
import pandas as pd

from topoppi_affinity.attribution_features import permuted_weights, split_blocks, weighted_pool
from topoppi_affinity.benchmark_models import LoadedFeatures
from topoppi_affinity.grouped_fitting import fit_nested, learning_subset


def test_blocks_keep_shared_sequence_once_and_separate_area_from_shape():
    ids = np.array(['a', 'b'])
    seq = LoadedFeatures(ids, np.array([[1., 2.], [3., 4.]]), np.array(['sequence.0', 'sequence.1']))
    names = ['plm.global.sum_absdiff.0', 'plm.contact.sum_absdiff.0']
    names += [f'surface.mean.sum_absdiff.{i}' for i in range(30)]
    names += ['surface.area_weighted_plm.sum_absdiff.0']
    native = LoadedFeatures(ids, np.arange(2*len(names)).reshape(2, -1), np.array(names))
    legacy = LoadedFeatures(ids[::-1], np.array([[20., 99.], [10., 98.]]),
                            np.array(['contact.count', 'surface.legacy']))
    blocks = split_blocks(seq, native, legacy)
    np.testing.assert_array_equal(blocks['S'].X, seq.X)
    np.testing.assert_array_equal(blocks['C'].X[:, 0], [10., 20.])
    assert blocks['A'].X.shape[1] == 6
    assert blocks['K'].X.shape[1] == 24
    assert not any('plm.global' in name for block in blocks.values() for name in block.feature_names)
    assert not set(blocks['A'].feature_names) & set(blocks['K'].feature_names)


def test_weight_controls_preserve_requested_strata_and_partner_swap():
    weights = np.arange(1., 17.)
    side = np.repeat([0, 1], 8)
    contact = np.tile(np.repeat([False, True], 4), 2)
    aa = np.tile([0, 0, 1, 1], 4)
    shuffled = permuted_weights(weights, side, contact, aa, 11, condition_aa=True)
    assert not np.array_equal(weights, shuffled)
    for partner in (0, 1):
        for flag in (False, True):
            for identity in (0, 1):
                mask = (side == partner) & (contact == flag) & (aa == identity)
                np.testing.assert_array_equal(np.sort(weights[mask]), np.sort(shuffled[mask]))
    embedding = np.arange(48.).reshape(16, 3)
    np.testing.assert_array_equal(weighted_pool(embedding, side, weights),
                                  weighted_pool(embedding, 1-side, weights))
    np.testing.assert_allclose(weighted_pool(embedding, side, np.ones(16))[:3],
                               embedding[:8].mean(axis=0)+embedding[8:].mean(axis=0))


def test_nested_fit_ignores_outer_labels_and_scaling_uses_only_training(tmp_path):
    rng = np.random.default_rng(7)
    ids = np.array([f'r{i}' for i in range(16)])
    x = rng.normal(size=(16, 4))
    x[12:] += 100.
    feature = LoadedFeatures(ids, x, np.array(['a', 'b', 'c', 'd']))
    frame = pd.DataFrame(dict(sample_id=ids, complex_id=ids, family_group=np.repeat(np.arange(8), 2),
                              pkd=5.+x[:, 0]))
    train, heldout = frame.iloc[:12], frame.iloc[12:].copy()
    heldout['pkd'] = 'forbidden outer target'
    result, info = fit_nested(train, heldout, feature, tmp_path/'first', alphas=(10., 100.), n_splits=3)
    heldout['pkd'] = -1e20
    repeated, _ = fit_nested(train, heldout, feature, tmp_path/'repeat', alphas=(10., 100.), n_splits=3)
    np.testing.assert_allclose(result.predicted, repeated.predicted, rtol=0, atol=0)
    model = joblib.load(tmp_path/'first/model.joblib')['model']
    np.testing.assert_allclose(model.scaler.mean_, x[:12].mean(axis=0))
    assignment = pd.read_csv(tmp_path/'first/inner_assignment.csv')
    assert assignment.groupby('family_group').inner_fold.nunique().max() == 1
    assert info['fit_count'] == 7 and info['heldout_targets_read'] is False


def test_learning_subsets_are_nested_and_keep_whole_families():
    frame = pd.DataFrame({'sample_id': range(80), 'family_group': np.repeat(np.arange(20), 4)})
    subsets = [learning_subset(frame, f, 29) for f in (.25, .5, .75, 1)]
    for left, right in zip(subsets, subsets[1:], strict=False):
        assert set(left.sample_id) <= set(right.sample_id)
        assert left.groupby('family_group').size().eq(4).all()
    assert set(subsets[-1].sample_id) == set(frame.sample_id)







