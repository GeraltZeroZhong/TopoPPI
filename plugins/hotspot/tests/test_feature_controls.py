"""Independent synthetic checks for reusable numerical primitives."""
import numpy as np
import pandas as pd
from topoppi_hotspot.feature_controls import control_features, native_raw_columns, partner_permutation
from topoppi_hotspot.features import FEATURE_SETS, NEIGHBOR_FEATURES, neighbor_means
from topoppi_hotspot.partner_features import feature_specifications, prepare_partner_features


def context_fixture():
    frame = pd.DataFrame(0., index=np.arange(8), columns=FEATURE_SETS['self'])
    frame['complex_id'], frame['residue_key'] = 'c', list('abcdefgh')
    frame['side'], frame['wt'] = [0] * 4 + [1] * 4, ['G', 'A', 'G', 'G'] * 2
    frame['contact.partner_distance'] = [2, 3, 4, 9] * 2
    frame['chemical.hydropathy'] = [1, 2, 9, 500, 10, 20, 90, 5000]
    frame['contact.degree'] = np.arange(8) + 1
    native = np.asarray([[0, 1], [1, 2], [4, 5], [5, 6]])
    spatial = np.asarray([[0, 1], [0, 2], [1, 2], [4, 5], [4, 6], [5, 6]])
    for kind, edges in [('topology', native), ('spatial', spatial)]:
        means, count = neighbor_means(frame[NEIGHBOR_FEATURES], edges)
        frame[f'{kind}.neighbor_count'] = count
        for index, source in enumerate(NEIGHBOR_FEATURES):
            frame[f'{kind}.mean.{source}'] = means[:, index]
    # Reverse graph serialization order, retaining noninterface nodes: the
    # transform must align residue keys rather than trust row positions.
    graphs = {'c': {'residue_keys': np.asarray(list('hgfedcba')),
                    'native_neighbor_edges': 7 - native}}
    return prepare_partner_features(frame), graphs

def test_native_source_shuffle_preserves_focal_spatial_and_counts_with_key_alignment():
    frame, graphs = context_fixture()
    result = control_features(frame, graphs, 'native_source_shuffle', 3)
    pd.testing.assert_frame_equal(result[feature_specifications()['relative_spatial']],
                                  frame[feature_specifications()['relative_spatial']])
    np.testing.assert_array_equal(result['topology.neighbor_count'], frame['topology.neighbor_count'])
    source = frame.iloc[partner_permutation(frame, 3)]['chemical.hydropathy'].to_numpy()
    expected = np.asarray([source[1], (source[0] + source[2]) / 2, source[1],
                           source[4], (source[3] + source[5]) / 2, source[4]])
    np.testing.assert_allclose(result['topology.mean.chemical.hydropathy'], expected)
    assert result['topology.mean.chemical.hydropathy'].max() <= 90
    # The unlabelled/native-A row remains in each three-node query population.
    assert len(frame) == 6 and frame.wt.eq('A').sum() == 2
    assert not np.array_equal(result['topology.mean.chemical.hydropathy'], frame['topology.mean.chemical.hydropathy'])

def test_append_shuffle_preserves_joint_full_WT_native_block_and_partner_relative_values():
    frame, graphs = context_fixture()
    result = control_features(frame, graphs, 'native_append_shuffle', 3)
    fixed = feature_specifications()['relative_spatial']
    pd.testing.assert_frame_equal(result[fixed], frame[fixed])
    permutation = partner_permutation(frame, 3)
    native = native_raw_columns()
    np.testing.assert_array_equal(result[native], frame.iloc[permutation][native])
    for positions in frame.groupby('side').indices.values():
        assert set(permutation[positions]) == set(positions)
        np.testing.assert_allclose(result.iloc[positions][native].cov(), frame.iloc[positions][native].cov())
    for column in native:
        for role in ('zscore', 'percentile'):
            key = f'relative.partner.{role}.{column}'
            np.testing.assert_allclose(result[key], frame.iloc[permutation][key], atol=1e-12)

def test_query_permutation_is_identical_when_another_complex_is_added_to_batch():
    frame, graphs = context_fixture()
    preceding = frame.assign(complex_id='before')
    batch = pd.concat([preceding, frame], ignore_index=True)
    batch_graphs = {**graphs, 'before': graphs['c']}
    for kind in ('native_source_shuffle', 'native_append_shuffle'):
        single = control_features(frame, graphs, kind, 3)
        together = control_features(batch, batch_graphs, kind, 3)
        selected = together.loc[together.complex_id.eq('c')].reset_index(drop=True)
        pd.testing.assert_frame_equal(single, selected)
