import numpy as np
import pandas as pd
from topoppi_hotspot.features import AA, CHEMICAL
from topoppi_hotspot.paired_surface_features import cross_partner_context, paired_vertex_fields


def test_paired_fields_rigid_invariance_and_facing_surfaces():
    vertices = np.array([[0., 0., 0.], [2., 0., 0.]])
    normals = np.tile([0., 0., 1.], (2, 1))
    other = vertices + [0., 0., 1.]
    expected = paired_vertex_fields(vertices, normals, other, -normals)
    np.testing.assert_allclose(expected[:, 0], 1.)
    np.testing.assert_allclose(expected[:, 1], 1.)
    np.testing.assert_allclose(expected[:, 2], np.exp(-.5))
    np.testing.assert_allclose(expected[:, 3], 1.)
    rotation = np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]])
    actual = paired_vertex_fields(vertices @ rotation + 5, normals @ rotation,
                                  other @ rotation + 5, -normals @ rotation)
    np.testing.assert_allclose(actual, expected)


def test_opposite_contact_mean_uses_full_keys_and_atom_weights():
    columns = [f'aa.{a}' for a in AA + 'X'] + list(CHEMICAL)
    frame = pd.DataFrame(0., index=range(3), columns=columns)
    frame['complex_id'] = 'toy'
    frame['residue_key'] = ['c', 'a', 'b']
    frame['side'] = [1, 0, 1]
    channel = list(CHEMICAL)[0]
    frame[channel] = [4., 10., 2.]
    graph = {'toy': {'residue_keys': np.array(['a', 'b', 'c']),
                     'cross_edges': np.array([[0, 1], [0, 2]]),
                     'cross_weights': np.array([1., 3.])}}
    result, _ = cross_partner_context(frame, graph)
    np.testing.assert_allclose(result[f'contact_context.mean.{channel}'], [10., 3.5, 10.])
    np.testing.assert_allclose(result[f'contact_context.product.{channel}'], [40., 35., 20.])
