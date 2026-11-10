import numpy as np
import pandas as pd
import trimesh
from topoppi.mesh.provenance import OPTCUTS_GEOMETRY_VERTEX_IDS, SOURCE_ATOM_IDS

from topoppi_affinity.uv_layout_features import (
    _mass_check,
    chart_ids,
    intrinsic_distances,
    local_correlations,
    sample_surface,
    uv_distances,
)


def square():
    mesh = trimesh.Trimesh(vertices=[[0., 0., 0.], [4., 0., 0.], [4., 4., 0.], [0., 4., 0.]],
                           faces=[[0, 1, 2], [0, 2, 3]], process=False)
    mesh.metadata[SOURCE_ATOM_IDS] = np.arange(4)
    mesh.metadata[OPTCUTS_GEOMETRY_VERTEX_IDS] = np.arange(4)
    return mesh, np.asarray(mesh.vertices)[np.asarray(mesh.faces), :2]


def test_w_mass_is_cut_invariant_but_spatial_readout_consumes_uv():
    mesh, uv = square()
    labels = np.array(['A:ALA:1', 'A:GLY:2', 'A:ALA:1', 'A:TYR:3'])
    expected = {'A:ALA:1': 32/3, 'A:GLY:2': 8/3, 'A:TYR:3': 8/3}
    cut = uv.copy()
    cut[1] += [12, 0]
    assert len(np.unique(chart_ids(mesh, uv))) == 1
    assert len(np.unique(chart_ids(mesh, cut))) == 2
    before = _mass_check(mesh, uv, labels, expected)
    after = _mass_check(mesh, cut, labels, expected)
    assert before['piece_count'] < after['piece_count']
    assert after['max_residue_mass_error'] < 1e-12
    samples = sample_surface(mesh, 128, 42)
    channels = np.column_stack((samples['xyz'][:, 0], samples['xyz'][:, 1]**2))
    original = local_correlations(channels, uv_distances(mesh, uv, samples)[0], samples['mass'], [2., 4.])
    changed = local_correlations(channels, uv_distances(mesh, cut, samples)[0], samples['mass'], [2., 4.])
    assert np.max(np.abs(changed-original)) > .01


def test_uv_readout_ignores_rigid_motion_scale_and_island_packing():
    mesh, uv = square()
    cut = uv.copy()
    cut[1] += 20.
    samples = sample_surface(mesh, 128, 42)
    reference = uv_distances(mesh, cut, samples)[0]
    rotation = np.array([[0., -1.], [1., 0.]])
    moved = 7.3 * (cut @ rotation)
    moved[0] += [270., -55.]
    moved[1] += [-700., 870.]
    actual = uv_distances(mesh, moved, samples)[0]
    np.testing.assert_allclose(actual, reference, atol=1e-10)
    same_chart = samples['face_indices'][:, None] == samples['face_indices'][None, :]
    assert np.isinf(actual[~same_chart]).all()


def test_geodesic_stitches_diskification_copies():
    mesh, _ = square()
    sample = sample_surface(mesh, 80, 12)
    original = intrinsic_distances(mesh, sample)
    duplicated = trimesh.Trimesh(vertices=np.asarray(mesh.triangles).reshape(-1, 3),
                                 faces=np.arange(6).reshape(2, 3), process=False)
    duplicated.metadata[OPTCUTS_GEOMETRY_VERTEX_IDS] = np.asarray(mesh.faces).reshape(-1)
    np.testing.assert_allclose(intrinsic_distances(duplicated, sample), original)
    assert np.isfinite(original).all()
    np.testing.assert_allclose(original, original.T)
    np.testing.assert_array_equal(np.diag(original), np.zeros(len(original)))


def test_uv_permutation_preserves_mass_and_channels_but_breaks_correspondence():
    mesh, uv = square()
    samples = sample_surface(mesh, 100, 42)
    z = np.c_[samples['xyz'][:, :2], samples['xyz'][:, 0]**2]
    before = z.copy()
    original, points, charts = uv_distances(mesh, uv, samples)
    changed, moved, new_charts = uv_distances(mesh, uv, samples, permutation_seed=17)
    np.testing.assert_array_equal(z, before)
    np.testing.assert_array_equal(np.sort(points, axis=0), np.sort(moved, axis=0))
    np.testing.assert_array_equal(np.sort(charts), np.sort(new_charts))
    assert not np.allclose(original, changed)
    a, b = [local_correlations(z, d, samples['mass'], [2., 4.]) for d in (original, changed)]
    assert not np.allclose(a, b)


def fake_cohort(n=30):
    return pd.DataFrame(dict(sample_id=[f's{i}' for i in range(n)], complex_id=[f'c{i}' for i in range(n)],
                              source='fixture', family_group=[f'g{i//2}' for i in range(n)],
                              structure_path='/fixture.pdb', chains_a='A', chains_b='B', pdb_id='fixture',
                              pkd=np.arange(n)))








