import numpy as np
import pytest
import trimesh
from topoppi_hotspot.mesh_diffusion import (
    cotangent_operators,
    geometry_channels,
    implicit_resolvent,
    native_footprint_pooling,
    pool_channels,
)

from topoppi.mesh.provenance import initialize_provenance, tracked_submesh


def _mesh():
    mesh = trimesh.creation.icosphere(subdivisions=1, radius=2.)
    initialize_provenance(mesh)
    mesh.metadata["source_atom_indices"] = np.arange(len(mesh.vertices)) % 2
    return mesh


def test_implicit_diffusion_preserves_constants_and_mass_and_dissipates_energy():
    mesh = _mesh()
    stiffness, mass, _ = cotangent_operators(mesh.vertices, mesh.faces)
    values = np.column_stack([np.full(len(mass), 7.), np.arange(len(mass)) % 5])
    output = implicit_resolvent(stiffness, mass, values)
    original_energy = float(values[:, 1] @ stiffness @ values[:, 1])
    for current in output.values():
        np.testing.assert_allclose(current[:, 0], 7., atol=1e-12)
        np.testing.assert_allclose(mass @ current, mass @ values, atol=1e-10)
        assert float(current[:, 1] @ stiffness @ current[:, 1]) <= original_energy + 1e-10


def test_geometry_and_resolvent_are_invariant_to_rigid_rotation_and_translation():
    mesh = _mesh()
    vertices, faces = np.asarray(mesh.vertices), np.asarray(mesh.faces)
    rotation = trimesh.transformations.rotation_matrix(.8, [1, 2, 3])[:3, :3]
    moved = vertices @ rotation.T + [8, -7, 3]
    stiffness, mass, _ = cotangent_operators(vertices, faces)
    other_stiffness, other_mass, _ = cotangent_operators(moved, faces)
    np.testing.assert_allclose(stiffness.toarray(), other_stiffness.toarray(), atol=1e-13)
    np.testing.assert_allclose(mass, other_mass, atol=1e-13)
    channels = geometry_channels(vertices, faces, stiffness, mass)
    other = geometry_channels(moved, faces, other_stiffness, other_mass)
    np.testing.assert_allclose(channels, other, atol=1e-13)
    solved = implicit_resolvent(stiffness, mass, channels)
    moved_solved = implicit_resolvent(other_stiffness, other_mass, other)
    for heat_time in solved:
        np.testing.assert_allclose(solved[heat_time], moved_solved[heat_time], atol=1e-12)


def test_native_footprint_pooling_uses_source_provenance_and_physical_area():
    mesh = trimesh.Trimesh(vertices=[[0., 0., 0.], [2., 0., 0.], [0., 1., 0.]], faces=[[0, 1, 2]], process=False)
    initialize_provenance(mesh)
    mesh.metadata["source_atom_indices"] = np.asarray([0, 1, 1])
    patch = tracked_submesh(mesh, [0], stage="test")
    # Deliberately reorder patch storage; source identities must still map exactly.
    patch.vertices = np.asarray(patch.vertices)[[2, 0, 1]]
    patch.faces = np.asarray([[1, 2, 0]])
    patch.metadata["source_vertex_ids"] = np.asarray([2, 0, 1])
    pooling = native_footprint_pooling(mesh, [patch], ["R0", "R1"], ["R0", "R1", "missing"])
    mean, std, area = pool_channels(pooling, np.asarray([[1.], [2.], [4.]]))
    np.testing.assert_allclose(area, [1 / 3, 2 / 3, 0])
    np.testing.assert_allclose(mean.ravel(), [1, 3, 0])
    np.testing.assert_allclose(std.ravel(), [0, 1, 0])
    assert area.sum() == pytest.approx(mesh.area)


def test_wavelet_modulus_before_pooling_retains_opposite_responses():
    from scipy import sparse

    pooling = sparse.csr_matrix([[1., 1.]])
    difference = np.asarray([[-2.], [2.]])
    before, _, _ = pool_channels(pooling, np.abs(difference))
    after, _, _ = pool_channels(pooling, difference)
    assert before.item() == 2.
    assert abs(after.item()) == 0.
