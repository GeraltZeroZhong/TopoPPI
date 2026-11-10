import numpy as np
import trimesh
from topoppi_hotspot.surface_rim import rim_distance_fields


def test_planar_rim_depth_has_physical_scale_and_dimensionless_normalization():
    vertices = np.array([[x, y, 0.] for y in range(3) for x in range(3)])
    faces = []
    for y in range(2):
        for x in range(2):
            a = 3 * y + x
            faces.extend([[a, a + 1, a + 3], [a + 1, a + 4, a + 3]])
    fields, report = rim_distance_fields(vertices, faces)
    assert report['n_boundary_vertices'] == 8
    np.testing.assert_allclose(fields[4, :4], [1., 1., .5, .5])
    scaled, _ = rim_distance_fields(vertices * 3 + [7, -5, 4], faces)
    np.testing.assert_allclose(scaled[:, :2], fields[:, :2] * 3)
    np.testing.assert_allclose(scaled[:, 2:4], fields[:, 2:4])


def test_closed_component_is_explicitly_flagged_not_infinite_depth():
    mesh = trimesh.creation.icosphere(subdivisions=0)
    fields, report = rim_distance_fields(mesh.vertices, mesh.faces)
    assert report['n_closed_components'] == 1
    np.testing.assert_allclose(fields[:, :-1], 0.)
    np.testing.assert_allclose(fields[:, -1], 1.)


def test_vertex_only_touch_does_not_combine_component_area_or_rim():
    first = np.array([[x, y, 0.] for y in range(3) for x in range(3)])
    faces = []
    for y in range(2):
        for x in range(2):
            a = 3 * y + x
            faces.extend([[a, a + 1, a + 3], [a + 1, a + 4, a + 3]])
    faces = np.asarray(faces)
    # A second square of four times the area touches only the first's corner.
    second = first * 2 + [2., 2., 0.]
    vertices = np.vstack([first, second[1:]])
    mapping = np.r_[8, np.arange(9, 17)]
    union = np.vstack([faces, mapping[faces]])
    fields, report = rim_distance_fields(vertices, union)
    assert report['n_face_components'] == 2
    assert report['n_shared_component_vertices'] == 1
    np.testing.assert_allclose(fields[4, :4], [1., 1., .5, .5])
    np.testing.assert_allclose(fields[mapping[4], :4], [2., 2., .5, .5])
    # Splitting the shared vertex into equal-coordinate copies changes no field.
    separated, _ = rim_distance_fields(np.vstack([first, second]), np.vstack([faces, faces + 9]))
    np.testing.assert_allclose(fields[:9], separated[:9])
    np.testing.assert_allclose(fields[mapping], separated[9:])
