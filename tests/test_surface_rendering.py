"""Geometric and visibility invariants of native three-dimensional views."""

from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.colors import to_rgba
from mpl_toolkits.mplot3d import proj3d
from test_atlas_io import make_atlas

from topoppi.atlas.footprints import mesh_vertex_residue_labels
from topoppi.visualization.export import save_figure
from topoppi.visualization.surface_rendering import ProjectedSurface, surface_geometry


def test_surface_cells_partition_actual_3d_faces_and_uv_seams_remain_unique():
    patches, viz = make_atlas()
    mesh = patches[0]
    mesh.vertices[:, 2] = [.1, .3, .8, -.2]
    vertices, uv = mesh.vertices.copy(), mesh.metadata["uv_global"].copy()
    labels = mesh_vertex_residue_labels(mesh, viz.source_residue_labels_A)
    geometry = surface_geometry(mesh, labels)
    cells = geometry["cells"]
    areas = sum(np.linalg.norm(np.cross(cells[:, i] - cells[:, 0], cells[:, i + 1] - cells[:, 0]), axis=1) / 2
                for i in (1, 2))
    np.testing.assert_allclose(areas.reshape(3, -1), np.tile(mesh.area_faces / 3, (3, 1)))
    assert len(geometry["boundary_segments"]) == 4
    assert len(geometry["seam_segments"]) == 1
    np.testing.assert_allclose(geometry["seam_segments"][0], vertices[[0, 2]])
    np.testing.assert_array_equal(mesh.vertices, vertices)
    np.testing.assert_array_equal(mesh.metadata["uv_global"], uv)


@pytest.mark.parametrize("projection", ["ortho", "persp"])
def test_projected_lookup_selects_front_face_and_preserves_perspective_weights(projection):
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.set_proj_type(projection)
    ax.view_init(elev=90, azim=-90)
    ax.set(xlim=(-2, 2), ylim=(-2, 2), zlim=(-2, 2))
    triangles = np.array([[[-1., -1., 0.], [1., -1., .1], [0., 1., .2]],
                          [[-1., -1., 1.], [1., -1., 1.1], [0., 1., 1.2]]])
    weights = np.array([.2, .3, .5])
    point = weights @ triangles[1]
    matrix = ax.get_proj()
    projected = np.array(proj3d.proj_transform(*point, matrix))
    index = ProjectedSurface(triangles, matrix)
    face, depth, restored = index.query([projected[:2], [100, 100]])
    assert face.tolist() == [1, -1]
    np.testing.assert_allclose(restored[0], weights, atol=1e-12)
    assert depth[0] == pytest.approx(projected[2])
    plt.close(fig)


def test_surface_pick_follows_visible_residue_cells_and_misses_outside_patch():
    patches, viz = make_atlas()
    fig = viz.plot_patches(patches, show=False, style_config={"view": "surface", "map_style": "footprints"})
    scene = fig._topoppi_surface
    cell = scene["cells"][0].mean(axis=0)
    xy = proj3d.proj_transform(*cell, scene["axis"].get_proj())[:2]
    event = SimpleNamespace(inaxes=scene["axis"], xdata=xy[0], ydata=xy[1])
    assert scene["pick_residue"](event) == scene["cell_gids"][0]
    event.xdata, event.ydata = 100, 100
    assert scene["pick_residue"](event) is None
    plt.close(fig)


def test_csv_residue_colors_match_2d_and_ignore_saved_manual_overrides():
    patches, viz = make_atlas()
    style = {"map_style": "footprints", "annotation_values": {"A:GLY:1": 1.2, "A:ALA:2": None},
             "residue_color_overrides": {"A:GLY:1": "green"}}
    flat = viz.plot_patches(patches, show=False, style_config=style)
    expected = {uid: obj["collection"].get_facecolors()[0] for uid, obj in viz.artist_map.items()}
    plt.close(flat)
    fig = viz.plot_patches(patches, show=False, style_config={**style, "view": "surface"})
    scene = fig._topoppi_surface
    for uid, color in expected.items():
        chosen = scene["base_colors"][scene["cell_gids"] == uid]
        np.testing.assert_allclose(chosen, np.tile(color, (len(chosen), 1)))
    assert viz.last_report["missing_value_residue_count"] == 1
    plt.close(fig)


@pytest.mark.parametrize("projection", ["orthographic", "perspective"])
@pytest.mark.parametrize("size", [(4.78, 2.81), (7.04, 5.08)])
def test_surface_units_remain_visible_on_compact_canvas_and_export(tmp_path, projection, size):
    patches, viz = make_atlas()
    annotations = {"A:GLY:1": -.5, "A:ALA:2": None}
    fig = viz.plot_patches(patches, show=False, style_config={
        "view": "surface", "map_style": "footprints", "surface_projection": projection,
        "annotation_values": annotations, "annotation_label": "Measured log2 enrichment",
        "value_min": -1., "value_max": 1.,
    })
    observed = []
    original_colors = fig._topoppi_surface["base_colors"].copy()

    def record_bounds(event):
        color_axis = fig.axes[-1]
        texts = [color_axis.xaxis.label, *color_axis.get_xticklabels(), *fig.texts]
        observed.append((fig.bbox.frozen(), [text.get_window_extent(event.renderer).frozen()
                                            for text in texts if text.get_visible() and text.get_text()]))

    try:
        fig.canvas.mpl_connect("draw_event", record_bounds)
        fig.set_size_inches(*size)
        fig.canvas.draw()
        target = tmp_path / "compact-surface.png"
        save_figure(fig, target)
        assert len(observed) >= 2
        for canvas, text_bounds in observed:
            for bounds in text_bounds:
                assert bounds.x0 >= canvas.x0 and bounds.x1 <= canvas.x1
                assert bounds.y0 >= canvas.y0 and bounds.y1 <= canvas.y1
        assert target.stat().st_size > 0
        assert viz.last_style["annotation_values"] == annotations
        np.testing.assert_array_equal(fig._topoppi_surface["base_colors"], original_colors)
    finally:
        plt.close(fig)


def test_global_residue_color_applies_to_all_2d_marker_pieces():
    patches, viz = make_atlas(style="markers")
    style = {"residue_color_overrides": {"A:GLY:1": "#00ff00"}, "color_by_type": False}
    fig = viz.plot_patches(patches, show=False, style_config={**style, "view": "surface"})
    np.testing.assert_allclose(viz.artist_map["1_Gly1"]["scatter"].get_facecolors()[0], to_rgba("#00ff00"))
    plt.close(fig)
    fig = viz.plot_patches(patches, show=False, style_config=style)
    pieces = [obj for uid, obj in viz.artist_map.items() if uid.startswith("1_Gly1")]
    assert len(pieces) == 2
    for obj in pieces:
        np.testing.assert_allclose(obj["scatter"].get_facecolors()[0], to_rgba("#00ff00"))
    plt.close(fig)


def test_hidden_labels_reappear_after_camera_returns_to_visible_side():
    patches, viz = make_atlas()
    front = patches[0].copy()
    front.vertices[:, 2] += .3
    front.metadata["source_atom_indices"] = np.full(4, 4)
    patches[0].metadata["source_atom_indices"] = np.full(4, 0)
    fig = viz.plot_patches([patches[0], front], show=False,
                           style_config={"view": "surface", "surface_elevation": 90., "show_mesh": False})
    scene = fig._topoppi_surface
    ax = scene["axis"]
    first = {uid: obj["text"].get_visible() for uid, obj in viz.artist_map.items()}
    assert sum(first.values()) == 1
    ax.view_init(elev=-90, azim=-90)
    fig.canvas.draw()
    second = {uid: obj["text"].get_visible() for uid, obj in viz.artist_map.items()}
    assert sum(second.values()) == 1
    assert second != first
    ax.view_init(elev=90, azim=-90)
    fig.canvas.draw()
    assert {uid: obj["text"].get_visible() for uid, obj in viz.artist_map.items()} == first
    plt.close(fig)
