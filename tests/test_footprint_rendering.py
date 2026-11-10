"""Geometric coverage and annotation/export behavior of native footprint maps."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest
import trimesh
from matplotlib.colors import to_rgba
from PIL import Image
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union

from topoppi.atlas.uv import set_uv_layout
from topoppi.config import VisualizationConfig
from topoppi.io.io_loader import PDBLoader
from topoppi.visualization.export import save_figure
from topoppi.visualization.footprint_rendering import (
    _place_labels,
    footprint_anchor,
    footprint_geometry,
    read_residue_annotations,
    resolve_annotation_values,
)
from topoppi.visualization.visualizer import InterfaceVisualizer

FIXTURE = Path(__file__).parent / "fixtures" / "tiny_complex.pdb"


def square_patch(cut=False, single_residue=False):
    mesh = trimesh.Trimesh(
        vertices=[[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]],
        faces=[[0, 1, 2], [0, 2, 3]], process=False,
    )
    mesh.metadata["source_atom_indices"] = np.asarray([0, 1, 2, 3] if single_residue else [0, 1, 4, 5])
    corners = mesh.vertices[:, :2][mesh.faces].copy()
    if cut:
        corners[1] += [3, 0]
    set_uv_layout(mesh, corners)
    set_uv_layout(mesh, corners, key="uv_global")
    return mesh


def visualizer(**kwargs):
    loader = PDBLoader(FIXTURE)
    coords_a, atoms_a = loader.get_chain_data("A")
    coords_b, atoms_b = loader.get_chain_data("B")
    return InterfaceVisualizer(atoms_a, coords_a, coords_b, atoms_b, chain_a_id="A", chain_b_id="B",
                               config=VisualizationConfig(residue_scope="patch", color_by_interaction_type=False),
                               **kwargs)


def test_equal_area_quadrilaterals_cover_triangles_without_changing_uv():
    mesh = square_patch(cut=True)
    uv = mesh.metadata["uv_global"].copy()
    geometry = footprint_geometry(mesh, uv, ["A:GLY:1", "A:GLY:1", "A:ALA:2", "A:ALA:2"])
    polygons = [Polygon(points) for points in geometry["polygons"]]
    assert len(polygons) == 3 * len(mesh.faces)
    np.testing.assert_allclose([polygon.area for polygon in polygons], np.repeat(1 / 6, 6))
    assert sum(polygon.area for polygon in polygons) == pytest.approx(unary_union(polygons).area)
    expected = unary_union([Polygon(triangle) for triangle in uv])
    assert expected.symmetric_difference(unary_union(polygons)).area < 1e-14
    np.testing.assert_array_equal(uv, mesh.metadata["uv_global"])
    assert len(geometry["boundary_segments"]) == 4
    assert len(geometry["seam_segments"]) == 2
    assert sorted(segment[:, 0].mean() for segment in geometry["seam_segments"]) == [.5, 3.5]
    uncut = square_patch()
    assert len(footprint_geometry(uncut, uncut.metadata["uv"], ["r"] * 4)["seam_segments"]) == 0


def test_disconnected_regions_are_preserved_and_anchor_lies_inside_concave_region():
    mesh = square_patch(cut=True, single_residue=True)
    viz = visualizer()
    figure = viz.plot_patches([mesh], show=False, style_config={"map_style": "footprints"})
    try:
        item = viz.artist_map["1_Gly1"]
        shape = unary_union([Polygon(path.vertices) for path in item["collection"].get_paths()])
        assert shape.geom_type == "MultiPolygon"
        assert len(shape.geoms) == 2
        assert shape.contains(Point(item["anchor"]))
        assert len(item["collection"].get_paths()) == 6
    finally:
        plt.close(figure)
    # The ordinary centroid of this U shape lies in the empty notch.
    shape_parts = np.asarray([[[0, 0], [1, 0], [1, 4], [0, 4]],
                              [[1, 0], [3, 0], [3, 1], [1, 1]],
                              [[3, 0], [4, 0], [4, 4], [3, 4]]])
    assert unary_union([Polygon(part) for part in shape_parts]).contains(Point(footprint_anchor(shape_parts)))


@pytest.mark.parametrize("encoding", ["utf-8", "utf-8-sig"])
def test_csv_author_ids_missing_values_and_alias_errors(tmp_path, encoding):
    keys = {"A:GLY:1", "A:ALA:2", "A:SER:3A"}
    path = tmp_path / "values.csv"
    path.write_text("residue,value\nA:GLY:1,2.3\nA:2,NA\n3A,-0.2\n", encoding=encoding)
    assert read_residue_annotations(path, keys) == {"A:GLY:1": 2.3, "A:ALA:2": None, "A:SER:3A": -.2}
    with pytest.raises(ValueError, match="Duplicate"):
        resolve_annotation_values([("A:GLY:1", 1), ("1", 2)], keys)
    with pytest.raises(ValueError, match="Unknown"):
        resolve_annotation_values([("A:GLY:9", 1)], keys)
    with pytest.raises(ValueError, match="Ambiguous"):
        resolve_annotation_values([("1", 1)], keys | {"B:GLY:1"})
    with pytest.raises(ValueError, match="finite"):
        resolve_annotation_values([("1", "inf")], keys)


def test_annotations_are_self_contained_and_scope_does_not_remove_regions(tmp_path):
    viz = visualizer(interaction_partner_map={"A:GLY:1": {"B:GLY:1": 1}})
    # This fixture has explicit partner evidence but no ProLIF type classification.
    viz.interaction_residue_source = "geometric"
    viz.interaction_partner_map = {"A:GLY:1": {"B:GLY:1": 1}}
    path = tmp_path / "values.csv"
    path.write_text("residue,value\n1,2\n2,NA\n")
    mesh = square_patch()
    figure = viz.plot_patches([mesh], show=False, style_config={
        "map_style": "footprints", "residue_scope": "interaction", "annotation_file": str(path),
    })
    try:
        assert viz.last_report["displayed_residue_count"] == 2
        assert viz.last_report["scope_eligible_residue_count"] == 1
        assert viz.last_report["footprint_polygon_count"] == 6
        assert viz.artist_map["1_Ala2"]["text"] is None
        np.testing.assert_allclose(viz.artist_map["1_Ala2"]["collection"].get_facecolors()[0], to_rgba("#DCE8EF"))
        assert viz.last_style["annotation_values"] == {"A:GLY:1": 2., "A:ALA:2": None}
        saved_style = dict(viz.last_style)
    finally:
        plt.close(figure)
    path.unlink()
    saved_style["residue_scope"] = "patch"
    saved_style["residue_color_overrides"] = {"A:ALA:2": "#ff0000"}
    figure = viz.plot_patches([mesh], show=False, style_config=saved_style)
    try:
        np.testing.assert_allclose(viz.artist_map["1_Ala2"]["collection"].get_facecolors()[0], to_rgba("#D9D9D9"))
        assert viz.last_style["resolved_value_min"] == -2
        assert viz.last_style["resolved_value_max"] == 2
        assert viz.last_report["color_override_count_applied"] == 0
        assert viz.last_style["residue_color_overrides"]["A:ALA:2"] == "#ff0000"
    finally:
        plt.close(figure)


def test_whole_chain_annotation_rows_are_reported_outside_interface():
    viz = visualizer()
    figure = viz.plot_patches([square_patch(single_residue=True)], show=False, style_config={
        "map_style": "footprints", "annotation_values": {"A:GLY:1": 1, "A:ALA:2": 4},
    })
    try:
        assert viz.last_report["outside_domain_residue_count"] == 1
        assert viz.last_style["resolved_value_max"] == 1
        assert viz.last_style["annotation_values"]["A:ALA:2"] == 4
    finally:
        plt.close(figure)
    with pytest.raises(ValueError, match="No annotation"):
        viz.plot_patches([square_patch(single_residue=True)], show=False, style_config={
            "map_style": "footprints", "annotation_values": {"A:ALA:2": 4},
        })


@pytest.mark.parametrize("limits", [(0, 100), (-5, 0), (-1, 3)])
def test_explicit_value_limits_support_signed_and_single_sided_annotations(limits):
    viz = visualizer()
    figure = viz.plot_patches([square_patch()], show=False, style_config={
        "map_style": "footprints", "annotation_values": {"A:GLY:1": 1, "A:ALA:2": 0},
        "value_min": limits[0], "value_max": limits[1],
    })
    try:
        assert viz.last_style["resolved_value_min"] == limits[0]
        assert viz.last_style["resolved_value_max"] == limits[1]
    finally:
        plt.close(figure)


@pytest.mark.parametrize("values, expected", [
    ((0., 1.), (0, 0, "neither")),
    ((-3., 0.), (1, 0, "min")),
    ((1., 3.), (0, 1, "max")),
    ((-3., 3.), (1, 1, "both")),
    ((None, 3.), (0, 1, "max")),
])
def test_colorbar_marks_only_values_beyond_its_limits(values, expected):
    viz = visualizer()
    annotations = dict(zip(("A:GLY:1", "A:ALA:2"), values, strict=True))
    figure = viz.plot_patches([square_patch()], show=False, style_config={
        "map_style": "footprints", "annotation_values": annotations, "value_min": -1., "value_max": 1.,
    })
    try:
        below, above, extend = expected
        assert figure.axes[-1]._colorbar.extend == extend
        assert viz.last_report["below_scale_residue_count"] == below
        assert viz.last_report["above_scale_residue_count"] == above
        assert viz.last_report["colorbar_extend"] == extend
        assert viz.last_style["annotation_values"] == annotations
    finally:
        plt.close(figure)


def test_colorbar_extensions_exclude_residues_outside_annotation_scope():
    viz = visualizer(interaction_partner_map={"A:GLY:1": {"B:GLY:1": 1}})
    viz.interaction_residue_source = "geometric"
    figure = viz.plot_patches([square_patch()], show=False, style_config={
        "map_style": "footprints", "residue_scope": "interaction",
        "annotation_values": {"A:GLY:1": 0., "A:ALA:2": 3.}, "value_min": -1., "value_max": 1.,
    })
    try:
        assert figure.axes[-1]._colorbar.extend == "neither"
        assert viz.last_report["above_scale_residue_count"] == 0
    finally:
        plt.close(figure)


@pytest.mark.parametrize("size", [(7.1, 5.1), (6., 3.3)])
@pytest.mark.parametrize("use_atlas", [True, False])
def test_numeric_footer_stays_inside_resized_canvas_and_png_export(tmp_path, size, use_atlas):
    viz = visualizer()
    annotations = {"A:GLY:1": -.5}
    figure = viz.plot_patches([square_patch()], show=False, style_config={
        "map_style": "footprints", "annotation_values": annotations, "use_uv_atlas": use_atlas,
        "annotation_label": "Measured log2 enrichment", "value_min": -1., "value_max": 1.,
    })
    observed = []

    def record_bounds(event):
        color_axis = figure.axes[-1]
        texts = [color_axis.xaxis.label, *color_axis.get_xticklabels(), *figure.texts]
        observed.append((figure.bbox.frozen(), [text.get_window_extent(event.renderer).frozen()
                                               for text in texts if text.get_visible() and text.get_text()]))

    try:
        figure.canvas.mpl_connect("draw_event", record_bounds)
        figure.set_size_inches(*size)
        figure.canvas.draw()
        destination = tmp_path / "resized-assay.png"
        save_figure(figure, destination)
        assert len(observed) >= 2, "Check both the onscreen-sized draw and the export renderer"
        for canvas, text_bounds in observed:
            for bounds in text_bounds:
                assert bounds.x0 >= canvas.x0 and bounds.x1 <= canvas.x1
                assert bounds.y0 >= canvas.y0 and bounds.y1 <= canvas.y1
        assert destination.stat().st_size > 0
        assert viz.last_style["annotation_values"] == annotations
        assert (viz.last_style["resolved_value_min"], viz.last_style["resolved_value_max"]) == (-1., 1.)
    finally:
        plt.close(figure)


def test_highlights_custom_offsets_and_publication_exports(tmp_path):
    viz = visualizer()
    figure = viz.plot_patches([square_patch()], show=False, style_config={
        "map_style": "footprints", "highlight_residues": ["A:1"], "footprint_labels": "highlighted",
        "label_offsets": {"1_Gly1": [.1, .15]}, "residue_color_overrides": {"A:ALA:2": "#123456"},
    })
    try:
        gly, ala = viz.artist_map["1_Gly1"], viz.artist_map["1_Ala2"]
        np.testing.assert_allclose(gly["collection"].get_facecolors()[0], to_rgba("#A64D79"))
        np.testing.assert_allclose(ala["collection"].get_facecolors()[0], to_rgba("#123456"))
        np.testing.assert_allclose(gly["text"].get_position(), gly["anchor"] + [.1, .15])
        assert ala["text"] is None
        for suffix in ("svg", "pdf", "png", "tiff"):
            save_figure(figure, tmp_path / f"map.{suffix}")
        assert "<text" in (tmp_path / "map.svg").read_text()
        assert b"/FontFile2" in (tmp_path / "map.pdf").read_bytes()
        with Image.open(tmp_path / "map.png") as raster:
            assert raster.info["dpi"][0] == pytest.approx(300, abs=.1)
        with Image.open(tmp_path / "map.tiff") as raster:
            assert raster.info["dpi"] == (600, 600)
            assert raster.info["compression"] == "tiff_lzw"
    finally:
        plt.close(figure)


def test_crowded_label_layout_reports_hidden_labels_without_removing_regions():
    # A supported 14-point label setting on a compact map can exceed its capacity.
    figure, ax = plt.subplots(figsize=(178 / 25.4, 65 / 25.4))
    ax.set(xlim=(0, 1), ylim=(0, 1))
    records = []
    for index in range(80):
        anchor = np.array([.5, .5])
        text = ax.text(*anchor, f"Glu{index}", fontsize=14, ha="center", va="center")
        (connector,) = ax.plot([.5, .5], [.5, .5])
        records.append({"uid": str(index), "residue_key": f"A:GLU:{index}", "anchor": anchor,
                        "text": text, "connector": connector, "highlighted": False})
    try:
        hidden = _place_labels(figure, records, {"avoid_label_overlap": True})
        assert hidden
        assert len(records) == 80
        assert sum(not record["text"].get_visible() for record in records) == len(hidden)
        assert all(np.array_equal(record["anchor"], [.5, .5]) for record in records)
    finally:
        plt.close(figure)
