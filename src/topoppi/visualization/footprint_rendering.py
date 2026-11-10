"""Residue regions on the unmodified, per-corner TopoPPI UV atlas.

Each corner owns the quadrilateral between that corner, its two incident edge
midpoints, and the triangle barycentre. The three regions partition a triangle
into equal areas; disconnected residue regions remain separate paths.
"""

from __future__ import annotations

import csv
import logging
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox, blended_transform_factory
from shapely.geometry import Polygon
from shapely.ops import polylabel, unary_union

from topoppi.atlas.footprints import mesh_vertex_residue_labels
from topoppi.atlas.seams import uv_seam_topology
from topoppi.atlas.uv import as_corner_uv
from topoppi.mesh.provenance import OPTCUTS_GEOMETRY_VERTEX_IDS, SOURCE_VERTEX_IDS
from topoppi.visualization.export import save_figure

logger = logging.getLogger("Visualizer")
VALUE_CMAP = colors.LinearSegmentedColormap.from_list(
    "topoppi_value", ["#2B6F9C", "#FAF9F6", "#996029"]
)


def resolve_residue_keys(values, residue_keys):
    """Resolve full author keys or unambiguous chain:sequence / sequence aliases."""
    keys = {str(key) for key in residue_keys}
    aliases = defaultdict(set)
    for key in keys:
        chain, _name, token = key.split(":", 2)
        for alias in (key, f"{chain}:{token}", token):
            aliases[alias].add(key)
    result = []
    for value in values:
        token = str(value).strip()
        matches = aliases.get(token, set())
        if not matches:
            raise ValueError(f"Unknown residue {token!r}; use an author key such as A:GLU:37.")
        if len(matches) != 1:
            raise ValueError(f"Ambiguous residue {token!r}; use its full chain:name:sequence key.")
        result.append(next(iter(matches)))
    return result


def _annotation_value(value, residue):
    if value is None or str(value).strip().lower() in {"", "na", "nan", "n/a"}:
        return None
    try:
        number = float(value)
    except (ValueError, TypeError) as error:
        raise ValueError(f"Annotation for {residue} must be a number or NA, got {value!r}.") from error
    if not np.isfinite(number):
        raise ValueError(f"Annotation for {residue} must be finite or NA, got {value!r}.")
    return number


def resolve_annotation_values(items, residue_keys):
    """Validate annotation rows while preserving explicit missing values."""
    result = {}
    items = list(items)
    resolved = resolve_residue_keys([key for key, _value in items], residue_keys)
    for key, (_input_key, value) in zip(resolved, items, strict=True):
        if key in result:
            raise ValueError(f"Duplicate annotation for residue {key}.")
        result[key] = _annotation_value(value, key)
    return result


def read_residue_annotations(path, residue_keys=None):
    """Read a CSV with residue,value columns and resolve source-chain author IDs."""
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not {"residue", "value"}.issubset(reader.fieldnames):
            raise ValueError("Annotation CSV must contain residue,value columns.")
        items = [(row["residue"], row["value"]) for row in reader]
    if not items:
        raise ValueError("Annotation CSV contains no residue rows.")
    if residue_keys is None:
        result = {}
        for key, value in items:
            key = str(key).strip()
            if not key:
                raise ValueError("Annotation CSV contains an empty residue key.")
            if key in result:
                raise ValueError(f"Duplicate annotation for residue {key}.")
            result[key] = _annotation_value(value, key)
        return result
    return resolve_annotation_values(items, residue_keys)


def footprint_geometry(mesh, uv, vertex_labels):
    """Return exact quadrilaterals, residue borders and both sides of true seams."""
    corners = as_corner_uv(mesh, uv)
    if not np.isfinite(corners).all():
        raise ValueError("Footprint rendering requires finite UV coordinates.")
    labels = np.asarray(vertex_labels, dtype=object)[np.asarray(mesh.faces)]
    centre = corners.mean(axis=1)
    polygons, polygon_labels, borders = [], [], []
    for k in range(3):
        next_k, prev_k = (k + 1) % 3, (k + 2) % 3
        midpoint = (corners[:, k] + corners[:, next_k]) / 2
        polygons.append(np.stack([corners[:, k], midpoint, centre,
                                  (corners[:, k] + corners[:, prev_k]) / 2], axis=1))
        polygon_labels.extend(labels[:, k])
        different = labels[:, k] != labels[:, next_k]
        borders.extend(np.stack([midpoint[different], centre[different]], axis=1))

    topology = uv_seam_topology(mesh, corners)
    identity = OPTCUTS_GEOMETRY_VERTEX_IDS if OPTCUTS_GEOMETRY_VERTEX_IDS in mesh.metadata else SOURCE_VERTEX_IDS
    source_ids = np.asarray(mesh.metadata.get(identity, np.arange(len(mesh.vertices))))
    boundary_edges = {tuple(edge) for edge in topology.source_edges[topology.boundary_mask]}
    seam_edges = {tuple(edge) for edge in topology.source_edges[topology.seam_mask]}
    boundary_segments, seam_segments = [], []
    for face, triangle in zip(np.asarray(mesh.faces), corners, strict=True):
        for k in range(3):
            next_k = (k + 1) % 3
            edge = tuple(sorted((int(source_ids[face[k]]), int(source_ids[face[next_k]]))))
            if edge in boundary_edges:
                boundary_segments.append(triangle[[k, next_k]])
            if edge in seam_edges:
                seam_segments.append(triangle[[k, next_k]])
    return {
        "polygons": np.concatenate(polygons),
        "labels": np.asarray(polygon_labels, dtype=object),
        "residue_borders": np.asarray(borders).reshape(-1, 2, 2),
        "boundary_segments": np.asarray(boundary_segments).reshape(-1, 2, 2),
        "seam_segments": np.asarray(seam_segments).reshape(-1, 2, 2),
    }


def footprint_anchor(polygons):
    """Place an anchor inside the largest connected region, including concave ones."""
    shape = unary_union([Polygon(polygon) for polygon in polygons])
    if shape.geom_type != "Polygon":
        shape = max((part for part in shape.geoms if part.geom_type == "Polygon"), key=lambda part: part.area)
    extent = max(shape.bounds[2] - shape.bounds[0], shape.bounds[3] - shape.bounds[1])
    point = polylabel(shape, tolerance=max(extent / 200, np.finfo(float).eps))
    return np.asarray([point.x, point.y])


def _resolve_style(visualizer, style, domain):
    style = dict(style)
    keys = visualizer.residue_metadata_A
    highlights = style.get("highlight_residues", ())
    if isinstance(highlights, str):
        highlights = highlights.replace(",", " ").split()
    style["highlight_residues"] = sorted(set(resolve_residue_keys(highlights, keys)))
    if style.get("annotation_values") is not None:
        values = resolve_annotation_values(style["annotation_values"].items(), keys)
    elif style.get("annotation_file"):
        values = read_residue_annotations(style["annotation_file"], keys)
    else:
        values = None
    if values is not None and not domain.intersection(values):
        raise ValueError("No annotation CSV residues occur on this interface atlas.")
    style["annotation_values"] = values
    labels = style.get("footprint_labels", "all")
    if labels not in {"all", "highlighted", "none"}:
        raise ValueError("footprint_labels must be all, highlighted, or none.")
    norm = None
    if values is not None:
        finite = [value for key, value in values.items() if key in domain and value is not None]
        extent = max((abs(value) for value in finite), default=1.0) or 1.0
        lower = float(style["value_min"]) if style.get("value_min") is not None else -extent
        upper = float(style["value_max"]) if style.get("value_max") is not None else extent
        if not np.isfinite([lower, upper]).all() or not lower < upper:
            raise ValueError("Value limits must be finite with value_min < value_max.")
        norm = (colors.TwoSlopeNorm(vcenter=0, vmin=lower, vmax=upper) if lower < 0 < upper
                else colors.Normalize(vmin=lower, vmax=upper))
        style["resolved_value_min"], style["resolved_value_max"] = lower, upper
    return style, norm


def _place_labels(fig, records, style):
    """Greedily separate labels in display units, retaining anchors and custom offsets."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    occupied = defaultdict(list)
    hidden = []
    offsets = style.get("label_offsets", {})
    # Reserve explicit placements first, then highlighted labels.
    records = sorted(records, key=lambda rec: (rec["uid"] not in offsets, not rec["highlighted"]))
    for record in records:
        text, ax = record["text"], record["text"].axes
        anchor = record["anchor"]
        anchor_display = ax.transData.transform(anchor)
        bbox = text.get_window_extent(renderer=renderer)
        if not style.get("avoid_label_overlap", True) or record["uid"] in offsets:
            occupied[ax].append(bbox.expanded(1.06, 1.12))
            continue
        width, height = bbox.width * 1.06, bbox.height * 1.14
        candidates = [anchor_display]
        for radius in (1, 1.8, 2.7, 4, 5.5):
            candidates.extend(anchor_display + [np.cos(angle) * width * radius,
                                                 np.sin(angle) * height * radius]
                              for angle in np.linspace(0, 2 * np.pi, 12, endpoint=False))
        placed = False
        for point in candidates:
            candidate = Bbox.from_bounds(point[0] - width / 2, point[1] - height / 2, width, height)
            if not ax.bbox.contains(candidate.x0, candidate.y0) or not ax.bbox.contains(candidate.x1, candidate.y1):
                continue
            if any(candidate.overlaps(other) for other in occupied[ax]):
                continue
            text.set_position(ax.transData.inverted().transform(point))
            occupied[ax].append(candidate)
            placed = True
            break
        if not placed:
            text.set_visible(False)
            record["connector"].set_visible(False)
            hidden.append(record["residue_key"])
            continue
        position = text.get_position()
        record["connector"].set_data([anchor[0], position[0]], [anchor[1], position[1]])
        record["connector"].set_visible(np.linalg.norm(ax.transData.transform(position) - anchor_display) > height * .7)
    if hidden:
        logger.info("Label spacing hid %d labels; all residue regions remain visible.", len(hidden))
    return hidden


def plot_footprints(visualizer, patches, style, output_file=None, show=True):
    """Render complete source-domain regions, independent of annotation scope."""
    use_atlas = bool(style.get("use_uv_atlas", True))
    uv_key = "uv_global" if use_atlas else "uv"
    patch_geometry, domain = [], set()
    for patch in patches:
        uv = as_corner_uv(patch, key=uv_key)
        geometry = footprint_geometry(patch, uv, mesh_vertex_residue_labels(patch, visualizer.source_residue_labels_A))
        patch_geometry.append(geometry)
        domain.update(geometry["labels"])
    style, norm = _resolve_style(visualizer, style, domain)
    visualizer.last_style = style
    values = style["annotation_values"]
    highlights = set(style["highlight_residues"])
    interaction_domain = domain.intersection(visualizer.interaction_partner_map)
    eligible = interaction_domain if style["residue_scope"] == "interaction" else domain
    below_scale = above_scale = 0
    scale_extend = "neither"
    if values is not None:
        displayed_values = [values[key] for key in eligible if values.get(key) is not None]
        below_scale = sum(value < norm.vmin for value in displayed_values)
        above_scale = sum(value > norm.vmax for value in displayed_values)
        scale_extend = "both" if below_scale and above_scale else "min" if below_scale else "max" if above_scale else "neither"

    points = np.concatenate([geometry["polygons"].reshape(-1, 2) for geometry in patch_geometry])
    span = np.ptp(points, axis=0)
    height_mm = np.clip(178 * span[1] / max(span[0], np.finfo(float).eps), 65, 240)
    if values is not None:
        height_mm += 14
    if use_atlas:
        fig, axis = plt.subplots(figsize=(178 / 25.4, height_mm / 25.4))
        axes = [axis]
    else:
        cols = min(3, len(patches))
        rows = int(np.ceil(len(patches) / cols))
        fig, axes = plt.subplots(rows, cols, figsize=(178 / 25.4, 65 * rows / 25.4), squeeze=False)
        axes = axes.ravel().tolist()
    fig.patch.set_facecolor("white")
    records = []
    for patch_id, geometry in enumerate(patch_geometry, start=1):
        ax = axes[0] if use_atlas else axes[patch_id - 1]
        ax.set_facecolor("white")
        for key in sorted(set(geometry["labels"])):
            polygons = geometry["polygons"][geometry["labels"] == key]
            color = style.get("footprint_color", "#DCE8EF")
            if key in eligible:
                if values is not None:
                    value = values.get(key)
                    color = style.get("missing_color", "#D9D9D9") if value is None else VALUE_CMAP(norm(value))
                elif key in highlights:
                    color = style.get("highlight_color", "#A64D79")
            # Preserve the quantitative colour scale when switching a previously
            # hand-coloured atlas to CSV annotations. Saved edits remain available
            # when the user clears the numeric annotations.
            if values is None:
                color = style.get("residue_color_overrides", {}).get(key, color)
            collection = PolyCollection(polygons, facecolors=color, edgecolors="none", antialiaseds=False,
                                        zorder=1, picker=True)
            ax.add_collection(collection)
            label = visualizer._residue_label(key)
            uid = f"{patch_id}_{label}"
            collection.set_gid(uid)
            anchor = footprint_anchor(polygons)
            text, connector = None, None
            visible_label = (key in eligible and style.get("show_labels", True)
                             and style["footprint_labels"] != "none"
                             and (style["footprint_labels"] == "all" or key in highlights))
            if visible_label:
                offset = np.asarray(style.get("label_offsets", {}).get(uid, (0, 0)))
                label_text = visualizer._label_for_residue(style.get("label_mode", "chain_a"), key)
                text = ax.text(*(anchor + offset), label_text, fontsize=style.get("font_size", 8),
                               fontfamily=style.get("font_family", "sans-serif"), ha="center", va="center",
                               color="#202B33", zorder=5, picker=True,
                               bbox={"facecolor": "white", "edgecolor": "none", "alpha": .72, "pad": .16})
                text.set_gid(uid)
                (connector,) = ax.plot([anchor[0], anchor[0] + offset[0]], [anchor[1], anchor[1] + offset[1]],
                                       color="#52616B", linewidth=.45, zorder=4)
                records.append({"uid": uid, "residue_key": key, "anchor": anchor, "text": text,
                                "connector": connector, "highlighted": key in highlights})
            visualizer.artist_map[uid] = {"residue_key": key, "anchor": anchor, "collection": collection,
                                          "scatter": None, "text": text, "connector": connector}
        if style.get("show_residue_borders", True):
            borders = LineCollection(geometry["residue_borders"], colors="#7D929F", linewidths=.3, zorder=2)
            borders.set_gid(f"{patch_id}_residue_borders")
            ax.add_collection(borders)
        boundary = LineCollection(geometry["boundary_segments"], colors="#536873", linewidths=.55, zorder=3)
        boundary.set_gid(f"{patch_id}_source_boundary")
        ax.add_collection(boundary)
        if style.get("show_seams", True):
            seams = LineCollection(geometry["seam_segments"], colors="#202B33", linewidths=.65, zorder=3)
            seams.set_gid(f"{patch_id}_seams")
            ax.add_collection(seams)
        ax.autoscale_view()
        ax.margins(.045)
        ax.set_aspect("equal")
        ax.axis("off")
    for ax in axes[len(patches):] if not use_atlas else []:
        ax.axis("off")
    canvas_height_mm = fig.get_figheight() * 25.4
    fig.subplots_adjust(left=.025, right=.975, top=.975,
                        bottom=.025 + 14 / canvas_height_mm if values is not None else .025,
                        wspace=.045, hspace=.045)
    if values is not None:
        color_axis = fig.add_axes([.34, 8 / canvas_height_mm, .32, 1.5 / canvas_height_mm])
        # GUI canvases resize the figure after creation. Keep the footer in
        # physical units, like its fonts, during both redraws and file export.
        def locate_colorbar(_axes, _renderer):
            current_height_mm = fig.get_figheight() * 25.4
            return Bbox.from_bounds(.34, 8 / current_height_mm, .32, 1.5 / current_height_mm)

        color_axis.set_axes_locator(locate_colorbar)
        colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap=VALUE_CMAP), cax=color_axis,
                               orientation="horizontal", extend=scale_extend)
        colorbar.set_label(style.get("annotation_label", "Value"), fontsize=8, labelpad=2)
        colorbar.ax.tick_params(labelsize=7, length=2, width=.5, pad=2)
        colorbar.outline.set_linewidth(.45)
        if any(values.get(key) is None for key in eligible):
            footer_transform = blended_transform_factory(fig.transFigure, fig.dpi_scale_trans)
            fig.add_artist(Rectangle((.70, 8 / 25.4), 3 / 178, 1.5 / 25.4,
                                     transform=footer_transform, facecolor=style.get("missing_color", "#D9D9D9"),
                                     edgecolor="#536873", linewidth=.4))
            fig.text(.726, 8.75 / 25.4, "NA", transform=footer_transform,
                     fontsize=7, va="center", color="#202B33")
    hidden = _place_labels(fig, records, style)
    visualizer.last_report = {
        "status": "ok", "map_style": "footprints", "patch_count": len(patches),
        "residue_scope": style["residue_scope"], "patch_residue_count": len(domain),
        "displayed_residue_count": len(domain), "scope_eligible_residue_count": len(eligible),
        "displayed_marker_count": 0, "displayed_label_count": len(records) - len(hidden),
        "hidden_label_count": len(hidden), "hidden_label_residues": hidden,
        "footprint_polygon_count": sum(len(g["polygons"]) for g in patch_geometry),
        "seam_segment_count": sum(len(g["seam_segments"]) for g in patch_geometry),
        "boundary_segment_count": sum(len(g["boundary_segments"]) for g in patch_geometry),
        "chain_interaction_residue_count": len(visualizer.interaction_partner_map),
        "patch_interaction_residue_count": len(interaction_domain),
        "interaction_residue_retention_ratio": len(interaction_domain) / len(visualizer.interaction_partner_map)
        if visualizer.interaction_partner_map else 0.,
        "interaction_residue_source": visualizer.interaction_residue_source,
        "interaction_type_source": visualizer.interaction_type_source,
        "color_by_interaction_type": False,
        "annotation_residue_count": len(domain.intersection(values)) if values is not None else 0,
        "outside_domain_residue_count": len(set(values).difference(domain)) if values is not None else 0,
        "missing_value_residue_count": sum(values.get(key) is None for key in eligible) if values is not None else 0,
        "below_scale_residue_count": below_scale,
        "above_scale_residue_count": above_scale,
        "colorbar_extend": scale_extend,
        "color_override_count_applied": len(domain.intersection(style.get("residue_color_overrides", {})))
        if values is None else 0,
    }
    if output_file:
        save_figure(fig, output_file)
    if show:
        plt.show()
    return fig
