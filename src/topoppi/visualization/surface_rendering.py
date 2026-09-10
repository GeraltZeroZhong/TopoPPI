"""Native three-dimensional residue surfaces with shared atlas annotations."""

from __future__ import annotations

from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib import patches as mpatches
from matplotlib.cm import ScalarMappable
from matplotlib.transforms import Bbox
from mpl_toolkits.mplot3d import Axes3D, art3d, proj3d

from topoppi.atlas.footprints import mesh_vertex_residue_labels
from topoppi.atlas.seams import uv_seam_topology
from topoppi.atlas.uv import as_corner_uv
from topoppi.mesh.provenance import OPTCUTS_GEOMETRY_VERTEX_IDS, SOURCE_VERTEX_IDS
from topoppi.visualization.export import save_figure
from topoppi.visualization.footprint_rendering import VALUE_CMAP, _resolve_style


def _unit(vectors):
    lengths = np.linalg.norm(vectors, axis=-1, keepdims=True)
    return np.divide(vectors, lengths, out=np.zeros_like(vectors, dtype=float), where=lengths > 0)


def surface_geometry(mesh, vertex_labels):
    """Partition the source triangles into exact three-dimensional residue cells."""
    vertices, faces = np.asarray(mesh.vertices), np.asarray(mesh.faces)
    triangles = vertices[faces]
    labels = np.asarray(vertex_labels)[faces]
    centre = triangles.mean(axis=1)
    area_normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    identity = OPTCUTS_GEOMETRY_VERTEX_IDS if OPTCUTS_GEOMETRY_VERTEX_IDS in mesh.metadata else SOURCE_VERTEX_IDS
    ids = np.asarray(mesh.metadata.get(identity, np.arange(len(vertices))))
    _ids, inverse = np.unique(ids, return_inverse=True)
    normals = np.zeros((len(_ids), 3))
    for k in range(3):
        np.add.at(normals, inverse[faces[:, k]], area_normals)
    normals = _unit(normals)[inverse][faces]
    cells, cell_normals, borders = [], [], []
    for k in range(3):
        j, p = (k + 1) % 3, (k + 2) % 3
        midpoint = (triangles[:, k] + triangles[:, j]) / 2
        cells.append(np.stack([triangles[:, k], midpoint, centre,
                               (triangles[:, k] + triangles[:, p]) / 2], axis=1))
        cell_normals.append(_unit(normals[:, k] + _unit(normals[:, k] + normals[:, j])
                                  + _unit(normals.sum(axis=1)) + _unit(normals[:, k] + normals[:, p])))
        different = labels[:, k] != labels[:, j]
        borders.extend(np.stack([midpoint[different], centre[different]], axis=1))
    uv_key = "uv_global" if "uv_global" in mesh.metadata else "uv"
    topology = uv_seam_topology(mesh, as_corner_uv(mesh, key=uv_key))
    return {
        "triangles": triangles, "corner_labels": labels,
        "cells": np.concatenate(cells), "labels": labels.T.reshape(-1),
        "normals": np.concatenate(cell_normals),
        "residue_borders": np.asarray(borders).reshape(-1, 2, 3),
        "mesh_segments": vertices[topology.edges],
        "boundary_segments": vertices[topology.edges[topology.boundary_mask]],
        "seam_segments": vertices[topology.edges[topology.seam_mask]],
    }


def surface_frame(points, partner_direction):
    """Choose one rigid display frame for all patches, facing the partner."""
    centre = (points.min(axis=0) + points.max(axis=0)) / 2
    basis = np.linalg.svd(points - points.mean(axis=0), full_matrices=False)[2]
    normal, up = basis[-1].copy(), basis[0].copy()
    sign = np.dot(normal, partner_direction)
    if abs(sign) < 1e-12:
        sign = normal[np.argmax(np.abs(normal))]
    if sign < 0:
        normal *= -1
    if up[np.argmax(np.abs(up))] < 0:
        up *= -1
    right = _unit(np.cross(up, normal))
    up = _unit(np.cross(normal, right))
    return centre, np.asarray([right, up, normal])


class ProjectedSurface:
    """Screen-space triangle lookup for visible labels, edges and residue picking."""

    def __init__(self, triangles, matrix):
        self.triangles = triangles
        xyz = triangles.reshape(-1, 3)
        homogeneous = np.column_stack([xyz, np.ones(len(xyz))]) @ matrix.T
        self.w = homogeneous[:, 3].reshape(-1, 3)
        self.projected = (homogeneous[:, :3] / homogeneous[:, 3:]).reshape(-1, 3, 3)
        self.lo = self.projected[:, :, :2].min(axis=(0, 1))
        self.span = np.maximum(np.ptp(self.projected[:, :, :2], axis=(0, 1)), 1e-12)
        self.bins = defaultdict(list)
        self.size = 32
        lower = self._bins(self.projected[:, :, :2].min(axis=1))
        upper = self._bins(self.projected[:, :, :2].max(axis=1))
        for i, (a, b) in enumerate(zip(lower, upper, strict=True)):
            for x in range(a[0], b[0] + 1):
                for y in range(a[1], b[1] + 1):
                    self.bins[x * self.size + y].append(i)

    def _bins(self, points):
        return np.clip(((points - self.lo) / self.span * self.size).astype(int), 0, self.size - 1)

    def query(self, points):
        """Return frontmost face, projected depth and perspective-correct weights."""
        points = np.asarray(points).reshape(-1, 2)
        groups = self._bins(points)
        codes = groups[:, 0] * self.size + groups[:, 1]
        face = np.full(len(points), -1, dtype=int)
        depth = np.full(len(points), np.inf)
        weights = np.zeros((len(points), 3))
        for code in np.unique(codes):
            indices = np.flatnonzero(codes == code)
            candidates = np.asarray(self.bins.get(int(code), []), dtype=int)
            if not len(candidates):
                continue
            t = self.projected[candidates]
            a, b, c = t[:, 0, :2], t[:, 1, :2], t[:, 2, :2]
            denominator = (b[:, 1] - c[:, 1]) * (a[:, 0] - c[:, 0]) + (c[:, 0] - b[:, 0]) * (a[:, 1] - c[:, 1])
            usable = np.abs(denominator) > 1e-18
            denominator = np.where(usable, denominator, 1.)
            for start in range(0, len(indices), 128):
                ii = indices[start:start + 128]
                p = points[ii, None, :] - c[None, :, :]
                u = ((b[:, 1] - c[:, 1]) * p[:, :, 0] + (c[:, 0] - b[:, 0]) * p[:, :, 1]) / denominator
                v = ((c[:, 1] - a[:, 1]) * p[:, :, 0] + (a[:, 0] - c[:, 0]) * p[:, :, 1]) / denominator
                barycentric = np.stack([u, v, 1 - u - v], axis=-1)
                inside = (barycentric.min(axis=-1) >= -1e-8) & usable
                z = np.sum(barycentric * t[None, :, :, 2], axis=-1)
                z = np.where(inside, z, np.inf)
                choice = z.argmin(axis=1)
                selected_depth = z[np.arange(len(ii)), choice]
                hit = np.isfinite(selected_depth)
                selected_faces = candidates[choice[hit]]
                face[ii[hit]] = selected_faces
                depth[ii[hit]] = selected_depth[hit]
                corrected = barycentric[np.arange(len(ii))[hit], choice[hit]] / self.w[selected_faces]
                weights[ii[hit]] = corrected / corrected.sum(axis=1, keepdims=True)
        return face, depth, weights


class SurfaceAxes3D(Axes3D):
    """Refresh occlusion and camera state before every display or export draw."""

    def draw(self, renderer):
        update = getattr(self, "update_surface", None)
        if update is not None:
            update(renderer)
        super().draw(renderer)


def plot_surface(visualizer, patches, style, output_file=None, show=True):
    """Draw all retained patches in their shared three-dimensional coordinate frame."""
    from topoppi.config import VisualizationConfig

    VisualizationConfig(**{key: value for key, value in style.items()
                           if key in VisualizationConfig.__dataclass_fields__}).validate()
    geometries = [surface_geometry(patch, mesh_vertex_residue_labels(patch, visualizer.source_residue_labels_A))
                  for patch in patches]
    domain = set(np.concatenate([g["labels"] for g in geometries]))
    style, norm = _resolve_style(visualizer, style, domain)
    visualizer.last_style = style
    values = style["annotation_values"]
    highlights = set(style["highlight_residues"])
    interaction_colors = {**visualizer.interaction_colors, **style.get("interaction_colors", {})}
    style["interaction_colors"] = interaction_colors
    typed_markers = style["map_style"] == "markers" and style.get("color_by_type", False) and values is None
    geometric_types = (visualizer._geometric_interaction_types()
                       if typed_markers and visualizer.prolif_data is None and visualizer.tree_B is not None else {})
    used_types = set()
    eligible = domain.intersection(visualizer.interaction_partner_map) if style["residue_scope"] == "interaction" else domain
    raw_points = np.concatenate([np.asarray(patch.vertices) for patch in patches])
    direction = (visualizer.coords_B.mean(axis=0) - visualizer.coords_A.mean(axis=0)
                 if len(visualizer.coords_B) and len(visualizer.coords_A) else np.zeros(3))
    centre, rotation = surface_frame(raw_points, direction)

    def transform(xyz):
        return (xyz - centre) @ rotation.T

    triangles = np.concatenate([transform(g["triangles"]) for g in geometries])
    cells = np.concatenate([transform(g["cells"]) for g in geometries])
    normals = np.concatenate([g["normals"] @ rotation.T for g in geometries])
    triangle_gids, cell_gids, records = [], [], []
    fig = plt.figure(figsize=(7, 6.5 if values is not None else 6))
    ax = SurfaceAxes3D(fig, [.025, .13 if values is not None else .025, .95, .84 if values is not None else .95],
                       computed_zorder=False)
    fig.add_axes(ax)
    ax.set_axis_off()
    ax.set_box_aspect((1, 1, 1))
    ax.set_proj_type("ortho" if style["surface_projection"] == "orthographic" else "persp")
    coords = transform(raw_points)
    middle = (coords.min(axis=0) + coords.max(axis=0)) / 2
    radius = max(np.ptp(coords, axis=0).max(), 1e-3) * .53
    default_limits = np.column_stack([middle - radius, middle + radius])

    def set_camera(camera):
        elevation = float(camera.get("elevation", style["surface_elevation"]))
        azimuth = float(camera.get("azimuth", style["surface_azimuth"]))
        roll = float(camera.get("roll", 0.))
        limits = np.asarray(camera.get("limits", default_limits), dtype=float)
        if not np.isfinite([elevation, azimuth, roll]).all() or limits.shape != (3, 2) or not np.isfinite(limits).all() or np.any(limits[:, 1] <= limits[:, 0]):
            raise ValueError("Surface camera requires finite angles and three increasing coordinate limits.")
        if hasattr(ax, "roll"):
            ax.view_init(elev=elevation, azim=azimuth, roll=roll)
        else:
            ax.view_init(elev=elevation, azim=azimuth)
        ax.set_xlim3d(*limits[0])
        ax.set_ylim3d(*limits[1])
        ax.set_zlim3d(*limits[2])

    def capture_camera():
        camera = {"elevation": float(ax.elev), "azimuth": float(ax.azim),
                  "limits": [list(ax.get_xlim3d()), list(ax.get_ylim3d()), list(ax.get_zlim3d())]}
        if hasattr(ax, "roll"):
            camera["roll"] = float(ax.roll)
        visualizer.last_style["surface_camera"] = camera
        return camera

    def reset_camera():
        set_camera({})
        capture_camera()
        fig.canvas.draw_idle()

    set_camera({})
    for _ in range(2):
        projected_points = np.column_stack(proj3d.proj_transform(*coords.T, ax.get_proj()))
        screen = ax.transData.transform(projected_points[:, :2])
        scale = max(np.ptp(screen[:, 0]) / ax.bbox.width, np.ptp(screen[:, 1]) / ax.bbox.height) / .86
        default_limits = middle[:, None] + (default_limits - middle[:, None]) * scale
        set_camera({})
    default_limits = middle[:, None] + (default_limits - middle[:, None]) / style["surface_zoom"]
    set_camera(style.get("surface_camera") or {})
    rgba = []
    for patch_id, g in enumerate(geometries, start=1):
        lookup = {}
        for key in sorted(set(g["labels"])):
            metadata = visualizer.residue_metadata_A[key]
            uid = f"{patch_id}_{visualizer._format_residue_label(metadata['residue_name'], metadata['residue_token'])}"
            lookup[key] = uid
            color = style.get("color", "red") if style["map_style"] == "markers" else style["footprint_color"]
            marker_eligible = key in eligible
            if typed_markers and marker_eligible:
                types = (visualizer.prolif_data.get(metadata["residue_token"], set())
                         if visualizer.prolif_data is not None else geometric_types.get(key, set()))
                active = set(style.get("active_types", visualizer.interaction_types))
                ranked = [(visualizer.interaction_rank[t], t) for t in types if t in active and t in visualizer.interaction_rank]
                if ranked:
                    best_type = min(ranked)[1]
                    color = interaction_colors.get(best_type, color)
                    used_types.add(best_type)
                elif types:
                    marker_eligible = False
            if key in eligible:
                if values is not None:
                    color = style["missing_color"] if values.get(key) is None else VALUE_CMAP(norm(values[key]))
                elif key in highlights:
                    color = style["highlight_color"]
            if values is None:
                if style["map_style"] == "markers":
                    color = style.get("marker_color_overrides", {}).get(uid, color)
                color = style.get("residue_color_overrides", {}).get(key, color)
            group = transform(g["cells"][g["labels"] == key]).mean(axis=1)
            positions = group[np.argsort(np.linalg.norm(group - group.mean(axis=0), axis=1))]
            partners = {visualizer.residue_metadata_B[p]["residue_token"]: n
                        for p, n in visualizer.interaction_partner_map.get(key, {}).items()
                        if p in visualizer.residue_metadata_B}
            label = visualizer._build_label_text(style.get("label_mode", "chain_a"), metadata["residue_token"],
                                                metadata["residue_name"], partners)
            wanted = (key in eligible and style.get("show_labels", True) and style["footprint_labels"] != "none"
                      and (style["footprint_labels"] == "all" or key in highlights))
            if style["map_style"] == "markers":
                wanted = wanted and marker_eligible
            text = ax.text2D(0, 0, label, fontsize=style.get("font_size", 8),
                             fontfamily=style.get("font_family", "sans-serif"), ha="center", va="center",
                             color="#202B33", zorder=5,
                             bbox={"facecolor": "white", "edgecolor": "none", "alpha": .8, "pad": .18}) if wanted else None
            scatter = None
            if style["map_style"] == "markers" and marker_eligible:
                scatter = ax.scatter(*positions[0], color=color, s=26, depthshade=False, edgecolors="white", linewidths=.5, zorder=4)
            records.append({"uid": uid, "key": key, "positions": positions, "text": text, "scatter": scatter,
                            "color": color, "highlighted": key in highlights})
            visualizer.artist_map[uid] = {"residue_key": key, "anchor": positions[0], "text": text,
                                         "scatter": scatter, "connector": None}
        triangle_gids.extend([[lookup[key] for key in row] for row in g["corner_labels"]])
        cell_gids.extend([lookup[key] for key in g["labels"]])
        palette = {record["key"]: record["color"] for record in records}
        rgba.extend([colors.to_rgba(palette[key] if style["map_style"] == "footprints" else style["footprint_color"])
                     for key in g["labels"]])
    rgba = np.asarray(rgba)
    light = _unit(np.array([-.35, -.45, 1.]))
    brightness = .72 + .28 * np.abs(normals @ light)
    lit_colors = rgba.copy()
    lit_colors[:, :3] *= brightness[:, None]
    surface = art3d.Poly3DCollection(cells, facecolors=lit_colors, edgecolors="none", antialiaseds=False, zorder=1)
    ax.add_collection3d(surface)
    for item in visualizer.artist_map.values():
        item["collection"] = surface
    triangle_gids = np.asarray(triangle_gids)
    line_records = []
    for name, enabled, color, width, alpha in (
        ("mesh_segments", style["show_mesh"], "#7D929F", .25, .55),
        ("residue_borders", style["show_residue_borders"], "#748895", .4, .85),
        ("boundary_segments", True, "#536873", .65, 1.),
        ("seam_segments", style["show_seams"], "#202B33", .85, 1.),
    ):
        if enabled:
            segments = np.concatenate([transform(g[name]) for g in geometries])
            # Short subsegments allow folds to occlude only the hidden part of an edge.
            ends = segments[:, :1] + np.linspace(0, 1, 4)[None, :, None] * (segments[:, 1:] - segments[:, :1])
            pieces = np.stack([ends[:, :-1], ends[:, 1:]], axis=2).reshape(-1, 2, 3)
            collection = art3d.Line3DCollection([], colors=color, linewidths=width, alpha=alpha, zorder=2)
            collection.set_gid(name)
            ax.add_collection(collection)
            line_records.append((collection, pieces))
    state = {"projection": None, "index": None}

    def projected(xyz, matrix):
        return np.column_stack(proj3d.proj_transform(*np.asarray(xyz).T, matrix))

    def refresh(renderer=None):
        matrix = ax.get_proj()
        if state["projection"] is None or not np.array_equal(matrix, state["projection"]):
            state["projection"] = matrix.copy()
            state["index"] = ProjectedSurface(triangles, matrix)
            index = state["index"]
            for collection, segments in line_records:
                xyz = projected(segments.mean(axis=1), matrix)
                _faces, depths, _weights = index.query(xyz[:, :2])
                collection.set_segments(segments[xyz[:, 2] <= depths + 2e-5])
            for record in records:
                xyz = projected(record["positions"], matrix)
                _faces, depths, _weights = index.query(xyz[:, :2])
                visible = np.flatnonzero(xyz[:, 2] <= depths + 2e-5)
                record["visible"] = bool(len(visible))
                if len(visible):
                    selected = visible[0]
                    record["screen"] = xyz[selected, :2]
                    if record["scatter"] is not None:
                        record["scatter"]._offsets3d = tuple(record["positions"][selected:selected + 1].T)
                if record["scatter"] is not None:
                    record["scatter"].set_visible(record["visible"])
            capture_camera()
        occupied, hidden = [], []
        for record in sorted(records, key=lambda r: not r["highlighted"]):
            text = record["text"]
            if text is None:
                continue
            visible = record.get("visible", False)
            if visible:
                text.set_visible(True)
                text.set_position(record["screen"])
                if renderer is not None:
                    bbox = text.get_window_extent(renderer).expanded(1.06, 1.12)
                    inside = Bbox.intersection(bbox, ax.bbox)
                    visible = inside is not None and inside.width >= bbox.width * .99 and inside.height >= bbox.height * .99
                    if style.get("avoid_label_overlap", True) and any(bbox.overlaps(old) for old in occupied):
                        visible = False
                    if visible:
                        occupied.append(bbox)
            text.set_visible(visible)
            if not visible:
                hidden.append(record["key"])
        visualizer.last_report.update(displayed_label_count=sum(r["text"] is not None for r in records) - len(hidden),
                                      hidden_label_count=len(hidden), hidden_label_residues=hidden)

    def pick_residue(event):
        if event.inaxes is not ax or event.xdata is None or event.ydata is None:
            return None
        refresh()
        face, _depth, weights = state["index"].query([[event.xdata, event.ydata]])
        return str(triangle_gids[face[0], weights[0].argmax()]) if face[0] >= 0 else None

    def zoom(event):
        if event.inaxes is not ax:
            return
        factor = 1.15 ** (-event.step)
        for getter, setter in ((ax.get_xlim3d, ax.set_xlim3d), (ax.get_ylim3d, ax.set_ylim3d), (ax.get_zlim3d, ax.set_zlim3d)):
            limits = np.asarray(getter())
            setter(*(limits.mean() + (limits - limits.mean()) * factor))
        capture_camera()
        fig.canvas.draw_idle()

    visualizer.last_report = {
        "status": "ok", "view": "surface", "map_style": style["map_style"], "patch_count": len(patches),
        "residue_scope": style["residue_scope"], "patch_residue_count": len(domain),
        "displayed_residue_count": len(domain), "scope_eligible_residue_count": len(eligible),
        "displayed_marker_count": sum(r["scatter"] is not None for r in records),
        "footprint_polygon_count": len(cells),
        "seam_segment_count": sum(len(g["seam_segments"]) for g in geometries),
        "boundary_segment_count": sum(len(g["boundary_segments"]) for g in geometries),
        "annotation_residue_count": len(domain.intersection(values)) if values is not None else 0,
        "outside_domain_residue_count": len(set(values).difference(domain)) if values is not None else 0,
        "missing_value_residue_count": sum(values.get(key) is None for key in eligible) if values is not None else 0,
        "color_by_interaction_type": bool(typed_markers), "projection": style["surface_projection"],
        "chain_interaction_residue_count": len(visualizer.interaction_partner_map),
        "patch_interaction_residue_count": len(domain.intersection(visualizer.interaction_partner_map)),
        "interaction_residue_retention_ratio": len(domain.intersection(visualizer.interaction_partner_map)) / len(visualizer.interaction_partner_map)
        if visualizer.interaction_partner_map else 0.,
        "interaction_type_source": visualizer.interaction_type_source,
        "interaction_residue_source": visualizer.interaction_residue_source,
        "color_override_count_applied": len(domain.intersection(style.get("residue_color_overrides", {})))
        if values is None else 0,
        "below_scale_residue_count": 0, "above_scale_residue_count": 0, "colorbar_extend": "neither",
    }
    if values is not None:
        color_axis = fig.add_axes([.32, .075, .36, .012])
        finite = [values[key] for key in eligible if values.get(key) is not None]
        below, above = any(v < norm.vmin for v in finite), any(v > norm.vmax for v in finite)
        extend = "both" if below and above else "min" if below else "max" if above else "neither"
        visualizer.last_report.update(below_scale_residue_count=sum(v < norm.vmin for v in finite),
                                      above_scale_residue_count=sum(v > norm.vmax for v in finite), colorbar_extend=extend)
        bar = fig.colorbar(ScalarMappable(norm=norm, cmap=VALUE_CMAP), cax=color_axis, orientation="horizontal", extend=extend)
        bar.set_label(style.get("annotation_label", "Value"), fontsize=8)
        bar.ax.tick_params(labelsize=7, length=2)
        if any(values.get(key) is None for key in eligible):
            fig.text(.73, .078, "■ NA", color=style["missing_color"], fontsize=8)
    elif used_types:
        ax.set_position([.025, .10, .95, .875])
        fig.legend(handles=[mpatches.Patch(color=interaction_colors[t], label=t)
                            for t in sorted(used_types, key=lambda t: visualizer.interaction_rank[t])],
                   loc="lower center", ncol=min(4, len(used_types)), frameon=False, fontsize=7)
    ax.update_surface = refresh
    visualizer.capture_surface_camera = capture_camera
    fig._topoppi_surface = {"axis": ax, "capture_camera": capture_camera, "reset_camera": reset_camera,
                            "pick_residue": pick_residue, "rotation": rotation, "centre": centre,
                            "source_triangles": triangles, "cells": cells, "base_colors": rgba,
                            "cell_gids": np.asarray(cell_gids)}
    fig.canvas.mpl_connect("scroll_event", zoom)
    fig.canvas.draw()
    if output_file:
        save_figure(fig, output_file)
    if show:
        plt.show()
    return fig
