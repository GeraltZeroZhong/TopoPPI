"""Render UV patches with residue-level interaction annotations."""

import logging
from typing import Mapping

import matplotlib.patches as mpatches
import matplotlib.patheffects as patheffects
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from scipy.spatial import cKDTree

from topoppi.atlas.footprints import (
    atom_residue_label,
    geometric_contact_partner_map,
    mesh_vertex_residue_labels,
    residue_footprint_pieces,
    source_atom_residue_labels,
)
from topoppi.atlas.uv import as_corner_uv
from topoppi.config import DEFAULT_CONTACT_DISTANCE_ANGSTROM, VisualizationConfig
from topoppi.interactions.interaction_engine import (
    load_prolif_document,
    normalize_interaction_name,
    residue_sequence_token,
)
from topoppi.interactions.metadata import INTERACTION_COLORS, INTERACTION_TYPES
from topoppi.visualization.export import save_figure
from topoppi.visualization.footprint_rendering import plot_footprints

logger = logging.getLogger("Visualizer")


def select_patches_for_display(patches, visualizer, *, map_style=None, min_points=None, view=None):
    """Select visible patches and return interaction counts for the complete atlas."""
    patches = list(patches)
    map_style = visualizer.config.map_style if map_style is None else map_style
    min_points = visualizer.config.min_points if min_points is None else min_points
    counts = [int(visualizer.count_patch_interaction_residues(patch)) for patch in patches]
    view = getattr(visualizer.config, "view", "atlas") if view is None else view
    if map_style == "footprints" or view == "surface":
        return patches, counts
    displayed = [patch for patch, count in zip(patches, counts, strict=True) if count >= min_points]
    if patches and not displayed:
        raise ValueError(
            f"No interface patch meets the marker display threshold ({min_points} interaction residues). "
            "Lower the threshold or select Residue footprints to display the complete atlas."
        )
    return displayed, counts


class InterfaceVisualizer:
    charged_pos = frozenset({"ARG", "LYS", "HIS"})
    charged_neg = frozenset({"ASP", "GLU"})
    aromatic = frozenset({"PHE", "TYR", "TRP", "HIS"})
    cation_atoms = frozenset({"NZ", "NH1", "NH2", "ND1", "NE2"})
    anion_atoms = frozenset({"OD1", "OD2", "OE1", "OE2", "OXT"})
    polar_atoms = frozenset({"N", "O", "S", "F"})
    backbone_atoms = frozenset({"CA", "C", "O", "N"})

    def __init__(
        self,
        chain_A_atoms,
        chain_A_coords,
        chain_B_coords,
        chain_B_atoms=None,
        chain_a_id=None,
        chain_b_id=None,
        structure_label=None,
        prolif_file=None,
        config: VisualizationConfig | None = None,
        interaction_partner_map: Mapping[str, Mapping[str, int]] | None = None,
        contact_distance_angstrom: float | None = None,
    ):
        """
        Args:
            chain_A_atoms: List of Bio.PDB Atom objects for Chain A
            chain_A_coords: Numpy array (N,3) for Chain A
            chain_B_coords: Numpy array (M,3) for Chain B
            chain_B_atoms: List of Bio.PDB Atom objects for Chain B
            chain_a_id: Chain ID string for A (e.g. "A")
            chain_b_id: Chain ID string for B (e.g. "B")
            structure_label: Short structure name displayed in the atlas title
            prolif_file: Path to the ProLIF interaction JSON
            interaction_partner_map: Resolved ProLIF Chain-A-to-B residue map.
                When omitted, it is reconstructed from ``prolif_file``.
            contact_distance_angstrom: Heavy-atom cutoff used only by the
                explicitly enabled geometric fallback.
        """
        self.config = config or VisualizationConfig()
        self.contact_distance_angstrom = float(
            DEFAULT_CONTACT_DISTANCE_ANGSTROM if contact_distance_angstrom is None else contact_distance_angstrom
        )
        self.atoms_A = chain_A_atoms
        self.coords_A = chain_A_coords
        self.coords_B = chain_B_coords
        self.atoms_B = chain_B_atoms
        self.chain_a_id = str(chain_a_id).strip() if chain_a_id else None
        self.chain_b_id = str(chain_b_id).strip() if chain_b_id else None
        self.structure_label = str(structure_label).strip() if structure_label else None
        self.residue_lookup_B = self._build_residue_lookup(self.atoms_B)
        self.source_residue_labels_A = source_atom_residue_labels(self.atoms_A)
        self.residue_metadata_A = self._build_residue_metadata(self.atoms_A)
        self.residue_metadata_B = self._build_residue_metadata(self.atoms_B)
        self._geometric_types_cache = None
        self.artist_map = {}
        self.capture_surface_camera = None
        self.last_report = {}
        self.last_style = {}

        # --- Interaction categories (canonicalized from ProLIF output) ---
        self.interaction_types = list(INTERACTION_TYPES)
        self.interaction_rank = {name: index for index, name in enumerate(self.interaction_types)}
        self.interaction_colors = dict(INTERACTION_COLORS)
        self.patch_fill_palette = [
            "#5A8DEE",
            "#58C4A3",
            "#F6C667",
            "#E78AC3",
            "#7CC6FE",
            "#B8DE6F",
            "#FFA06B",
            "#A78BFA",
            "#7BDFF2",
            "#FF9AA2",
        ]

        # Data store: {PDB residue token: set([interaction_types]) }
        self.prolif_data = None
        self.prolif_partners = {}

        if prolif_file:
            self._load_prolif_data(prolif_file)

        if self.prolif_data is not None:
            if interaction_partner_map is None:
                interaction_partner_map = self._prolif_partner_label_map()
            if not interaction_partner_map:
                raise ValueError(
                    "ProLIF did not yield any interaction residue pairs that resolve against the selected chains."
                )
        elif self.config.use_geometric_interaction_fallback:
            if interaction_partner_map is None and self.atoms_B:
                interaction_partner_map = geometric_contact_partner_map(
                    self.coords_A,
                    self.atoms_A,
                    self.coords_B,
                    self.atoms_B,
                    distance_cutoff=self.contact_distance_angstrom,
                )
        else:
            interaction_partner_map = {}
        self.interaction_partner_map = interaction_partner_map or {}

        self.tree_B = (
            cKDTree(self.coords_B)
            if self.prolif_data is None and self.config.use_geometric_interaction_fallback
            else None
        )
        self.interaction_type_source = (
            "prolif" if self.prolif_data is not None else "geometric_fallback" if self.tree_B is not None else "none"
        )
        self.interaction_residue_source = self.interaction_type_source

    def _build_residue_lookup(self, atoms):
        lookup = {}
        if not atoms:
            return lookup
        for atom in atoms:
            parent = atom.get_parent()
            token = self._residue_token(parent.get_id())
            if token not in lookup:
                lookup[token] = parent.get_resname()
        return lookup

    def _map_title(self):
        context = None
        if self.chain_a_id and self.chain_b_id:
            context = f"surface {self.chain_a_id} / partner {self.chain_b_id}"
        if self.structure_label and context:
            return f"{self.structure_label} - {context}"
        return self.structure_label or context or "Protein interface map"

    def _build_residue_metadata(self, atoms):
        metadata = {}
        for atom in atoms or ():
            residue = atom.get_parent()
            label = atom_residue_label(atom)
            metadata.setdefault(
                label,
                {
                    "residue_id": residue.get_id(),
                    "residue_token": self._residue_token(residue.get_id()),
                    "residue_name": residue.get_resname(),
                },
            )
        return metadata

    def _prolif_partner_label_map(self):
        labels_a = {data["residue_token"]: label for label, data in self.residue_metadata_A.items()}
        labels_b = {data["residue_token"]: label for label, data in self.residue_metadata_B.items()}
        partner_map = {}
        for token_a, partners in self.prolif_partners.items():
            label_a = labels_a.get(token_a)
            if label_a is None:
                continue
            resolved = {
                labels_b[token_b]: int(count)
                for token_b, count in partners.items()
                if token_b in labels_b and int(count) > 0
            }
            if resolved:
                partner_map[label_a] = resolved
        return partner_map

    @staticmethod
    def _residue_token(residue_id):
        sequence = int(residue_id[1])
        insertion = str(residue_id[2]).strip().upper()
        return f"{sequence}{insertion}"

    def _format_residue_label(self, res_name, res_seq):
        if res_name:
            return f"{str(res_name).title()}{res_seq}"
        return str(res_seq)

    def _build_label_text(self, mode, res_a_seq, res_a_name, partner_counts):
        a_label = self._format_residue_label(res_a_name, res_a_seq)
        if mode == "chain_a":
            return a_label
        if not partner_counts:
            return a_label if mode == "pair" else "N/A"
        best_b_seq = min(
            partner_counts.items(),
            key=lambda item: (-int(item[1]), str(item[0])),
        )[0]
        b_name = self.residue_lookup_B.get(best_b_seq, "UNK")
        b_label = self._format_residue_label(b_name, best_b_seq)
        if mode == "chain_b":
            return b_label
        return f"{a_label}-{b_label}"

    @staticmethod
    def _normalize_residue_scope(value):
        scope = str(value).strip().lower()
        if scope not in {"interaction", "patch"}:
            raise ValueError("residue_scope must be 'interaction' or 'patch'.")
        return scope

    def _load_prolif_data(self, json_path):
        payload, _raw_types, _raw_partners = load_prolif_document(json_path)
        engine = str(payload.get("engine") or "").strip().lower()
        expected_a = str(payload.get("chain_a", "")).strip()
        expected_b = str(payload.get("chain_b", "")).strip()
        if engine and engine != "prolif":
            raise ValueError(f"Interaction JSON engine must be 'prolif', got {engine!r}.")
        if expected_a and self.chain_a_id and expected_a != self.chain_a_id:
            raise ValueError(f"ProLIF chain_a mismatch: {expected_a} != {self.chain_a_id}")
        if expected_b and self.chain_b_id and expected_b != self.chain_b_id:
            raise ValueError(f"ProLIF chain_b mismatch: {expected_b} != {self.chain_b_id}")
        valid_a = {data["residue_token"] for data in self.residue_metadata_A.values()}
        valid_b = {data["residue_token"] for data in self.residue_metadata_B.values()}
        residue_types = {}
        partners = {}
        for item in payload.get("interactions", []):
            if not isinstance(item, dict):
                continue
            token_a = residue_sequence_token(item.get("res_a_seq"))
            token_b = residue_sequence_token(item.get("res_b_seq"))
            if token_a not in valid_a or token_b not in valid_b:
                continue
            interaction_type = normalize_interaction_name(item.get("interaction"))
            if interaction_type:
                residue_types.setdefault(token_a, set()).add(interaction_type)
            counts = partners.setdefault(token_a, {})
            counts[token_b] = counts.get(token_b, 0) + 1
        self.prolif_data = residue_types
        self.prolif_partners = partners
        logger.debug(
            "Loaded ProLIF interactions for %d Chain-%s residues.",
            len(self.prolif_data),
            self.chain_a_id,
        )

    def _get_interaction_type_heuristic(self, atom_A, atom_B, dist):
        if dist > max(
            self.config.vdw_distance,
            self.config.ionic_distance,
            self.config.polar_contact_distance,
        ):
            return None

        res_A = atom_A.get_parent().get_resname()
        res_B = atom_B.get_parent().get_resname()
        name_A = atom_A.get_name()
        name_B = atom_B.get_name()
        elem_A = atom_A.element.upper()
        elem_B = atom_B.element.upper()

        is_ani_A = (res_A in self.charged_neg and name_A in self.anion_atoms) or name_A == "OXT"
        is_ani_B = (res_B in self.charged_neg and name_B in self.anion_atoms) or name_B == "OXT"
        is_cat_A = res_A in self.charged_pos and name_A in self.cation_atoms
        is_cat_B = res_B in self.charged_pos and name_B in self.cation_atoms
        is_aro_A = res_A in self.aromatic and name_A not in self.backbone_atoms
        is_aro_B = res_B in self.aromatic and name_B not in self.backbone_atoms
        opposite_charge = (is_cat_A and is_ani_B) or (is_ani_A and is_cat_B)

        if dist < self.config.ionic_distance and opposite_charge:
            return "Ionic"
        if dist < self.config.ionic_distance and ((is_aro_A and is_cat_B) or (is_cat_A and is_aro_B)):
            return "PiCation"

        polar_pair = elem_A in self.polar_atoms and elem_B in self.polar_atoms
        if dist < self.config.polar_contact_distance and polar_pair:
            return "PolarContact"

        if dist < self.config.vdw_distance:
            nonpolar_A = elem_A in {"C", "S"} and name_A not in self.backbone_atoms
            nonpolar_B = elem_B in {"C", "S"} and name_B not in self.backbone_atoms
            if nonpolar_A and nonpolar_B:
                return "Hydrophobic"
            return "VdWContact"

        return None

    def _geometric_interaction_types(self):
        if self._geometric_types_cache is not None:
            return self._geometric_types_cache
        residue_types = {}
        if self.tree_B is not None and self.atoms_B:
            for atom_index, atom_a in enumerate(self.atoms_A):
                residue_label = atom_residue_label(atom_a)
                nearby_partner_indices = self.tree_B.query_ball_point(
                    self.coords_A[atom_index],
                    r=self.contact_distance_angstrom,
                )
                for atom_b_index in nearby_partner_indices:
                    distance = np.linalg.norm(self.coords_A[atom_index] - self.coords_B[atom_b_index])
                    interaction_type = self._get_interaction_type_heuristic(
                        atom_a,
                        self.atoms_B[atom_b_index],
                        distance,
                    )
                    if interaction_type:
                        residue_types.setdefault(residue_label, set()).add(interaction_type)
        self._geometric_types_cache = residue_types
        return residue_types

    def plot_patches(self, patches, output_file=None, show=True, style_config=None):
        if not patches:
            self.last_report = {
                "status": "empty",
                "patch_count": 0,
                "residue_scope": self._normalize_residue_scope(self.config.residue_scope),
                "interaction_type_source": self.interaction_type_source,
                "interaction_residue_source": self.interaction_residue_source,
                "color_by_interaction_type": bool(self.config.color_by_interaction_type),
            }
            return None
        self.artist_map = {}
        style = {
            "color": "red",
            "font_family": "sans-serif",
            "font_size": 9,
            "color_by_type": self.config.color_by_interaction_type,
            "active_types": self.interaction_types,
            "show_labels": True,
            "label_mode": "chain_a",
            "residue_scope": self.config.residue_scope,
            "avoid_label_overlap": True,
            "use_uv_atlas": True,
            "label_offsets": {},
            "marker_color_overrides": {},
            "mesh_fill_alpha": self.config.mesh_fill_alpha,
            "mesh_line_alpha": self.config.mesh_line_alpha,
        }
        for key, default in {
            "map_style": "markers", "highlight_residues": (), "annotation_file": "",
            "annotation_label": "Value", "value_min": None, "value_max": None,
            "footprint_labels": "all", "show_seams": True, "show_residue_borders": True,
            "footprint_color": "#DCE8EF", "highlight_color": "#A64D79", "missing_color": "#D9D9D9",
            "view": "atlas", "show_mesh": True, "surface_projection": "orthographic",
            "surface_elevation": 73., "surface_azimuth": -90., "surface_zoom": 1.,
        }.items():
            style[key] = getattr(self.config, key, default)
        if style_config:
            style.update(style_config)
        style["residue_scope"] = self._normalize_residue_scope(style["residue_scope"])
        if style["map_style"] not in {"markers", "footprints"}:
            raise ValueError("map_style must be 'markers' or 'footprints'.")
        if style["view"] not in {"atlas", "surface"}:
            raise ValueError("view must be 'atlas' or 'surface'.")
        if style["view"] == "surface":
            from topoppi.visualization.surface_rendering import plot_surface

            return plot_surface(self, patches, style, output_file=output_file, show=show)
        if style["map_style"] == "footprints":
            if self.interaction_residue_source == "none" and style["residue_scope"] == "interaction":
                raise ValueError("Interaction scope requires interaction data; use residue_scope='patch' for all residues.")
            if not style_config or "font_size" not in style_config:
                style["font_size"] = 8
            return plot_footprints(self, patches, style, output_file=output_file, show=show)
        if self.interaction_residue_source == "none" and (
            style["color_by_type"] or style["residue_scope"] == "interaction"
        ):
            raise ValueError(
                "Interaction-residue rendering requires a ProLIF JSON. "
                "Geometric inference is disabled by default and must be enabled explicitly."
            )
        interaction_colors = dict(self.interaction_colors)
        interaction_colors.update(style.get("interaction_colors") or {})
        style["interaction_colors"] = interaction_colors
        self.last_style = style

        n_patches = len(patches)
        use_uv_atlas = bool(style.get("use_uv_atlas", True))
        if use_uv_atlas:
            fig, axes = plt.subplots(1, 1, figsize=(10, 8))
            axes = np.asarray([axes], dtype=object)
        else:
            n_cols = min(3, n_patches)
            n_rows = int(np.ceil(n_patches / n_cols))
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.8 * n_cols, 4.2 * n_rows))
            axes = np.atleast_1d(axes).flatten()
        fig.patch.set_facecolor("white")

        logger.debug(
            "Visualizing %d patches in %s mode.",
            n_patches,
            "atlas" if use_uv_atlas else "separate",
        )

        used_interactions = set()
        patch_residues = set()
        interaction_residues = set()
        eligible_residues = set()
        displayed_residues = set()
        marker_count = 0
        label_count = 0
        for i, patch in enumerate(patches):
            ax = axes[0] if use_uv_atlas else axes[i]
            ax.set_facecolor("#fcfcfd")
            found, annotation = self._draw_single_patch(
                ax,
                patch,
                i + 1,
                style,
                use_uv_atlas=use_uv_atlas,
            )
            used_interactions.update(found)
            patch_residues.update(annotation["patch_residues"])
            interaction_residues.update(annotation["interaction_residues"])
            eligible_residues.update(annotation["eligible_residues"])
            displayed_residues.update(annotation["displayed_residues"])
            marker_count += int(annotation["marker_count"])
            label_count += int(annotation["label_count"])
            if not use_uv_atlas:
                ax.set_title(f"Patch {i + 1}", fontsize=10, weight="semibold", color="#2f3640")
                ax.set_aspect("equal")
                ax.axis("off")

        if not use_uv_atlas:
            for extra_ax in axes[n_patches:]:
                extra_ax.axis("off")

        has_interaction_legend = bool(style["color_by_type"] and used_interactions)
        if has_interaction_legend:
            legend_handles = [
                mpatches.Patch(color=style["interaction_colors"].get(name, "gray"), label=name)
                for name in self.interaction_types
                if name in used_interactions
            ]
            legend = fig.legend(
                handles=legend_handles,
                loc="upper center",
                bbox_to_anchor=(0.5, 0.99),
                ncol=min(len(legend_handles), 5),
                frameon=False,
            )
            # ``tight_layout`` can otherwise push a figure-level legend above
            # the canvas for tall UV atlases. Reserve a deterministic header
            # band for both the legend and the axes title.
            legend.set_in_layout(False)

        if use_uv_atlas:
            axes[0].set_title(self._map_title(), fontsize=13, weight="semibold", color="#2f3640", pad=14)
            axes[0].set_aspect("equal")
            axes[0].axis("off")

        fig.tight_layout(rect=[0, 0, 1, 0.87 if has_interaction_legend else 0.95])
        if output_file:
            save_figure(fig, output_file)
        if show:
            plt.show()
        self.last_report = {
            "status": "ok",
            "map_style": "markers",
            "patch_count": int(n_patches),
            "residue_scope": style["residue_scope"],
            "interaction_residue_source": self.interaction_residue_source,
            "interaction_residue_definition": (
                "Chain A residue present in at least one resolved ProLIF interaction record with Chain B"
                if self.interaction_residue_source == "prolif"
                else (
                    "Chain A residue with any Chain A/B heavy-atom pair at distance <= "
                    f"{self.contact_distance_angstrom:g} Angstrom "
                    "(explicit geometric fallback)"
                )
            ),
            "patch_residue_count": int(len(patch_residues)),
            "chain_interaction_residue_count": int(len(self.interaction_partner_map)),
            "patch_interaction_residue_count": int(len(interaction_residues)),
            "interaction_residue_retention_ratio": (
                float(len(interaction_residues) / len(self.interaction_partner_map))
                if self.interaction_partner_map
                else 0.0
            ),
            "scope_eligible_residue_count": int(len(eligible_residues)),
            "displayed_residue_count": int(len(displayed_residues)),
            "displayed_marker_count": int(marker_count),
            "displayed_label_count": int(label_count),
            "color_by_interaction_type": bool(style["color_by_type"]),
            "interaction_type_source": self.interaction_type_source,
        }
        logger.debug(
            "Rendered %d residue(s) from %d %s interaction residue(s) on a %d-residue patch domain.",
            len(displayed_residues),
            len(interaction_residues),
            self.interaction_residue_source,
            len(patch_residues),
        )
        return fig

    def _draw_single_patch(self, ax, patch, patch_id, style, use_uv_atlas=True):
        found_types = set()
        key = "uv_global" if use_uv_atlas else "uv"
        corners = as_corner_uv(patch, key=key)

        patch_color = self.patch_fill_palette[(patch_id - 1) % len(self.patch_fill_palette)]
        fill = PolyCollection(
            corners,
            facecolors=patch_color,
            edgecolors="none",
            alpha=float(style.get("mesh_fill_alpha", 0.22)),
            zorder=0,
        )
        ax.add_collection(fill)
        edges = np.concatenate(
            [corners[:, [0, 1]], corners[:, [1, 2]], corners[:, [2, 0]]],
            axis=0,
        )
        ax.add_collection(
            LineCollection(
                edges,
                colors="#2f3b52",
                alpha=float(style.get("mesh_line_alpha", 0.60)),
                linewidths=0.70,
                zorder=1,
            )
        )

        all_residue_data = self._collect_patch_residue_data(
            patch,
            corners,
            include_types=style["color_by_type"],
        )
        patch_residues = set(all_residue_data)
        interaction_residues = {label for label, data in all_residue_data.items() if data["is_interaction"]}
        if style["residue_scope"] == "interaction":
            residue_data = {label: data for label, data in all_residue_data.items() if data["is_interaction"]}
        else:
            residue_data = all_residue_data
        eligible_residues = set(residue_data)
        displayed_residues = set()
        marker_count = 0

        label_records = []

        for residue_label, data in residue_data.items():
            types = data["types"]

            best_type = None
            if style["color_by_type"]:
                active_types = set(style.get("active_types", []))
                if types:
                    _, best_type = min(
                        (
                            (self.interaction_rank[name], name)
                            for name in types
                            if name in active_types and name in self.interaction_rank
                        ),
                        default=(len(self.interaction_rank), None),
                    )
                    if best_type is None:
                        continue
                    final_color = style["interaction_colors"].get(best_type, "#FFC0CB")
                    found_types.add(best_type)
                else:
                    final_color = style["color"]
            else:
                final_color = style["color"]

            displayed_residues.add(residue_label)

            res_token = data["residue_token"]
            res_name = data["residue_name"]
            label_text = self._build_label_text(
                style.get("label_mode", "chain_a"), res_token, res_name, data.get("partners", {})
            )
            res_name_for_id = self._format_residue_label(res_name, res_token)
            pieces = data["pieces"]
            for piece_index, piece in enumerate(pieces):
                u_center, v_center = piece["uv_centroid"]
                uid = f"{patch_id}_{res_name_for_id}"
                if len(pieces) > 1:
                    uid += f"__piece_{piece_index + 1}"
                residue_color = style.get("residue_color_overrides", {}).get(residue_label, final_color)
                marker_color = style.get("marker_color_overrides", {}).get(uid, residue_color)

                sc = ax.scatter(
                    u_center,
                    v_center,
                    c=marker_color,
                    s=80,
                    edgecolors="white",
                    zorder=10,
                    picker=5,
                )
                sc.set_gid(uid)
                txt = None
                connector = None
                if style.get("show_labels", True):
                    offsets = style.get("label_offsets", {})
                    if uid in offsets:
                        offset = offsets[uid]
                        text_pos = (u_center + offset[0], v_center + offset[1])
                    else:
                        text_pos = (u_center, v_center + self.config.label_offset)
                    txt = ax.text(
                        text_pos[0],
                        text_pos[1],
                        label_text,
                        fontsize=style["font_size"],
                        fontname=style["font_family"],
                        ha="center",
                        fontweight="bold",
                        color="#2f3640",
                        zorder=20,
                    )
                    txt.set_path_effects([patheffects.withStroke(linewidth=2.5, foreground="white")])
                    txt.set_gid(uid)
                    (connector,) = ax.plot(
                        [u_center, text_pos[0]],
                        [v_center, text_pos[1]],
                        linestyle=(0, (2, 2)),
                        color="dimgray",
                        lw=0.8,
                        alpha=0.8,
                        zorder=6,
                    )
                    label_records.append({"text": txt, "scatter": sc, "connector": connector})
                self.artist_map[uid] = {
                    "residue_key": residue_label,
                    "anchor": np.asarray([u_center, v_center]),
                    "collection": None,
                    "scatter": sc,
                    "text": txt,
                    "connector": connector,
                }
                marker_count += 1

        if style.get("show_labels", True) and style.get("avoid_label_overlap", True) and label_records:
            self._relax_labels(ax, label_records)

        center = corners.reshape(-1, 2).mean(axis=0)
        ax.text(center[0], center[1], f"P{patch_id}", fontsize=8, color="#4b5563", alpha=0.35)
        return found_types, {
            "patch_residues": patch_residues,
            "interaction_residues": interaction_residues,
            "eligible_residues": eligible_residues,
            "displayed_residues": displayed_residues,
            "marker_count": marker_count,
            "label_count": len(label_records),
        }

    def _collect_patch_residue_data(self, patch, uv, include_types=True):
        vertex_labels = mesh_vertex_residue_labels(
            patch,
            self.source_residue_labels_A,
        )
        pieces_by_label = residue_footprint_pieces(
            patch,
            as_corner_uv(patch, uv),
            vertex_labels,
        )
        residue_data = {}
        for label, pieces in pieces_by_label.items():
            metadata = self.residue_metadata_A[label]
            interaction_partners = self.interaction_partner_map.get(str(label), {})
            partner_tokens = {}
            for partner_label, interaction_count in interaction_partners.items():
                partner_metadata = self.residue_metadata_B.get(partner_label)
                if partner_metadata is None:
                    continue
                partner_tokens[partner_metadata["residue_token"]] = int(interaction_count)
            residue_data[str(label)] = {
                "pieces": pieces,
                "types": set(),
                "partners": partner_tokens,
                "is_interaction": bool(interaction_partners),
                "residue_token": metadata["residue_token"],
                "residue_name": metadata["residue_name"],
            }

        if self.prolif_data is not None:
            for data in residue_data.values():
                residue_token = data["residue_token"]
                if include_types:
                    data["types"].update(self.prolif_data.get(residue_token, ()))
            return residue_data

        if not include_types or not self.atoms_B or self.tree_B is None:
            return residue_data

        geometric_types = self._geometric_interaction_types()
        for residue_label, data in residue_data.items():
            data["types"].update(geometric_types.get(residue_label, ()))
        return residue_data

    def count_patch_interaction_residues(self, patch):
        """Count authoritative interaction residues represented on one patch."""

        try:
            key = "uv_global" if "uv_global" in patch.metadata else "uv"
            uv = as_corner_uv(patch, key=key)
        except ValueError:
            return 0
        residue_data = self._collect_patch_residue_data(patch, uv, include_types=False)
        return sum(bool(data["is_interaction"]) for data in residue_data.values())

    def _relax_labels(self, ax, label_records, max_iter=80):
        fig = ax.figure
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()

        for _ in range(max_iter):
            moved = False
            bboxes = [rec["text"].get_window_extent(renderer=renderer).expanded(1.05, 1.1) for rec in label_records]
            for i in range(len(label_records)):
                for j in range(i + 1, len(label_records)):
                    b1, b2 = bboxes[i], bboxes[j]
                    if not b1.overlaps(b2):
                        continue
                    moved = True
                    dx = min(b1.x1 - b2.x0, b2.x1 - b1.x0) * 0.5
                    dy = min(b1.y1 - b2.y0, b2.y1 - b1.y0) * 0.5
                    sign_x = 1 if b1.x0 <= b2.x0 else -1
                    sign_y = 1 if b1.y0 <= b2.y0 else -1
                    self._move_text_display(ax, label_records[i]["text"], -sign_x * dx, -sign_y * dy)
                    self._move_text_display(ax, label_records[j]["text"], sign_x * dx, sign_y * dy)

            if not moved:
                break
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()

        for rec in label_records:
            pt = rec["scatter"].get_offsets()[0]
            tx, ty = rec["text"].get_position()
            if rec.get("connector") is not None:
                rec["connector"].set_data([pt[0], tx], [pt[1], ty])

    def _move_text_display(self, ax, txt, dx, dy):
        cur_data = txt.get_position()
        cur_disp = ax.transData.transform(cur_data)
        new_disp = (cur_disp[0] + dx, cur_disp[1] + dy)
        new_data = ax.transData.inverted().transform(new_disp)
        txt.set_position((new_data[0], new_data[1]))
