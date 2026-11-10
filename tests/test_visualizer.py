import json
import tempfile
import unittest
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import trimesh
from matplotlib.colors import to_rgba

from topoppi.atlas.uv import set_uv_layout
from topoppi.config import VisualizationConfig
from topoppi.interactions.metadata import INTERACTION_TYPES
from topoppi.io.io_loader import PDBLoader
from topoppi.visualization.visualizer import InterfaceVisualizer

FIXTURES = Path(__file__).parent / "fixtures"


class VisualizerTests(unittest.TestCase):
    @staticmethod
    def _two_residue_patch():
        mesh = trimesh.Trimesh(
            vertices=np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.asarray([[0, 1, 2], [0, 2, 3]]),
            process=False,
        )
        mesh.metadata["source_atom_indices"] = np.asarray([0, 1, 4, 5], dtype=np.int64)
        set_uv_layout(mesh, np.asarray(mesh.vertices[:, :2]))
        return mesh

    def test_residue_token_preserves_pdb_insertion_code(self):
        self.assertEqual(InterfaceVisualizer._residue_token((" ", 42, "A")), "42A")
        self.assertEqual(InterfaceVisualizer._residue_token((" ", -1, " ")), "-1")

    def test_specific_interactions_rank_before_generic_contacts(self):
        self.assertLess(INTERACTION_TYPES.index("HydrogenBond"), INTERACTION_TYPES.index("VdWContact"))
        self.assertLess(INTERACTION_TYPES.index("Ionic"), INTERACTION_TYPES.index("VdWContact"))

    def test_map_title_identifies_the_structure_and_chain_direction(self):
        visualizer = object.__new__(InterfaceVisualizer)
        visualizer.structure_label = "synthetic_complex"
        visualizer.chain_a_id = "A"
        visualizer.chain_b_id = "B"

        self.assertEqual(visualizer._map_title(), "synthetic_complex - surface A / partner B")

    def test_fragmented_residue_is_drawn_at_each_footprint_piece(self):
        loader = PDBLoader(FIXTURES / "tiny_complex.pdb")
        coords_a, atoms_a = loader.get_chain_data("A")
        coords_b, atoms_b = loader.get_chain_data("B")
        mesh = trimesh.Trimesh(
            vertices=np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.asarray([[0, 1, 2], [0, 2, 3]]),
            process=False,
        )
        mesh.metadata["source_atom_indices"] = np.arange(4, dtype=np.int64)
        corners = np.asarray(mesh.vertices[:, :2])[mesh.faces]
        corners[1] += np.asarray([3.0, 0.0])
        set_uv_layout(mesh, corners)
        visualizer = InterfaceVisualizer(
            chain_A_atoms=atoms_a,
            chain_A_coords=coords_a,
            chain_B_coords=coords_b,
            chain_B_atoms=atoms_b,
            chain_a_id="A",
            chain_b_id="B",
            config=VisualizationConfig(
                color_by_interaction_type=False,
                residue_scope="patch",
            ),
        )

        data = visualizer._collect_patch_residue_data(mesh, corners, include_types=False)
        self.assertEqual(len(data), 1)
        pieces = next(iter(data.values()))["pieces"]
        self.assertEqual(len(pieces), 2)

        figure = visualizer.plot_patches(
            [mesh],
            show=False,
            style_config={"use_uv_atlas": False, "avoid_label_overlap": False},
        )
        try:
            self.assertEqual(
                sorted(visualizer.artist_map),
                ["1_Gly1__piece_1", "1_Gly1__piece_2"],
            )
        finally:
            plt.close(figure)

        figure = visualizer.plot_patches(
            [mesh],
            show=False,
            style_config={
                "use_uv_atlas": False,
                "avoid_label_overlap": False,
                "marker_color_overrides": {"1_Gly1__piece_2": "#123456"},
            },
        )
        try:
            first = visualizer.artist_map["1_Gly1__piece_1"]["scatter"].get_facecolors()[0]
            second = visualizer.artist_map["1_Gly1__piece_2"]["scatter"].get_facecolors()[0]
            self.assertFalse(np.allclose(first, second))
            np.testing.assert_allclose(second, to_rgba("#123456"))
        finally:
            plt.close(figure)

    def test_default_scope_uses_prolif_interaction_residues(self):
        loader = PDBLoader(FIXTURES / "tiny_complex.pdb")
        coords_a, atoms_a = loader.get_chain_data("A")
        coords_b, atoms_b = loader.get_chain_data("B")
        mesh = self._two_residue_patch()
        visualizer = InterfaceVisualizer(
            chain_A_atoms=atoms_a,
            chain_A_coords=coords_a,
            chain_B_coords=coords_b,
            chain_B_atoms=atoms_b,
            chain_a_id="A",
            chain_b_id="B",
            prolif_file=str(FIXTURES / "prolif_interactions.json"),
            contact_distance_angstrom=6.0,
            config=VisualizationConfig(color_by_interaction_type=False),
        )

        figure = visualizer.plot_patches(
            [mesh],
            show=False,
            style_config={"use_uv_atlas": False, "avoid_label_overlap": False},
        )
        try:
            self.assertEqual(sorted(visualizer.artist_map), ["1_Gly1"])
            self.assertEqual(visualizer.last_report["patch_residue_count"], 2)
            self.assertEqual(visualizer.last_report["chain_interaction_residue_count"], 1)
            self.assertEqual(visualizer.last_report["patch_interaction_residue_count"], 1)
            self.assertEqual(
                visualizer.last_report["interaction_residue_retention_ratio"],
                1.0,
            )
            self.assertEqual(
                visualizer.last_report["interaction_residue_source"],
                "prolif",
            )
            self.assertEqual(visualizer.last_report["displayed_residue_count"], 1)
            self.assertEqual(visualizer.count_patch_interaction_residues(mesh), 1)
        finally:
            plt.close(figure)

        figure = visualizer.plot_patches(
            [mesh],
            show=False,
            style_config={
                "use_uv_atlas": False,
                "avoid_label_overlap": False,
                "residue_scope": "patch",
            },
        )
        try:
            self.assertEqual(sorted(visualizer.artist_map), ["1_Ala2", "1_Gly1"])
            self.assertEqual(visualizer.last_report["residue_scope"], "patch")
            self.assertEqual(visualizer.last_report["displayed_residue_count"], 2)
        finally:
            plt.close(figure)

    def test_type_coloring_keeps_untyped_residues_in_patch_scope(self):
        loader = PDBLoader(FIXTURES / "tiny_complex.pdb")
        coords_a, atoms_a = loader.get_chain_data("A")
        coords_b, atoms_b = loader.get_chain_data("B")
        visualizer = InterfaceVisualizer(
            chain_A_atoms=atoms_a,
            chain_A_coords=coords_a,
            chain_B_coords=coords_b,
            chain_B_atoms=atoms_b,
            chain_a_id="A",
            chain_b_id="B",
            prolif_file=str(FIXTURES / "prolif_interactions.json"),
        )

        figure = visualizer.plot_patches(
            [self._two_residue_patch()],
            show=False,
            style_config={
                "use_uv_atlas": False,
                "avoid_label_overlap": False,
                "residue_scope": "patch",
            },
        )
        try:
            self.assertEqual(sorted(visualizer.artist_map), ["1_Ala2", "1_Gly1"])
        finally:
            plt.close(figure)

    def test_geometric_residue_without_a_heuristic_type_remains_visible(self):
        loader = PDBLoader(FIXTURES / "tiny_complex.pdb")
        coords_a, atoms_a = loader.get_chain_data("A")
        _coords_b, atoms_b = loader.get_chain_data("B")
        coords_b = np.tile(np.asarray([0.0, 0.0, 5.5]), (len(atoms_b), 1))
        visualizer = InterfaceVisualizer(
            chain_A_atoms=atoms_a,
            chain_A_coords=coords_a,
            chain_B_coords=coords_b,
            chain_B_atoms=atoms_b,
            chain_a_id="A",
            chain_b_id="B",
            contact_distance_angstrom=5.6,
            config=VisualizationConfig(use_geometric_interaction_fallback=True),
        )

        figure = visualizer.plot_patches(
            [self._two_residue_patch()],
            show=False,
            style_config={"use_uv_atlas": False, "avoid_label_overlap": False},
        )
        try:
            self.assertEqual(sorted(visualizer.artist_map), ["1_Gly1"])
        finally:
            plt.close(figure)

    def test_saved_label_offset_is_relative_to_the_marker(self):
        loader = PDBLoader(FIXTURES / "tiny_complex.pdb")
        coords_a, atoms_a = loader.get_chain_data("A")
        coords_b, atoms_b = loader.get_chain_data("B")
        visualizer = InterfaceVisualizer(
            chain_A_atoms=atoms_a,
            chain_A_coords=coords_a,
            chain_B_coords=coords_b,
            chain_B_atoms=atoms_b,
            config=VisualizationConfig(
                color_by_interaction_type=False,
                residue_scope="patch",
            ),
        )

        figure = visualizer.plot_patches(
            [self._two_residue_patch()],
            show=False,
            style_config={
                "use_uv_atlas": False,
                "avoid_label_overlap": False,
                "label_offsets": {"1_Gly1": (0.2, 0.3)},
            },
        )
        try:
            artist = visualizer.artist_map["1_Gly1"]
            marker = artist["scatter"].get_offsets()[0]
            text = artist["text"].get_position()
            np.testing.assert_allclose(np.asarray(text) - marker, [0.2, 0.3])
        finally:
            plt.close(figure)

    def test_prolif_types_exclude_unresolved_partner_records(self):
        loader = PDBLoader(FIXTURES / "tiny_complex.pdb")
        coords_a, atoms_a = loader.get_chain_data("A")
        coords_b, atoms_b = loader.get_chain_data("B")
        with tempfile.TemporaryDirectory() as tmp:
            interactions = Path(tmp) / "interactions.json"
            interactions.write_text(
                json.dumps(
                    {
                        "engine": "prolif",
                        "chain_a": "A",
                        "chain_b": "B",
                        "interactions": [
                            {"res_a_seq": "1", "res_b_seq": "1", "interaction": "VdWContact"},
                            {"res_a_seq": "1", "res_b_seq": "999", "interaction": "HBAcceptor"},
                        ],
                    }
                ),
                encoding="utf-8",
            )
            visualizer = InterfaceVisualizer(
                chain_A_atoms=atoms_a,
                chain_A_coords=coords_a,
                chain_B_coords=coords_b,
                chain_B_atoms=atoms_b,
                chain_a_id="A",
                chain_b_id="B",
                prolif_file=str(interactions),
            )

            data = visualizer._collect_patch_residue_data(
                self._two_residue_patch(),
                np.asarray(self._two_residue_patch().vertices[:, :2]),
            )

        self.assertEqual(data["A:GLY:1"]["types"], {"VdWContact"})
        self.assertEqual(data["A:GLY:1"]["partners"], {"1": 1})

    def test_geometric_interaction_types_are_cached_across_redraws(self):
        loader = PDBLoader(FIXTURES / "tiny_complex.pdb")
        coords_a, atoms_a = loader.get_chain_data("A")
        coords_b, atoms_b = loader.get_chain_data("B")
        mesh = self._two_residue_patch()
        visualizer = InterfaceVisualizer(
            chain_A_atoms=atoms_a,
            chain_A_coords=coords_a,
            chain_B_coords=coords_b,
            chain_B_atoms=atoms_b,
            chain_a_id="A",
            chain_b_id="B",
            config=VisualizationConfig(use_geometric_interaction_fallback=True),
        )
        style = {
            "use_uv_atlas": False,
            "avoid_label_overlap": False,
            "color_by_type": True,
        }

        first_figure = visualizer.plot_patches([mesh], show=False, style_config=style)
        first_cache = visualizer._geometric_types_cache
        second_figure = visualizer.plot_patches([mesh], show=False, style_config=style)
        try:
            self.assertTrue(visualizer.artist_map)
            self.assertIs(visualizer._geometric_types_cache, first_cache)
        finally:
            plt.close(first_figure)
            plt.close(second_figure)

    def test_default_type_rendering_uses_prolif_source(self):
        loader = PDBLoader(FIXTURES / "tiny_complex.pdb")
        coords_a, atoms_a = loader.get_chain_data("A")
        coords_b, atoms_b = loader.get_chain_data("B")
        mesh = self._two_residue_patch()
        visualizer = InterfaceVisualizer(
            chain_A_atoms=atoms_a,
            chain_A_coords=coords_a,
            chain_B_coords=coords_b,
            chain_B_atoms=atoms_b,
            chain_a_id="A",
            chain_b_id="B",
            prolif_file=str(FIXTURES / "prolif_interactions.json"),
        )

        figure = visualizer.plot_patches(
            [mesh],
            show=False,
            style_config={"use_uv_atlas": False, "avoid_label_overlap": False},
        )
        try:
            self.assertEqual(sorted(visualizer.artist_map), ["1_Gly1"])
            self.assertTrue(visualizer.last_report["color_by_interaction_type"])
            self.assertEqual(
                visualizer.last_report["interaction_type_source"],
                "prolif",
            )
            self.assertEqual(len(figure.legends), 1)
            self.assertFalse(figure.legends[0].get_in_layout())
            self.assertLessEqual(figure.axes[0].get_position().y1, 0.87)
        finally:
            plt.close(figure)

    def test_type_rendering_does_not_silently_infer_geometry(self):
        loader = PDBLoader(FIXTURES / "tiny_complex.pdb")
        coords_a, atoms_a = loader.get_chain_data("A")
        coords_b, atoms_b = loader.get_chain_data("B")
        visualizer = InterfaceVisualizer(
            chain_A_atoms=atoms_a,
            chain_A_coords=coords_a,
            chain_B_coords=coords_b,
            chain_B_atoms=atoms_b,
            chain_a_id="A",
            chain_b_id="B",
        )

        with self.assertRaisesRegex(ValueError, "requires an interaction JSON"):
            visualizer.plot_patches(
                [self._two_residue_patch()],
                show=False,
                style_config={"use_uv_atlas": False},
            )
