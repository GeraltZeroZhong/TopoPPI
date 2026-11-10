import io
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest import mock

import matplotlib.pyplot as plt
import numpy as np
import trimesh

from topoppi import cli
from topoppi.atlas.footprints import mesh_vertex_residue_labels
from topoppi.atlas.uv import as_corner_uv, set_uv_layout
from topoppi.config import VisualizationConfig
from topoppi.io.io_loader import PDBLoader
from topoppi.visualization.atlas_io import load_atlas, save_atlas
from topoppi.visualization.visualizer import InterfaceVisualizer, select_patches_for_display

FIXTURES = Path(__file__).parent / "fixtures"


def make_atlas(*, insertion_code=False, style="footprints"):
    loader = PDBLoader(FIXTURES / "tiny_complex.pdb")
    coords_a, atoms_a = loader.get_chain_data("A")
    coords_b, atoms_b = loader.get_chain_data("B")
    if insertion_code:
        atoms_a[0].get_parent().id = (" ", 1, "A")
    mesh = trimesh.Trimesh(
        vertices=np.asarray([[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]]),
        faces=np.asarray([[0, 1, 2], [0, 2, 3]]), process=False,
    )
    mesh.metadata["source_atom_indices"] = np.asarray([0, 1, 4, 5], dtype=np.int64)
    mesh.metadata["source_vertex_ids"] = np.arange(4)
    mesh.metadata["source_face_ids"] = np.asarray([11, 19])
    corners = np.asarray(mesh.vertices[:, :2])[mesh.faces].copy()
    corners[1] += [3., 0.]
    for key in ("uv", "uv_optcuts", "uv_global"):
        set_uv_layout(mesh, corners, key=key)
    viz = InterfaceVisualizer(
        atoms_a, coords_a, coords_b, atoms_b, chain_a_id="A", chain_b_id="B", structure_label="test interface",
        config=VisualizationConfig(map_style=style, color_by_interaction_type=False,
                                   residue_scope="patch", use_geometric_interaction_fallback=True),
    )
    return [mesh], viz


class AtlasIOTests(unittest.TestCase):
    def test_roundtrip_preserves_all_corner_coordinates_provenance_and_author_ids(self):
        patches, viz = make_atlas(insertion_code=True)
        labels = mesh_vertex_residue_labels(patches[0], viz.source_residue_labels_A)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "saved.npz"
            save_atlas(path, patches, viz, run_metadata={"objective": {"weight": 20.0}})
            restored = load_atlas(path)
            with np.load(path, allow_pickle=False) as saved:
                self.assertTrue(all(saved[key].dtype.kind != "O" for key in saved.files))
        for key in ("uv", "uv_optcuts", "uv_global"):
            np.testing.assert_array_equal(as_corner_uv(patches[0], key=key), as_corner_uv(restored.patches[0], key=key))
        for key in ("source_atom_indices", "source_vertex_ids", "source_face_ids"):
            np.testing.assert_array_equal(patches[0].metadata[key], restored.patches[0].metadata[key])
        np.testing.assert_array_equal(labels, mesh_vertex_residue_labels(restored.patches[0],
                                                                       restored.visualizer.source_residue_labels_A))
        self.assertIn("A:GLY:1A", restored.visualizer.residue_metadata_A)
        self.assertEqual(restored.metadata["objective"]["weight"], 20.)

    def test_saved_annotation_values_survive_missing_csv_and_recolour(self):
        patches, viz = make_atlas()
        with tempfile.TemporaryDirectory() as tmp:
            csv = Path(tmp) / "effects.csv"
            csv.write_text("residue,value\nA:GLY:1,1.25\nA:ALA:2,NA\n")
            figure = viz.plot_patches(patches, show=False, style_config={"annotation_file": str(csv)})
            plt.close(figure)
            path = Path(tmp) / "atlas.npz"
            save_atlas(path, patches, viz)
            csv.unlink()
            with mock.patch("topoppi.optimization.optcuts.OptCutsUVOptimizer", side_effect=AssertionError("No solver")):
                restored = load_atlas(path)
                self.assertEqual(restored.style["annotation_values"]["A:GLY:1"], 1.25)
                self.assertFalse(restored.style["annotation_file"])
                figure = restored.visualizer.plot_patches(restored.patches, show=False, style_config=restored.style)
                plt.close(figure)

    def test_prolif_membership_and_marker_style_survive_offline_load(self):
        patches, _viz = make_atlas(style="markers")
        loader = PDBLoader(FIXTURES / "tiny_complex.pdb")
        coords_a, atoms_a = loader.get_chain_data("A")
        coords_b, atoms_b = loader.get_chain_data("B")
        viz = InterfaceVisualizer(atoms_a, coords_a, coords_b, atoms_b, chain_a_id="A", chain_b_id="B",
                                  prolif_file=str(FIXTURES / "prolif_interactions.json"))
        original = viz.plot_patches(patches, show=False)
        expected = sorted(viz.artist_map)
        plt.close(original)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "atlas.custom"
            save_atlas(path, patches, viz)
            self.assertFalse(Path(str(path) + ".npz").exists())
            restored = load_atlas(path)
        self.assertEqual(restored.visualizer.interaction_partner_map, viz.interaction_partner_map)
        self.assertEqual(restored.visualizer.prolif_data, viz.prolif_data)
        figure = restored.visualizer.plot_patches(restored.patches, show=False, style_config=restored.style)
        self.assertEqual(sorted(restored.visualizer.artist_map), expected)
        plt.close(figure)

    def test_rejects_a_scientific_npz_without_atlas_document(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "other.npz"
            np.savez(path, uv=np.zeros((1, 3, 2)))
            with self.assertRaisesRegex(ValueError, "not a TopoPPI atlas"):
                load_atlas(path)

    def test_invalid_and_incomplete_atlas_files_explain_how_to_recover(self):
        array_file = io.BytesIO()
        np.save(array_file, np.zeros((2, 3)))
        archive = io.BytesIO()
        np.savez_compressed(archive, document_json=np.asarray('{}'))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "broken.npz"
            for payload in (b"not an atlas", b"", array_file.getvalue(), archive.getvalue()[:100]):
                with self.subTest(payload=payload[:12]):
                    path.write_bytes(payload)
                    with self.assertRaises(ValueError) as result:
                        load_atlas(path)
                    message = str(result.exception)
                    self.assertIn("Save Atlas or --export-atlas", message)
                    self.assertIn("export it again or obtain a complete copy", message)
                    self.assertNotIn("pickle", message)

    def test_valid_archive_keeps_the_specific_unsupported_format_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "unsupported.npz"
            np.savez(path, document_json=np.asarray('{"format":"topoppi-atlas","version":99}'))
            with self.assertRaisesRegex(ValueError, "Unsupported TopoPPI atlas format"):
                load_atlas(path)

    def test_cli_render_reuses_uv_and_saves_updated_highlight_style(self):
        patches, viz = make_atlas()
        with tempfile.TemporaryDirectory() as tmp:
            source, target = Path(tmp) / "atlas.npz", Path(tmp) / "styled.npz"
            output = Path(tmp) / "map.pdf"
            save_atlas(source, patches, viz)
            with mock.patch("topoppi.cli.run_interface_mapping", side_effect=AssertionError("No recomputation")):
                result = cli.main(["render", str(source), "-o", str(output), "--map-style", "footprints",
                                   "--highlight", "A:GLY:1", "--export-atlas", str(target)])
            self.assertEqual(result, 0)
            self.assertTrue(output.read_bytes().startswith(b"%PDF"))
            document = load_atlas(target)
            self.assertIn("A:GLY:1", document.style["highlight_residues"])
            np.testing.assert_array_equal(as_corner_uv(patches[0]), as_corner_uv(document.patches[0]))


    def test_cli_marker_roundtrip_preserves_filter_and_keeps_hidden_geometry(self):
        patches, viz = make_atlas(style="markers")
        viz.config = replace(viz.config, min_points=1)
        viz.interaction_partner_map = {"A:GLY:1": {"B:GLY:1": 1}}
        hidden = patches[0].copy()
        hidden.metadata["source_atom_indices"] = np.asarray([4, 5, 4, 5])
        for key in ("uv", "uv_optcuts", "uv_global"):
            set_uv_layout(hidden, as_corner_uv(hidden, key=key) + [8., 0.], key=key)
        all_patches = [patches[0], hidden]
        selected, counts = select_patches_for_display(all_patches, viz)
        self.assertEqual(counts, [1, 0])
        original = viz.plot_patches(selected, show=False)
        self.assertEqual(viz.last_report["patch_count"], 1)
        plt.close(original)
        reports = []
        plot = InterfaceVisualizer.plot_patches

        def capture(instance, *args, **kwargs):
            figure = plot(instance, *args, **kwargs)
            reports.append(dict(instance.last_report))
            return figure

        with tempfile.TemporaryDirectory() as tmp:
            source, saved = Path(tmp) / "source.npz", Path(tmp) / "saved.npz"
            output = Path(tmp) / "map.svg"
            save_atlas(source, all_patches, viz)
            with mock.patch.object(InterfaceVisualizer, "plot_patches", capture):
                result = cli.main(["render", str(source), "-o", str(output), "--export-atlas", str(saved)])
                self.assertEqual(result, 0)
                self.assertEqual(reports[-1]["patch_count"], 1)
                restored = load_atlas(saved)
                self.assertEqual(len(restored.patches), 2)
                np.testing.assert_array_equal(as_corner_uv(hidden), as_corner_uv(restored.patches[1]))
                result = cli.main(["render", str(saved), "--map-style", "footprints", "-o", str(output)])
                self.assertEqual(result, 0)
                self.assertEqual(reports[-1]["patch_count"], 2)


if __name__ == "__main__":
    unittest.main()
