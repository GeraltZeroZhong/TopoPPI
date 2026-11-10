import tempfile
import unittest
from pathlib import Path

from topoppi.io.io_loader import PDBLoader
from topoppi.io.pdb_records import residue_plddt_values

FIXTURES = Path(__file__).parent / "fixtures"


class PDBLoaderTests(unittest.TestCase):
    def test_empty_structure_reports_coordinate_problem(self):
        with tempfile.TemporaryDirectory() as tmp:
            for name, text, message in [
                ("empty.pdb", "END\n", "no coordinate model"),
                ("empty.cif", "data_empty\n#\n", "required coordinate field"),
            ]:
                with self.subTest(name=name):
                    path = Path(tmp) / name
                    path.write_text(text)
                    with self.assertRaisesRegex(ValueError, message):
                        PDBLoader(path)

    def test_plddt_summary_gives_each_residue_one_vote(self):
        pdb_text = (
            "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 60.00           C  \n"
            "ATOM      2  N   ALA A   2       1.000   0.000   0.000  1.00100.00           N  \n"
            "ATOM      3  CA  ALA A   2       2.000   0.000   0.000  1.00100.00           C  \n"
            "ATOM      4  C   ALA A   2       3.000   0.000   0.000  1.00100.00           C  \n"
            "TER\nEND\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "prediction.pdb"
            path.write_text(pdb_text, encoding="ascii")
            _coordinates, atoms = PDBLoader(str(path)).get_chain_data("A")

            values = residue_plddt_values(atoms)

        self.assertEqual(values, [60.0, 100.0])
        self.assertEqual(sum(values) / len(values), 80.0)

    def test_loads_chains_and_filters_heteroatoms(self):
        loader = PDBLoader(str(FIXTURES / "tiny_complex.pdb"))

        self.assertEqual(loader.get_protein_chain_ids(), ["A", "B"])
        coords_a, atoms_a = loader.get_chain_data("A")
        coords_b, atoms_b = loader.get_chain_data("B")

        self.assertEqual(coords_a.shape, (6, 3))
        self.assertEqual(len(atoms_a), 6)
        self.assertEqual(coords_b.shape, (6, 3))
        self.assertEqual(len(atoms_b), 6)
        self.assertEqual(loader.get_chain_residue_count("A"), 2)

    def test_missing_chain_raises(self):
        loader = PDBLoader(str(FIXTURES / "tiny_complex.pdb"))
        with self.assertRaises(ValueError):
            loader.get_chain_data("Z")

    def test_synthetic_fixture_has_complete_heavy_atom_records(self):
        loader = PDBLoader(str(FIXTURES / "synthetic_complex.pdb"))

        self.assertEqual(loader.get_protein_chain_ids(), ["A", "B"])
        self.assertEqual(loader.get_chain_residue_count("A"), 2)
        coords_a, atoms_a = loader.get_chain_data("A")
        self.assertEqual(coords_a.shape, (17, 3))
        self.assertEqual(len(atoms_a), 17)

    def test_auto_contact_chain_selection_is_deterministic_and_audited(self):
        loader = PDBLoader(str(FIXTURES / "tiny_complex.pdb"))
        first = loader.select_contact_chain_pair(min_chain_residues=1, distance_cutoff=4.0)
        second = loader.select_contact_chain_pair(min_chain_residues=1, distance_cutoff=4.0)

        self.assertEqual(first, second)
        self.assertEqual(first[:2], ("A", "B"))
        self.assertGreater(first[2]["contact_residue_pair_count"], 0)
        self.assertGreater(first[2]["contact_atom_pair_count"], 0)

    def test_declared_partner_groups_preserve_receptor_ligand_orientation(self):
        loader = PDBLoader(str(FIXTURES / "tiny_complex.pdb"))

        chain_a, chain_b, details = loader.select_contact_chain_pair_between_groups(
            ["B"], ["A"], min_chain_residues=1, distance_cutoff=4.0
        )

        self.assertEqual((chain_a, chain_b), ("B", "A"))
        self.assertEqual(details["candidate_chain_pair_count"], 1)
        self.assertEqual(details["selected_residue_contact_fraction"], 1.0)

    def test_nucleic_acid_residues_are_not_reported_as_protein_chains(self):
        pdb_text = (
            "ATOM      1  P    DA A   1       0.000   0.000   0.000  1.00 20.00           P  \n"
            "ATOM      2  CA  ALA B   1       1.000   0.000   0.000  1.00 20.00           C  \n"
            "END\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "mixed.pdb"
            path.write_text(pdb_text, encoding="utf-8")
            loader = PDBLoader(str(path))

            self.assertEqual(loader.get_protein_chain_ids(), ["B"])
            self.assertEqual(loader.get_chain_residue_count("A"), 0)
            coords, atoms = loader.get_chain_data("A")
            self.assertEqual(coords.shape, (0, 3))
            self.assertEqual(atoms, [])

    def test_recognized_modified_amino_acid_is_kept_but_cofactor_is_not(self):
        pdb_text = (
            "HETATM    1  CA  MSE A   1       0.000   0.000   0.000  1.00 20.00           C  \n"
            "HETATM    2  SE  MSE A   1       1.000   0.000   0.000  1.00 20.00          SE  \n"
            "HETATM    3 ZN    ZN A   2       2.000   0.000   0.000  1.00 20.00          ZN  \n"
            "END\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "modified.pdb"
            path.write_text(pdb_text, encoding="utf-8")
            loader = PDBLoader(str(path))

            self.assertEqual(loader.get_protein_chain_ids(), ["A"])
            self.assertEqual(loader.get_chain_residue_count("A"), 1)
            coords, atoms = loader.get_chain_data("A")
            self.assertEqual(coords.shape, (2, 3))
            self.assertEqual({atom.get_parent().get_resname() for atom in atoms}, {"MSE"})
