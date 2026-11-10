"""Failures here would invalidate supervised representations or mutation labels."""
import tempfile
import unittest
from pathlib import Path

import numpy as np

from topoppi_affinity.features import extract_complex_features, mutation_features, parse_mutations

FIXTURES = Path(__file__).parent / "fixtures"


class AffinityFeatureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bundle = extract_complex_features(FIXTURES / "tiny_complex.pdb", "A", "B")

    def test_complete_partner_swap_invariant_and_surface_is_physical(self):
        swapped = extract_complex_features(FIXTURES / "tiny_complex.pdb", "B", "A")
        self.assertEqual(self.bundle["absolute_features"], swapped["absolute_features"])
        area = self.bundle["absolute_features"]["surface.interface_area_sum"]
        self.assertGreater(area, 0.)
        self.assertAlmostEqual(area, sum(r["footprint_area"] for r in self.bundle["residues"].values()))

    def test_only_explicit_physical_features_no_annotation_leakage(self):
        mutated_metadata = dict(self.bundle)
        mutated_metadata["annotation_values"] = {"A:GLY:1": 99999.}
        mutated_metadata["ddg"] = 99999.
        expected = mutation_features(self.bundle, "GA1W")
        self.assertEqual(expected, mutation_features(mutated_metadata, "GA1W"))
        self.assertFalse(any("annotation" in k or "ddg" in k for k in expected))

    def test_wt_numbering_and_identity_validation(self):
        features = mutation_features(self.bundle, "GA1W")
        self.assertEqual(features["mutation.wt_G"], 1.)
        self.assertEqual(features["mutation.mt_W"], 1.)
        with self.assertRaisesRegex(ValueError, "WT identity mismatch"):
            mutation_features(self.bundle, "WA1F")
        with self.assertRaisesRegex(ValueError, "absent"):
            mutation_features(self.bundle, "GA100W")
        self.assertEqual(parse_mutations("WA23AF")[0]["insertion"], "A")
        self.assertEqual(parse_mutations("GA-1W")[0]["position"], -1)

    def test_structure_insertion_numbering_and_multichain_completeness(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "inserted.pdb"
            lines = (FIXTURES / "tiny_complex.pdb").read_text().splitlines()
            lines = [line[:22] + "  23A" + line[27:]
                     if line.startswith("ATOM") and line[21] == "A" and int(line[22:26]) == 1
                     else line for line in lines]
            path.write_text("\n".join(lines) + "\n")
            bundle = extract_complex_features(path, "A", "B")
            self.assertIn("A:GLY:23A", bundle["residues"])
            self.assertEqual(mutation_features(bundle, "GA23AW")["mutation.wt_G"], 1.)
            with self.assertRaisesRegex(ValueError, "absent"):
                mutation_features(bundle, "GA23W")
        grouped = extract_complex_features(FIXTURES / "1bvk.pdb", "A,B", "C")
        self.assertEqual({r["chain"] for r in grouped["residues"].values()}, {"A", "B", "C"})
        self.assertEqual(grouped["absolute_features"]["sequence.residue_count_sum"], 353.)

    def test_noop_zero_and_joint_features_nonadditive(self):
        noop = mutation_features(self.bundle, "GA1G")
        self.assertTrue(all(value == 0 for value in noop.values()))
        joint = mutation_features(self.bundle, "GA1W,AA2D")
        first = mutation_features(self.bundle, "GA1W")
        second = mutation_features(self.bundle, "AA2D")
        self.assertEqual(joint["mutation.delta_charge"], first["mutation.delta_charge"] + second["mutation.delta_charge"])
        self.assertEqual(joint["epistasis.pair_count"], 1.)
        self.assertNotEqual(joint["epistasis.delta_hydropathy_product_proximity"], 0.)
        self.assertTrue(np.isfinite(list(joint.values())).all())

    def test_distal_substitution_reported_not_dropped(self):
        bundle = {**self.bundle, "residues": {k: dict(v) for k, v in self.bundle["residues"].items()}}
        bundle["residues"]["A:GLY:1"]["partner_distance"] = 25.
        bundle["residues"]["A:GLY:1"]["footprint_area"] = 0.
        features = mutation_features(bundle, "GA1W")
        self.assertEqual(features["mutation.n_mutations"], 1.)
        self.assertEqual(features["mutation.contact_distal_count"], 1.)
        self.assertEqual(features["mutation.surface_covered_count"], 0.)

    def test_cache_same_features_and_complete_groups(self):
        with tempfile.TemporaryDirectory() as cache:
            first = extract_complex_features(FIXTURES / "tiny_complex.pdb", "A", "B", cache)
            second = extract_complex_features(FIXTURES / "tiny_complex.pdb", "B", "A", cache)
            self.assertEqual(first, second)
        with self.assertRaisesRegex(ValueError, "overlap"):
            extract_complex_features(FIXTURES / "tiny_complex.pdb", "A", "A,B")


if __name__ == "__main__":
    unittest.main()
