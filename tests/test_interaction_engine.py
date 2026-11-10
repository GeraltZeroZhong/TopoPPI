import json
import logging
import tempfile
import unittest
import warnings
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
import pytest
from Bio.PDB import MMCIFIO, PDBParser

from topoppi.interactions.interaction_engine import (
    INTERACTION_SCHEMA_VERSION,
    generate_prolif_interactions,
    load_prolif_document,
    load_prolif_partner_map,
)
from topoppi.io.io_loader import PDBLoader
from topoppi_chemistry.backends.prolif import (
    _mda_to_prolif_with_explicit_hydrogen,
    _prolif_structure_input,
    _to_records,
)

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture(autouse=True)
def chemistry_provider(monkeypatch):
    """Source tests register the same provider and entry point as the wheel."""
    from importlib.metadata import EntryPoint

    from topoppi import plugins
    from topoppi_chemistry import api
    from topoppi_chemistry.plugin import plugin

    registry = plugins.PluginRegistry(entry_points=())
    registry.register_plugin(plugin)
    monkeypatch.setattr(plugins, "_default_registry", registry)
    monkeypatch.setattr(api, "_backend_entry_points", lambda: (
        EntryPoint(name="prolif", value="topoppi_chemistry.backends.prolif:detect", group=api.BACKEND_ENTRY_POINT_GROUP),
    ))


class InteractionEngineTests(unittest.TestCase):
    def test_single_chain_pair_generates_contacts_without_a_process_pool(self):
        import prolif as plf

        with (
            tempfile.TemporaryDirectory() as output_dir,
            mock.patch.object(
                plf.Fingerprint,
                "_run_iter_parallel",
                side_effect=AssertionError("A single chain pair must run in the calling process"),
            ),
        ):
            result = generate_prolif_interactions(
                FIXTURES / "synthetic_complex.pdb", "A", "B", output_dir=output_dir,
            )
            payload = json.loads(Path(result).read_text(encoding="utf-8"))
            self.assertTrue(payload["interactions"])
            self.assertEqual(payload["engine"], "prolif")

    def test_mmcif_suffixes_support_automatic_prolif_generation(self):
        with tempfile.TemporaryDirectory() as tmp:
            for suffix, chain_a, chain_b in ((".cif", "A", "B"), (".mmcif", "surface", "partner")):
                with self.subTest(suffix=suffix, chain_a=chain_a, chain_b=chain_b):
                    structure = PDBParser(QUIET=True).get_structure("complex", FIXTURES / "synthetic_complex.pdb")
                    if len(chain_a) > 1:
                        structure[0]["A"].id = chain_a
                        structure[0]["B"].id = chain_b
                    input_path = Path(tmp) / f"complex{suffix}"
                    output_dir = Path(tmp) / suffix.removeprefix(".")
                    output_dir.mkdir()
                    writer = MMCIFIO()
                    writer.set_structure(structure)
                    writer.save(str(input_path))

                    result = generate_prolif_interactions(
                        input_path,
                        chain_a,
                        chain_b,
                        output_dir=output_dir,
                    )

                    payload = json.loads(Path(result).read_text(encoding="utf-8"))
                    self.assertEqual(payload["chain_a"], chain_a)
                    self.assertEqual(payload["chain_b"], chain_b)
                    self.assertEqual(payload["engine"], "prolif")
                    self.assertTrue(payload["interactions"])

    def test_prolif_preparation_adds_hydrogens_without_mutating_source(self):
        import MDAnalysis as mda
        from rdkit import Chem

        universe = mda.Universe(FIXTURES / "synthetic_complex.pdb")
        source = universe.select_atoms("chainID A")
        source_coordinates = source.positions.copy()
        source_atom_count = source.n_atoms

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            import prolif as plf

            prepared = _mda_to_prolif_with_explicit_hydrogen(
                source,
                "A",
                plf,
                logging.getLogger(),
            )

        np.testing.assert_array_equal(source.positions, source_coordinates)
        self.assertEqual(source.n_atoms, source_atom_count)
        self.assertEqual(
            sum(atom.GetAtomicNum() != 1 for atom in prepared.GetAtoms()),
            source_atom_count,
        )
        hydrogens = [atom for atom in prepared.GetAtoms() if atom.GetAtomicNum() == 1]
        self.assertTrue(hydrogens)
        self.assertTrue(all(atom.GetPDBResidueInfo() is not None for atom in hydrogens))
        self.assertEqual(prepared.n_residues, len(source.residues))
        self.assertIn(
            Chem.BondType.DOUBLE,
            {bond.GetBondType() for bond in prepared.GetBonds()},
        )
        self.assertIn(
            Chem.BondType.AROMATIC,
            {bond.GetBondType() for bond in prepared.GetBonds()},
        )
        self.assertFalse(any("hydrogen atom" in str(item.message).lower() for item in caught))

    def test_selected_atom_preparation_invalidates_old_sidecars(self):
        self.assertEqual(INTERACTION_SCHEMA_VERSION, 4)

    def test_prolif_and_surface_use_identical_selected_conformers_for_pdb_and_mmcif(self):
        import MDAnalysis as mda
        import prolif as plf

        lines = (FIXTURES / "synthetic_complex.pdb").read_text().splitlines()
        index = next(i for i, line in enumerate(lines)
                     if line.startswith("ATOM") and line[21] == "A" and line[12:16].strip() == "CA")
        original = lines[index]
        x = float(original[30:38])
        a = original[:16] + "A" + original[17:54] + "  0.20" + original[60:]
        b = (original[:6] + "99999" + original[11:16] + "B" + original[17:30]
             + f"{x + .25:8.3f}" + original[38:54] + "  0.80" + original[60:])
        lines[index:index + 1] = [a, b]
        # A B-only atom must also survive RDKit's alternate-location handling.
        n_index = next(i for i, line in enumerate(lines)
                       if line.startswith("ATOM") and line[21] == "A" and line[12:16].strip() == "N")
        lines[n_index] = lines[n_index][:16] + "B" + lines[n_index][17:]
        with tempfile.TemporaryDirectory() as tmp:
            pdb = Path(tmp) / "alternate.pdb"
            pdb.write_text("\n".join(lines) + "\n")
            structure = PDBParser(QUIET=True).get_structure("input", pdb)
            structure[0]["A"].id, structure[0]["B"].id = "surface", "partner"
            cif = Path(tmp) / "alternate.mmcif"
            writer = MMCIFIO()
            writer.set_structure(structure)
            writer.save(str(cif))
            for path, chain_a, chain_b in ((pdb, "A", "B"), (cif, "surface", "partner")):
                with self.subTest(format=path.suffix):
                    original_bytes = path.read_bytes()
                    coords, atoms = PDBLoader(path).get_chain_data(chain_a)
                    expected = {(atom.get_name(), atom.get_parent().get_resname(), atom.get_parent().id[1],
                                 atom.get_parent().id[2].strip()): xyz for atom, xyz in zip(atoms, coords, strict=True)}
                    self.assertEqual(atoms[0].get_altloc(), "B")
                    with _prolif_structure_input(path, chain_a, chain_b) as (selected_path, selected_a, _selected_b):
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            group = mda.Universe(selected_path).select_atoms(f"chainID {selected_a}")
                            molecule = _mda_to_prolif_with_explicit_hydrogen(group, chain_a, plf, logging.getLogger())
                    observed = {}
                    for atom in molecule.GetAtoms():
                        if atom.GetAtomicNum() == 1:
                            continue
                        info = atom.GetPDBResidueInfo()
                        key = (info.GetName().strip(), info.GetResidueName().strip(),
                               info.GetResidueNumber(), info.GetInsertionCode().strip())
                        position = molecule.GetConformer().GetAtomPosition(atom.GetIdx())
                        observed[key] = [position.x, position.y, position.z]
                        self.assertFalse(info.GetAltLoc().strip())
                    self.assertEqual(set(observed), set(expected))
                    for key, xyz in expected.items():
                        np.testing.assert_allclose(observed[key], xyz, atol=1e-4, rtol=0)
                    self.assertEqual(path.read_bytes(), original_bytes)

    def test_previous_interaction_schema_is_regenerated_and_current_schema_is_reused(self):
        with tempfile.TemporaryDirectory() as tmp:
            sidecar = Path(tmp) / "synthetic_complex.A-B.prolif.json"
            sidecar.write_text(json.dumps({"source_sha256": "same-source", "interaction_schema_version": 3,
                                          "interactions": [{"obsolete": True}]}))
            output = generate_prolif_interactions(FIXTURES / "synthetic_complex.pdb", "A", "B",
                                                 output_dir=tmp, source_sha256="same-source")
            payload = json.loads(Path(output).read_text())
            self.assertEqual(payload["interaction_schema_version"], INTERACTION_SCHEMA_VERSION)
            self.assertTrue(payload["interactions"])
            self.assertNotIn("obsolete", payload["interactions"][0])
            with mock.patch("topoppi_chemistry.backends.prolif._prolif_structure_input",
                            side_effect=AssertionError("Current cache should be reused")):
                self.assertEqual(generate_prolif_interactions(FIXTURES / "synthetic_complex.pdb", "A", "B",
                                 output_dir=tmp, source_sha256="same-source"), output)

    def test_load_prolif_document_normalizes_records(self):
        _payload, residue_types, partners = load_prolif_document(str(FIXTURES / "prolif_interactions.json"))

        self.assertEqual(residue_types["1"], {"HydrogenBond"})
        self.assertEqual(residue_types["2A"], {"Hydrophobic"})
        self.assertEqual(partners["1"]["2"], 1)
        self.assertEqual(partners["2A"]["5"], 1)
        self.assertNotIn("bad", residue_types)

    def test_prolif_partner_map_resolves_authoritative_structure_labels(self):
        loader = PDBLoader(FIXTURES / "tiny_complex.pdb")
        _coords_a, atoms_a = loader.get_chain_data("A")
        _coords_b, atoms_b = loader.get_chain_data("B")

        partners = load_prolif_partner_map(
            FIXTURES / "prolif_interactions.json",
            atoms_a,
            atoms_b,
            expected_chain_a="A",
            expected_chain_b="B",
        )

        self.assertEqual(partners, {"A:GLY:1": {"B:ALA:2": 1}})

    def test_prolif_partner_map_rejects_a_declared_source_mismatch(self):
        loader = PDBLoader(FIXTURES / "tiny_complex.pdb")
        _coords_a, atoms_a = loader.get_chain_data("A")
        _coords_b, atoms_b = loader.get_chain_data("B")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "interactions.json"
            path.write_text(
                json.dumps(
                    {
                        "chain_a": "A",
                        "chain_b": "B",
                        "source_sha256": "a" * 64,
                        "interactions": [
                            {
                                "res_a_seq": "1",
                                "res_b_seq": "2",
                                "interaction": "Hydrophobic",
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "source checksum"):
                load_prolif_partner_map(
                    path,
                    atoms_a,
                    atoms_b,
                    expected_chain_a="A",
                    expected_chain_b="B",
                    expected_source_sha256="b" * 64,
                )

    def test_dataframe_records_use_all_three_column_levels(self):
        columns = pd.MultiIndex.from_tuples(
            [
                ("GLU5.B", "LYS10.A", "XBAcceptor"),
                ("PHE8.B", "TYR11.A", "FaceToFace"),
                ("ASP9.B", "ARG12.A", "Anionic"),
                ("GLY10.B", "ALA13.A", "FutureInteraction"),
            ],
            names=["ligand", "protein", "interaction"],
        )
        dataframe = pd.DataFrame([[True, True, True, True]], columns=columns)

        records = _to_records(dataframe)

        self.assertEqual(
            records,
            [
                {"res_a_seq": "10", "res_b_seq": "5", "interaction": "HalogenBond"},
                {"res_a_seq": "11", "res_b_seq": "8", "interaction": "PiStacking"},
                {"res_a_seq": "12", "res_b_seq": "9", "interaction": "Ionic"},
                {"res_a_seq": "13", "res_b_seq": "10", "interaction": "Other"},
            ],
        )

    def test_dataframe_records_recover_unique_insertion_codes(self):
        columns = pd.MultiIndex.from_tuples([("GLY5.B", "SER10.A", "HBDonor"), ("GLY6.B", "SER11.A", "HBAcceptor")])
        dataframe = pd.DataFrame([[True, True]], columns=columns)
        tokens_a = {
            ("SER", 10): {"10A"},
            (None, 10): {"10A"},
            ("SER", 11): {"11", "11A"},
            (None, 11): {"11", "11A"},
        }

        records = _to_records(dataframe, residue_tokens_a=tokens_a)

        self.assertEqual(records, [{"res_a_seq": "10A", "res_b_seq": "5", "interaction": "HydrogenBond"}])
