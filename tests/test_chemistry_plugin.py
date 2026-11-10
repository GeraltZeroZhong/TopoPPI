"""The default chemistry provider consumes Core objects and preserves records."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from importlib.metadata import EntryPoint
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from topoppi import compute_geometry, get_interface, load_structure, plugins
from topoppi.execution import ExecutionContext
from topoppi_chemistry import api
from topoppi_chemistry.plugin import plugin

FIXTURE = Path(__file__).parent / "fixtures/synthetic_complex.pdb"
EXPECTED_RECORDS = [
    {"res_a_seq": "1", "res_b_seq": "1", "interaction": "VdWContact"},
    {"res_a_seq": "1", "res_b_seq": "2", "interaction": "HydrogenBond"},
    {"res_a_seq": "1", "res_b_seq": "2", "interaction": "VdWContact"},
    {"res_a_seq": "2", "res_b_seq": "2", "interaction": "Hydrophobic"},
    {"res_a_seq": "2", "res_b_seq": "2", "interaction": "PiStacking"},
    {"res_a_seq": "2", "res_b_seq": "2", "interaction": "VdWContact"},
]


@pytest.fixture(autouse=True)
def chemistry_provider(monkeypatch):
    registry = plugins.PluginRegistry(entry_points=())
    registry.register_plugin(plugin)
    monkeypatch.setattr(plugins, "_default_registry", registry)
    monkeypatch.setattr(api, "_backend_entry_points", lambda: (
        EntryPoint(name="prolif", value="topoppi_chemistry.backends.prolif:detect", group=api.BACKEND_ENTRY_POINT_GROUP),
    ))


@pytest.mark.parametrize("with_geometry", [False, True])
def test_detection_reuses_core_atoms_after_source_removal(tmp_path, monkeypatch, with_geometry):
    path = tmp_path / "synthetic.pdb"
    shutil.copyfile(FIXTURE, path)
    interface = get_interface(load_structure(path), "A", "B")
    selected = compute_geometry(interface, sides=()) if with_geometry else interface
    original_coordinates = [partner.coords.copy() for partner in selected.partners]
    original_atoms = [[atom.get_coord().copy() for atom in partner.atoms] for partner in selected.partners]
    path.unlink()

    def unexpected(*args, **kwargs):
        raise AssertionError("An existing Core object must not reparse the original structure")

    monkeypatch.setattr("topoppi.io.io_loader.PDBLoader", unexpected)
    events = []
    result = plugins.run("chemistry.detect", selected, context=ExecutionContext(progress=events.append))
    assert result.engine == "prolif"
    assert (result.chain_a, result.chain_b) == ("A", "B")
    assert result.interactions == EXPECTED_RECORDS
    assert result.metadata["analysis_id"] == selected.analysis_id
    assert result.to_document()["interactions"] == EXPECTED_RECORDS
    assert [event["stage"] for event in events] == [
        "started", "chemistry.preparing", "chemistry.detecting", "completed"]
    for partner, coordinates, atoms in zip(selected.partners, original_coordinates, original_atoms, strict=True):
        np.testing.assert_array_equal(partner.coords, coordinates)
        np.testing.assert_array_equal([atom.get_coord() for atom in partner.atoms], atoms)


def test_backend_discovery_does_not_load_unselected_backends(monkeypatch):
    from topoppi.errors import PluginDependencyError
    from topoppi.interactions.records import InteractionResult

    def broken():
        raise ImportError("unselected backend dependency")

    selected = SimpleNamespace(name="toy", load=lambda: lambda interface, context: InteractionResult(
        "toy", "A", "B", [], {"analysis_id": interface.analysis_id}))
    unavailable = SimpleNamespace(name="broken", load=broken)
    monkeypatch.setattr(api, "_backend_entry_points", lambda: (unavailable, selected))
    assert plugins.run("chemistry.backends") == ("broken", "toy")
    interface = get_interface(load_structure(FIXTURE), "A", "B")
    result = plugins.run("chemistry.detect", interface, backend="toy")
    assert result.engine == "toy" and result.metadata["analysis_id"] == interface.analysis_id
    with pytest.raises(ValueError, match="not installed"):
        api.detect(interface, backend="absent")
    with pytest.raises(PluginDependencyError, match="unavailable dependency"):
        plugins.run("chemistry.detect", interface, backend="broken")
    assert plugins.run("chemistry.detect", interface, backend="toy").engine == "toy"


def test_prolif_rejects_ambiguous_grouped_chain_schema():
    from topoppi_chemistry.backends.prolif import detect

    interface = get_interface(load_structure(FIXTURE), "A", "B")
    interface.partners[0].chains = ("A", "C")
    with pytest.raises(ValueError, match="one author chain per partner"):
        detect(interface)


@pytest.mark.parametrize("format", ["pdb", "cif"])
def test_negative_author_residue_numbers_survive_chemical_perception(tmp_path, format):
    from Bio.PDB import MMCIFIO, PDBIO, PDBParser

    structure = PDBParser(QUIET=True).get_structure("synthetic", FIXTURE)
    chain_names = ("AA", "20") if format == "cif" else ("A", "B")
    for chain, name in zip(structure[0], chain_names, strict=True):
        chain.id = name
        for residue in chain:
            residue.id = (" ", residue.id[1] - 3, "B" if residue.id[1] == 2 else " ")
    path = tmp_path / f"negative.{format}"
    writer = MMCIFIO() if format == "cif" else PDBIO()
    writer.set_structure(structure)
    writer.save(str(path))
    interface = get_interface(load_structure(path), *chain_names)
    coordinates = [partner.coords.copy() for partner in interface.partners]
    residue_ids = [[atom.parent.id for atom in partner.atoms] for partner in interface.partners]

    result = api.detect(interface)

    assert (result.chain_a, result.chain_b) == chain_names
    tokens = {"1": "-2", "2": "-1B"}
    expected = [{**record, "res_a_seq": tokens[record["res_a_seq"]],
                 "res_b_seq": tokens[record["res_b_seq"]]} for record in EXPECTED_RECORDS]
    assert result.interactions == expected
    for partner, xyz, identities in zip(interface.partners, coordinates, residue_ids, strict=True):
        np.testing.assert_array_equal(partner.coords, xyz)
        assert [atom.parent.id for atom in partner.atoms] == identities


@pytest.mark.parametrize("missing", ["structure", "backend"])
def test_direct_chemistry_module_reports_recoverable_errors(tmp_path, missing):
    args = [str(tmp_path / "missing.pdb") if missing == "structure" else str(FIXTURE),
            "--chain-a", "A", "--chain-b", "B"]
    if missing == "backend":
        args.extend(["--backend", "missing"])
    result = subprocess.run([sys.executable, "-m", "topoppi_chemistry.cli", "detect", *args],
                            capture_output=True, text=True, env=dict(os.environ))
    assert result.returncode == 1
    assert result.stdout == ""
    assert "Traceback" not in result.stderr
    assert "missing" in result.stderr


def test_metadata_discovery_and_help_do_not_import_scientific_dependencies():
    source = Path(__file__).resolve().parents[1] / "src"
    script = '''
import importlib.abc, sys
from importlib.metadata import EntryPoint
class DenyScientificImports(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'numpy','pandas','sklearn','torch','MDAnalysis','prolif','rdkit','tkinter'}:
            raise AssertionError('Heavy dependency imported: ' + fullname)
sys.meta_path.insert(0, DenyScientificImports())
from topoppi.plugins import PluginRegistry
from topoppi_chemistry.plugin import plugin
from topoppi_chemistry import api
api._backend_entry_points = lambda: (EntryPoint(name='prolif', value='topoppi_chemistry.backends.prolif:detect', group=api.BACKEND_ENTRY_POINT_GROUP),)
registry=PluginRegistry(entry_points=())
registry.register_plugin(plugin)
assert registry.info('chemistry').name == 'chemistry'
assert registry.run('chemistry.backends') == ('prolif',)
assert registry.run_command('chemistry', ['--help']) == 0
'''
    result = subprocess.run([sys.executable, "-I", "-c", f"import sys; sys.path.insert(0, {str(source)!r})\n" + script],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "backends" in result.stdout
