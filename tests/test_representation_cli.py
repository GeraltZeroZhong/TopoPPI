"""Structure-to-table descriptor workflows use independent synthetic inputs."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from topoppi import cli, get_interface, load_structure
from topoppi_representation import physical_descriptors

FIXTURE = Path(__file__).parent / "fixtures/tiny_complex.pdb"


@pytest.mark.parametrize("format", ["pdb", "cif"])
def test_structure_cli_exports_matching_named_features_and_residue_tables(format, tmp_path, capsys):
    structure = FIXTURE
    if format == "cif":
        from Bio.PDB import MMCIFIO, PDBParser

        structure = tmp_path / "input.cif"
        writer = MMCIFIO()
        writer.set_structure(PDBParser(QUIET=True).get_structure("synthetic", FIXTURE))
        writer.save(str(structure))
    output = tmp_path / "descriptors"
    assert cli.main(["representation", "descriptors", "--structure", str(structure),
                     "--chains-a", "B", "--chains-b", "A", "--output-dir", str(output),
                     "--residues", "--contacts", "--n-points", "100", "--json"]) == 0
    result = json.loads(capsys.readouterr().out)
    expected = physical_descriptors(get_interface(load_structure(FIXTURE), "A", "B"), n_points=100)
    assert result["descriptors"] == expected
    exported = pd.read_csv(result["files"]["descriptors"])
    assert exported.shape == (1, 477)
    assert exported.columns.tolist() == list(expected)
    np.testing.assert_allclose(exported.iloc[0], list(expected.values()), atol=1e-12)
    residues = pd.read_csv(result["files"]["residues"])
    contacts = pd.read_csv(result["files"]["contacts"])
    assert len(residues) == result["n_residues"] == 4
    assert len(contacts) == result["n_contact_pairs"] == expected["physical.contact_pairs.total"]
    assert set(contacts.residue_a).issubset(set(residues.residue_key))
    assert residues.loc[residues.side.eq(0), "residue_key"].str.startswith("A:").all()
    metadata = json.loads(Path(result["files"]["metadata"]).read_text())
    assert metadata["partner_groups"] == [["A"], ["B"]]
    assert metadata["parameters"] == {"sasa_points": 100, "probe_radius_angstrom": 1.4,
                                      "contact_distance_angstrom": 5.5}


def test_default_export_and_module_entry_point(tmp_path):
    result = subprocess.run(
        [sys.executable, "-m", "topoppi_representation", "descriptors", "--structure", str(FIXTURE),
         "--chains-a", "A", "--chains-b", "B", "--output", str(tmp_path), "--json"],
        text=True, capture_output=True, env=dict(os.environ), check=True,
    )
    record = json.loads(result.stdout)
    assert record["parameters"]["sasa_points"] == 960
    assert record["n_descriptors"] == 477
    assert set(record["files"]) == {"descriptors", "metadata"}
    assert not (tmp_path / "residues.csv").exists()
    assert not (tmp_path / "contacts.csv").exists()


def test_representation_help_keeps_optional_dependencies_unloaded():
    program = """
import sys
from topoppi_representation.cli import main
try:
    main(['descriptors', '--help'])
except SystemExit as error:
    assert error.code == 0
assert not {'pandas', 'torch', 'sklearn', 'DracoPy', 'zstandard', 'pywt'} & sys.modules.keys()
"""
    result = subprocess.run([sys.executable, "-c", program], text=True, capture_output=True,
                            env=dict(os.environ), check=True)
    assert "usage: topoppi representation descriptors" in result.stdout


def test_missing_descriptor_dependency_has_install_instruction(tmp_path):
    program = """
import sys
sys.modules['pandas'] = None
from topoppi_representation.cli import main
raise SystemExit(main(sys.argv[1:]))
"""
    result = subprocess.run(
        [sys.executable, "-c", program, "descriptors", "--structure", str(FIXTURE),
         "--chains-a", "A", "--chains-b", "B", "--output", str(tmp_path / "results")],
        text=True, capture_output=True, env=dict(os.environ),
    )
    assert result.returncode == 1
    assert "topoppi[representation]" in result.stderr
    assert "pandas" in result.stderr
    assert "Traceback" not in result.stderr
    assert not (tmp_path / "results").exists()


def test_invalid_partner_can_be_corrected_without_partial_outputs(tmp_path, capsys):
    args = ["representation", "descriptors", "--structure", str(FIXTURE),
            "--chains-a", "A", "--chains-b", "missing", "--output", str(tmp_path / "results")]
    assert cli.main(args) == 1
    assert "Available protein chains" in capsys.readouterr().err
    assert not (tmp_path / "results").exists()
    args[args.index("missing")] = "B"
    assert cli.main(args) == 0
    assert (tmp_path / "results/descriptors.csv").is_file()


@pytest.mark.parametrize("option,value", [
    ("--n-points", "0"), ("--n-points", "-1"),
    ("--probe-radius", "nan"), ("--probe-radius", "inf"),
    ("--probe-radius", "0"), ("--probe-radius", "-1"),
])
def test_invalid_sasa_options_stop_before_loading_or_writing(option, value, tmp_path, capsys):
    from topoppi_representation.cli import main

    output = tmp_path / "results"
    with pytest.raises(SystemExit) as error:
        main(["descriptors", "--structure", str(tmp_path / "not_loaded.pdb"),
              "-A", "A", "-B", "B", "-o", str(output), option, value, "--json"])
    assert error.value.code == 2
    captured = capsys.readouterr()
    assert captured.out == ""
    assert option in captured.err and "positive" in captured.err
    assert "not_loaded.pdb" not in captured.err.splitlines()[-1]
    assert not output.exists()


def test_repeated_export_preserves_complete_results_and_can_use_new_directory(tmp_path, capsys):
    from topoppi_representation.cli import main

    output = tmp_path / "results"
    args = ["descriptors", "--structure", str(FIXTURE), "-A", "A", "-B", "B", "-o", str(output)]
    assert main([*args, "--residues", "--contacts", "--n-points", "100", "--json"]) == 0
    capsys.readouterr()
    before = {p.name: p.read_bytes() for p in output.iterdir()}
    assert main([*args, "--json"]) == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "Choose a new output directory" in captured.err
    assert {p.name: p.read_bytes() for p in output.iterdir()} == before
    fresh = tmp_path / "fresh"
    fresh.mkdir()
    (fresh / "notes.txt").write_text("User notes")
    args[-1] = str(fresh)
    assert main([*args, "--json"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["parameters"]["sasa_points"] == 960
    assert not (fresh / "residues.csv").exists()
    assert not (fresh / "contacts.csv").exists()
    assert (fresh / "notes.txt").read_text() == "User notes"
