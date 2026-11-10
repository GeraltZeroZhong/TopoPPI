"""Regression scenarios for mistaken destinations in actual CLI workflows."""

from pathlib import Path
from unittest.mock import patch

import pytest

from topoppi import cli


@pytest.mark.parametrize("command", ["interface", "geometry"])
def test_core_export_refuses_to_replace_the_input_structure(command, tmp_path, capsys):
    source = tmp_path / "complex.pdb"
    original = (Path(__file__).parent / "fixtures/tiny_complex.pdb").read_bytes()
    source.write_bytes(original)

    with patch("topoppi.load_structure", side_effect=AssertionError("Started before checking output")):
        code = cli.main([command, str(source), "-A", "A", "-B", "B", "-o", str(source)])

    assert code == 1
    assert "output path separate from the input structure" in capsys.readouterr().err
    assert source.read_bytes() == original


def test_conflicting_render_destinations_preserve_an_existing_figure(tmp_path):
    output = tmp_path / "published.png"
    original = b"previous figure"
    output.write_bytes(original)

    with patch("topoppi.visualization.atlas_io.load_atlas", side_effect=AssertionError("Loaded before validation")):
        code = cli.main([
            "render", str(tmp_path / "source.npz"), "-o", str(output), "--export-atlas", str(output),
        ])

    assert code == 1
    assert output.read_bytes() == original


def test_corrected_contact_destination_succeeds_after_conflict(tmp_path, capsys):
    source = tmp_path / "complex.pdb"
    original = (Path(__file__).parent / "fixtures/tiny_complex.pdb").read_bytes()
    source.write_bytes(original)
    assert cli.main(["interface", str(source), "-o", str(source)]) == 1
    capsys.readouterr()

    output = tmp_path / "results" / "contacts.json"
    assert cli.main(["interface", str(source), "-o", str(output)]) == 0
    assert '"contacts"' in output.read_text()
    assert source.read_bytes() == original


def test_chemistry_export_preserves_structure_and_recovers_with_new_destination(tmp_path, capsys):
    from topoppi_chemistry import cli as chemistry_cli

    source = tmp_path / "complex.pdb"
    original = (Path(__file__).parent / "fixtures/synthetic_complex.pdb").read_bytes()
    source.write_bytes(original)
    assert chemistry_cli.main(["detect", str(source), "--chain-a", "A", "--chain-b", "B", "--output", str(source)]) == 1
    assert "output path separate from the input" in capsys.readouterr().err
    assert source.read_bytes() == original
    output = tmp_path / "contacts.json"
    assert chemistry_cli.main(["detect", str(source), "--chain-a", "A", "--chain-b", "B", "--output", str(output)]) == 0
    assert output.is_file() and source.read_bytes() == original


def test_command_typo_suggests_existing_command(capsys):
    assert cli.main(["geometery", "--help"]) == 1
    message = capsys.readouterr().err
    assert "topoppi geometry" in message
    assert "pip install" not in message


def test_native_geometry_render_explains_atlas_creation(tmp_path, capsys):
    native = tmp_path / "native.npz"
    source = Path(__file__).parent / "fixtures/synthetic_complex.pdb"
    assert cli.main(["geometry", str(source), "-o", str(native)]) == 0
    capsys.readouterr()
    assert cli.main(["render", str(native), "-o", str(tmp_path / "figure.png")]) == 1
    from topoppi.visualization.atlas_io import load_atlas

    with pytest.raises(ValueError, match="--export-atlas"):
        load_atlas(native)
