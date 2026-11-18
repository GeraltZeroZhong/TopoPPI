import subprocess
import sys
from pathlib import Path

import pytest

from topoppi.examples import copy_example, main


def test_exported_example_matches_synthetic_structure(tmp_path):
    output = tmp_path / "example"
    path = copy_example(output, synthetic=True)
    reference = Path(__file__).parent / "fixtures" / "synthetic_complex.pdb"
    assert path == output / "synthetic_complex.pdb"
    assert path.read_bytes() == reference.read_bytes()
    protein_chains = {line[21] for line in path.read_text().splitlines() if line.startswith("ATOM  ")}
    assert protein_chains == {"A", "B"}


@pytest.mark.parametrize("synthetic", [False, True])
def test_unchanged_example_is_reusable(tmp_path, synthetic):
    first = copy_example(tmp_path, synthetic=synthetic)
    assert first.name == ("synthetic_complex.pdb" if synthetic else "1bvk.pdb")
    assert copy_example(tmp_path, synthetic=synthetic) == first
    assert list(tmp_path.iterdir()) == [first]


def test_edited_example_is_preserved(tmp_path):
    path = copy_example(tmp_path, synthetic=True)
    edited = path.read_bytes() + b"REMARK user annotation\n"
    path.write_bytes(edited)
    with pytest.raises(FileExistsError, match="Choose another output folder"):
        copy_example(tmp_path, synthetic=True)
    assert path.read_bytes() == edited


def test_example_cli_prints_exported_structure(tmp_path, monkeypatch, capsys):
    output = tmp_path / "cli example"
    monkeypatch.setattr("sys.argv", ["topoppi.examples", str(output), "--synthetic"])
    main()
    path = Path(capsys.readouterr().out.strip())
    assert path == output / "synthetic_complex.pdb"
    assert path.is_file()


@pytest.mark.parametrize("destination_is_file", [False, True])
def test_example_cli_preserves_existing_data_with_readable_error(tmp_path, destination_is_file):
    output = tmp_path / "my example"
    if destination_is_file:
        output.write_text("user data\n")
        preserved = output
    else:
        preserved = copy_example(output, synthetic=True)
        preserved.write_bytes(preserved.read_bytes() + b"REMARK user annotation\n")
    original = preserved.read_bytes()

    result = subprocess.run(
        [sys.executable, "-m", "topoppi.examples", str(output), "--synthetic"],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 1
    assert result.stdout == ""
    assert "topoppi.examples:" in result.stderr
    assert "Traceback" not in result.stderr
    assert preserved.read_bytes() == original
    if not destination_is_file:
        assert "Choose another output folder" in result.stderr
