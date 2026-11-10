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
