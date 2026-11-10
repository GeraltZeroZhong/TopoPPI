"""Preserve source data when figure and preflight destinations conflict."""

from pathlib import Path
from unittest.mock import patch

import pytest

from topoppi import cli
from topoppi.benchmark_cli import main as benchmark_main
from topoppi.config import TopoPPIRunConfig, VisualizationConfig
from topoppi.errors import ConfigurationError


@pytest.mark.parametrize("destination", ["report", "annotations", "interactions"])
def test_mapping_atlas_preserves_other_workflow_files(tmp_path, destination):
    structure = Path(__file__).parent / "fixtures/synthetic_complex.pdb"
    figure = tmp_path / "figure.png"
    paths = {"report": figure.with_suffix(".topoppi.json"),
             "annotations": tmp_path / "values.csv", "interactions": tmp_path / "interactions.json"}
    paths["annotations"].write_text("residue,value\nA:2,1.0\n")
    paths["interactions"].write_text('{}')
    paths["report"].write_text('Previous run report')
    original = {key: path.read_bytes() for key, path in paths.items()}
    config = TopoPPIRunConfig(pdb_file=str(structure), output_file=str(figure),
                             atlas_output=str(paths[destination]), prolif_file=str(paths["interactions"]),
                             visualization=VisualizationConfig(annotation_file=str(paths["annotations"])))
    with pytest.raises(ConfigurationError, match="atlas output path separate"):
        config.validate()
    assert {key: path.read_bytes() for key, path in paths.items()} == original
    assert not figure.exists()


def test_render_keeps_annotation_input_when_selected_as_atlas_destination(tmp_path, caplog):
    values = tmp_path / "values.csv"
    original = "residue,value\nA:2,1.0\n"
    values.write_text(original)
    assert cli.main(["render", str(tmp_path / "saved.npz"), "--annotation-file", str(values),
                     "--export-atlas", str(values), "-o", str(tmp_path / "image.png")]) == 1
    assert "annotation input and atlas output" in caplog.text
    assert values.read_text() == original
    assert not (tmp_path / "image.png").exists()


@pytest.mark.parametrize("command", ["preflight", "preflight-sensitivity"])
def test_preflight_report_keeps_its_input_file(command, tmp_path, capsys):
    source = tmp_path / "input.json"
    original = '{"experiment": "user settings"}'
    source.write_text(original)
    with patch("topoppi.benchmark_cli.load_benchmark_config", side_effect=AssertionError("Loaded before path check")), \
         patch("topoppi.benchmark_cli.SensitivityBenchmarkRunner", side_effect=AssertionError("Loaded before path check")):
        assert benchmark_main([command, str(source), "--output-json", str(source)]) == 2
    assert "report path separate from the input" in capsys.readouterr().err
    assert source.read_text() == original
