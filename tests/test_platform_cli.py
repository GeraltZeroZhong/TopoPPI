"""Platform CLI checks use synthetic inputs and isolated provider registrations."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from topoppi import cli, plugins
from topoppi.plugins import PluginMetadata, PluginRegistry


@pytest.fixture(autouse=True)
def isolated_registry(monkeypatch):
    registry = PluginRegistry([])
    monkeypatch.setattr(plugins, "_default_registry", registry)
    return registry


def test_plugin_list_and_unavailable_error(capsys):
    assert cli.main(["plugins", "list", "--json"]) == 0
    assert json.loads(capsys.readouterr().out) == []
    assert cli.main(["affinity"]) == 1
    assert "topoppi-affinity" in capsys.readouterr().err
    assert cli.main(["plugins", "info", "affinity"]) == 1
    assert "topoppi-affinity" in capsys.readouterr().err


def test_plugin_command_and_help_are_dynamic(isolated_registry, capsys):
    calls = []

    class Toy:
        metadata = PluginMetadata("toy", "1.0", "Toy CLI")

        def register(self, registry):
            registry.register_command("toy", lambda argv: calls.append(argv) or 0)

    assert "Installed plugin commands:" not in cli.build_parser().format_help()
    isolated_registry.register_plugin(Toy())
    assert "Installed plugin commands: toy" in " ".join(cli.build_parser().format_help().split())
    assert cli.main(["toy", "--synthetic"]) == 0
    assert calls == [["--synthetic"]]
    assert cli.main(["plugins", "info", "toy", "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["description"] == "Toy CLI"


@pytest.mark.parametrize("arguments,exit_code,stream", [(["--help"], 0, "out"), (["--unknown"], 2, "err")])
def test_plugin_usage_includes_namespace_and_restores_argv(isolated_registry, monkeypatch, capsys,
                                                          arguments, exit_code, stream):
    class Toy:
        metadata = PluginMetadata("toy", "1.0", "Toy CLI")

        def register(self, registry):
            registry.register_command("toy", lambda argv: argparse.ArgumentParser().parse_args(argv))

    original = ["/env/bin/topoppi", "toy", *arguments]
    monkeypatch.setattr(sys, "argv", original)
    isolated_registry.register_plugin(Toy())
    assert cli.main(["toy", *arguments]) == exit_code
    assert "usage: topoppi toy" in getattr(capsys.readouterr(), stream)
    assert sys.argv is original


def test_existing_extensionless_structure_path_still_uses_mapping(tmp_path):
    path = tmp_path / "structure"
    path.write_text("synthetic placeholder")
    with patch("topoppi.cli.run_interface_mapping") as run:
        assert cli.main([str(path)]) == 0
    assert run.call_args.args[0].pdb_file == str(path)


def test_core_interface_command_writes_contacts_from_synthetic_structure(tmp_path):
    fixture = Path(__file__).parent / "fixtures/tiny_complex.pdb"
    output = tmp_path / "contacts.json"
    assert cli.main(["interface", str(fixture), "-A", "A", "-B", "B", "-o", str(output)]) == 0
    value = json.loads(output.read_text())
    assert value["contacts"]
    assert all(key.startswith("A:") for key in value["contacts"])


def test_geometry_dispatch_forwards_core_settings_without_renderer(monkeypatch, tmp_path, capsys):
    import topoppi

    structure, interface = object(), object()
    geometry = SimpleNamespace(metadata=lambda: {"kind": "toy_geometry"})
    calls = []
    monkeypatch.setattr(topoppi, "load_structure", lambda path: structure)
    monkeypatch.setattr(topoppi, "get_interface", lambda value, a, b, **kwargs: interface)

    def compute(value, **settings):
        assert value is interface
        calls.append(settings)
        return geometry

    monkeypatch.setattr(topoppi, "compute_geometry", compute)
    monkeypatch.setattr(topoppi, "save_interface", lambda value, path: Path(path))
    output = tmp_path / "geometry.npz"
    with patch("topoppi.cli.run_interface_mapping", side_effect=AssertionError("mapping invoked")):
        assert cli.main(["geometry", "toy.pdb", "--res", "1.5", "--cutoff", "5", "-o", str(output)]) == 0
    assert calls[0]["surface"].grid_resolution == 1.5
    assert calls[0]["topology"].distance_cutoff == 5.
    assert json.loads(capsys.readouterr().out)["output"] == str(output)
