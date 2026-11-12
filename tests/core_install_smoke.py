"""Verify a normal TopoPPI installation, including the default Chemistry plugin.

Run once before and once after installing examples/plugin (with --plugin).
Use --without-chemistry-deps only for a separately provisioned environment with
geometry dependencies and a --no-deps installation. That checks Core isolation;
it is not the dependency set installed by a normal pip install topoppi command.
"""
import argparse
import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import topoppi


def command(*args, success=True):
    result = subprocess.run([sys.executable, "-I", "-m", "topoppi.cli", *map(str, args)],
                            capture_output=True, text=True)
    assert (result.returncode == 0) == success, result.stdout + result.stderr
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plugin", action="store_true")
    parser.add_argument("--without-chemistry-deps", action="store_true")
    parser.add_argument("--allow-scikit-learn", action="store_true",
                        help="Allow scikit-learn supplied by the desktop Conda environment")
    args = parser.parse_args()
    absent = ("sklearn", "torch", "topoppi_affinity", "topoppi_hotspot")
    if args.allow_scikit_learn:
        absent = tuple(name for name in absent if name != "sklearn")
    if args.without_chemistry_deps:
        absent += ("pandas", "prolif", "MDAnalysis", "rdkit")
    for name in absent:
        assert importlib.util.find_spec(name) is None, f"Unexpected optional dependency: {name}"
    startup_imports = {"numpy", "Bio", "trimesh", "matplotlib", "tkinter", "pandas", "prolif", "MDAnalysis", "rdkit"}
    assert not startup_imports.intersection(sys.modules)
    assert topoppi.plugins.available("chemistry")
    assert "chemistry.detect" in topoppi.plugins.info("chemistry").capabilities
    assert topoppi.run("chemistry.backends") == ("prolif",)
    assert not startup_imports.intersection(sys.modules)
    assert not topoppi.plugins.available("affinity")
    assert not topoppi.plugins.available("hotspot")
    from topoppi.errors import PluginUnavailableError

    try:
        topoppi.run("affinity.predict")
    except PluginUnavailableError as exc:
        assert "topoppi-affinity" in str(exc)
    else:
        raise AssertionError("Missing plugin did not raise a TopoPPI-level error")
    assert "topoppi-affinity" in command("affinity", success=False).stderr
    listed = json.loads(command("plugins", "list", "--json").stdout)
    assert {item["name"] for item in listed} == ({"chemistry", "stats"} if args.plugin else {"chemistry"})
    assert command("chemistry", "backends").stdout.strip() == "prolif"

    import numpy as np
    from Bio.PDB import MMCIFIO, PDBParser

    from topoppi.examples import copy_example

    with tempfile.TemporaryDirectory(prefix="topoppi-core-smoke-") as directory:
        root = Path(directory)
        path = copy_example(root, synthetic=True)
        structure = topoppi.load_structure(path)
        assert structure.chain_ids == ("A", "B")
        geometry = topoppi.compute_geometry(topoppi.get_interface(structure, "A", "B"))
        assert all(len(partner.mesh.faces) and partner.patches for partner in geometry.partners)
        assert geometry.contacts()
        native = root / "native.npz"
        topoppi.save_interface(geometry, native)
        restored = topoppi.load_interface(native)
        for first, second in zip(geometry.partners, restored.partners, strict=True):
            np.testing.assert_array_equal(first.mesh.vertices, second.mesh.vertices)
            np.testing.assert_array_equal(first.mesh.faces, second.mesh.faces)
        writer = MMCIFIO()
        writer.set_structure(PDBParser(QUIET=True).get_structure("synthetic", path))
        cif = root / "synthetic.cif"
        writer.save(str(cif))
        cif_geometry = topoppi.compute_geometry(topoppi.get_interface(topoppi.load_structure(cif), "A", "B"))
        assert cif_geometry.contacts() == geometry.contacts()
        command("interface", path, "-A", "A", "-B", "B", "-o", root / "contacts.json")
        command("geometry", cif, "-A", "A", "-B", "B", "-o", root / "cli.native.npz")
        if args.plugin:
            assert topoppi.plugins.info("stats").capabilities == ("stats.describe",)
            summary = topoppi.run("stats.describe", geometry)
            assert summary["atoms"] == [17, 17]
            assert json.loads(command("stats", path, "--chains-a", "A", "--chains-b", "B").stdout) == summary
        path.unlink()
        source_free = topoppi.compute_geometry(topoppi.get_interface(structure, "B", "A"))
        assert source_free.partners[0].chains == ("B",)
        if not args.without_chemistry_deps:
            interactions = topoppi.run("chemistry.detect", geometry)
            assert interactions.engine == "prolif"
            assert (interactions.chain_a, interactions.chain_b) == ("A", "B")
            assert len(interactions.interactions) == 6
            assert all({"res_a_seq", "res_b_seq", "interaction"} <= record.keys()
                       for record in interactions.interactions)
    assert not set(absent).intersection(sys.modules)
    print(json.dumps({"default_chemistry": not args.without_chemistry_deps, "example_plugin": args.plugin,
                      "chemistry_dependency_isolation": args.without_chemistry_deps,
                      "pdb_mmcif_geometry": "passed", "native_roundtrip": "passed",
                      "missing_plugin": "passed", "cli": "passed"}))


if __name__ == "__main__":
    main()
