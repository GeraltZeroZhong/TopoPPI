"""Analytical, independently constructed representation and descriptor inputs."""

from __future__ import annotations

import copy
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from topoppi_representation.schema import ComplexSource, Patch
from topoppi_representation.source import (
    canonical_source_bytes,
    load_source,
    native_surface,
    reference_quadrature,
    save_source,
)


def source():
    xy = np.array([[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]])
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    patches = []
    names = [("ASP", "VAL", "GLY", "TYR"), ("LYS", "PHE", "GLN", "UNK")]
    keys = [f"side{s}:A:{aa}:{i + 1}" for s, group in enumerate(names) for i, aa in enumerate(group)]
    for side in (0, 1):
        uv = xy[faces].copy()
        uv[1] += [5.0, 0.0]
        patches.append(
            Patch(
                side,
                0,
                np.column_stack([xy, np.full(4, 2.0 * side)]),
                faces.copy(),
                np.arange(4) + 4 * side,
                {"ordinary": uv, "aware": xy[faces]},
            )
        )
    return ComplexSource(sample_id="toy", residue_keys=keys, amino_acids=list("DVGYKFQX"), patches=patches)


def test_import_keeps_analysis_and_codecs_lazy():
    code = "import sys; import topoppi_representation; import topoppi_representation.api; assert not {'torch','sklearn','pandas','DracoPy','zstandard','pywt'} & sys.modules.keys()"
    subprocess.run([sys.executable, "-c", code], check=True, env=dict(os.environ))


def test_source_roundtrip_and_quadrature_mass(tmp_path):
    value = source()
    path = save_source(value, tmp_path / "source.npz")
    restored = load_source(path)
    assert canonical_source_bytes(restored) == canonical_source_bytes(value)
    a, b = native_surface(value), reference_quadrature(value, subdivision=3)
    np.testing.assert_allclose(a.mass.sum(axis=0), b.mass.sum(axis=0), atol=1e-12)
    assert a.weights.sum() == pytest.approx(8.0)


def test_physical_queries_partner_swap_and_translation():
    from topoppi_representation.descriptors.surface import physical_targets

    first = native_surface(source())
    other = copy.deepcopy(first)
    other.side = 1 - other.side
    other.xyz += [100.0, -30.0, 8.0]
    a, b = physical_targets(first), physical_targets(other)
    np.testing.assert_allclose(a["X"], b["X"], atol=1e-14)
    assert a["O_by_side"]["0"]["L_H"] == b["O_by_side"]["1"]["L_H"]
    assert sum(a["X"]) == pytest.approx(1.0)


def test_direct_preserves_missing_and_only_stored_queries():
    pytest.importorskip("zstandard")
    from topoppi_representation import decode_direct, encode_direct
    from topoppi_representation.codecs.direct import query_direct

    code = encode_direct({"a": 2.5, "missing": None, "nan": float("nan")})
    assert code.nbytes == len(code.blob)
    values = decode_direct(code.blob)
    assert values == {"a": 2.5, "missing": None, "nan": None}
    assert query_direct(values, "unrequested") == {"value": None, "available": False}


@pytest.mark.parametrize("layout", ["ordinary", "aware"])
def test_uv_carrier_area_and_standalone_decode(layout):
    pytest.importorskip("zstandard")
    pytest.importorskip("pywt")
    from topoppi_representation import decode_uv, encode_uv

    original = source()
    code = encode_uv(original, layout=layout, target_cells=16, quant_step=0)
    expected = native_surface(original).mass.sum(axis=0)
    original.patches.clear()
    decoded = decode_uv(code.blob)
    np.testing.assert_allclose(decoded.mass.sum(axis=0), expected, atol=1e-6)
    assert decoded.residue_keys == original.residue_keys
    assert np.isfinite(decoded.xyz).all()
    assert np.all(decoded.side[decoded.edges[:, 0]] == decoded.side[decoded.edges[:, 1]])


def test_native_codec_retains_mass_row_identity():
    pytest.importorskip("zstandard")
    pytest.importorskip("DracoPy")
    from topoppi_representation import decode_native, encode_native

    original = source()
    decoded = decode_native(encode_native(original).blob)
    np.testing.assert_allclose(decoded.mass.sum(axis=0), native_surface(original).mass.sum(axis=0), atol=1e-6)
    assert decoded.faces.shape == (4, 3)


def test_physical_dictionary_and_sasa_recovery_preserve_schema(monkeypatch):
    pytest.importorskip("pandas")
    import topoppi
    from topoppi_representation.descriptors import sasa
    from topoppi_representation.descriptors.physical import physical_descriptors, recover_from_sasa, residue_physics

    structure = topoppi.load_structure(Path(__file__).parent / "fixtures/tiny_complex.pdb")
    interface = topoppi.get_interface(structure, "A", "B")
    monkeypatch.setattr(sasa, "PDBLoader", lambda *args: pytest.fail("Prepared atoms must be reused"))
    values = physical_descriptors(interface, n_points=100)
    residues, pairs = residue_physics(interface, n_points=100)
    recovered = recover_from_sasa(residues)
    assert len(values) == 477
    assert len(recovered) == 430
    for name, value in recovered.items():
        assert value == pytest.approx(values[name])
    swapped = topoppi.get_interface(structure, "B", "A")
    assert values == physical_descriptors(swapped, n_points=100)
    assert len(pairs) == values["physical.contact_pairs.total"]


def test_adapter_reuses_core_identity_and_explicit_layouts():
    from types import SimpleNamespace

    import trimesh
    from Bio.PDB import Atom, Chain, Residue

    from topoppi_representation import from_geometry

    chain = Chain.Chain("A")
    residue = Residue.Residue((" ", 7, "B"), "VAL", " ")
    chain.add(residue)
    atom = Atom.Atom("CA", np.array([0.0, 0.0, 0.0]), 0.0, 1.0, " ", " CA ", 1, element="C")
    residue.add(atom)
    mesh = trimesh.Trimesh(vertices=[[0, 0, 0], [1, 0, 0], [0, 1, 0]], faces=[[0, 1, 2]], process=False)
    mesh.metadata["source_atom_indices"] = np.zeros(3, dtype=int)
    partner = SimpleNamespace(atoms=[atom], patches=[mesh], chains=("A",))
    geometry = SimpleNamespace(partners=[partner], analysis_id="toy")
    layout = np.array([[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]])
    adapted = from_geometry(geometry, layouts={(0, 0): {"ordinary": layout}})
    assert adapted.residue_keys == ["side0:A:VAL:7B"]
    np.testing.assert_array_equal(adapted.patches[0].uv["ordinary"], layout)
    assert adapted.family_id == ""
    assert mesh.metadata["source_atom_indices"].tolist() == [0, 0, 0]


def test_recovery_import_has_no_torch():
    pytest.importorskip("sklearn")
    pytest.importorskip("pandas")
    code = "import sys; from topoppi_representation.analysis.physical_recovery import PhysicalTargetReader; assert 'torch' not in sys.modules"
    subprocess.run([sys.executable, "-c", code], check=True, env=dict(os.environ))


@pytest.mark.parametrize("codec", ["native", "uv", "direct"])
def test_representation_comparison_script(tmp_path, codec):
    pytest.importorskip("zstandard")
    pytest.importorskip("DracoPy")
    pytest.importorskip("pywt")
    import json

    path = save_source(source(), tmp_path / "input.npz")
    script = Path(__file__).parents[1] / "scripts/analysis/compare_representations.py"
    output = tmp_path / codec
    subprocess.run(
        [sys.executable, str(script), str(path), "--codec", codec, "--cells", "16", "--output", str(output)], check=True
    )
    report = json.loads((output / "comparison.json").read_text())
    assert report["encoded_bytes"] == (output / "encoded.bin").stat().st_size
    assert report["canonical_source_bytes"] == len(canonical_source_bytes(source()))
    assert report["encode_cpu_seconds"] >= 0


def test_analysis_script_paths_from_explicit_synthetic_partitions(tmp_path):
    pytest.importorskip("torch")
    import pandas as pd

    from topoppi_representation.analysis.graphs import (
        TARGET_NAMES,
        make_probe_sample,
        save_probe_sample,
        targets_from_physics,
    )
    from topoppi_representation.descriptors.surface import physical_targets

    paths, targets, features, ids = [], [], [], []
    for i in range(9):
        item = source()
        item.sample_id = f"toy{i}"
        item.family_id = f"group{i}"
        item.patches[1].xyz[:, 2] += i * 0.1
        y, mask = targets_from_physics(physical_targets(native_surface(item)))
        sample = make_probe_sample(item, y, mask, k=2)
        path = tmp_path / f"sample{i}.npz"
        save_probe_sample(sample, path)
        split = ["train", "validation", "test"][i // 3]
        paths.append({"path": path.name, "split": split})
        targets.append(
            {"sample_id": item.sample_id, "group": item.family_id, "split": split, **dict(zip(TARGET_NAMES, y, strict=False))}
        )
        features.append([float(i), float(i * i)])
        ids.append(item.sample_id)
    pd.DataFrame(paths).to_csv(tmp_path / "graphs.csv", index=False)
    pd.DataFrame(targets).to_csv(tmp_path / "targets.csv", index=False)
    np.savez(tmp_path / "features.npz", X=features, sample_ids=ids, feature_names=["position", "position_squared"])
    scripts = Path(__file__).parents[1] / "scripts/analysis"
    subprocess.run(
        [
            sys.executable,
            str(scripts / "physical_recovery.py"),
            str(tmp_path / "features.npz"),
            str(tmp_path / "targets.csv"),
            "--output",
            str(tmp_path / "linear"),
        ],
        check=True,
    )
    subprocess.run(
        [
            sys.executable,
            str(scripts / "finite_readout.py"),
            str(tmp_path / "graphs.csv"),
            "--epochs",
            "2",
            "--width",
            "4",
            "--output",
            str(tmp_path / "graph"),
        ],
        check=True,
    )
    for folder in ["linear", "graph"]:
        result = pd.read_csv(tmp_path / folder / "predictions.csv")
        assert result.sample_id.tolist() == ["toy6", "toy7", "toy8"]
        assert np.isfinite(result[list(TARGET_NAMES)].to_numpy()).all()
