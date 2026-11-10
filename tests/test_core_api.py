"""Public structure, contact and geometry contracts on independent toy inputs."""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import topoppi
from topoppi import (
    GeometryRepresentation,
    PreparedInterface,
    compute_geometry,
    get_interface,
    load_interface,
    load_structure,
    residue_dual_boundaries,
    save_interface,
    vertex_area_masses,
)

FIXTURE = Path(__file__).parent / "fixtures" / "synthetic_complex.pdb"


def test_structure_and_contacts_load_without_surface_or_task_modules():
    script = """
import sys
import topoppi
assert not any(n in sys.modules for n in ('numpy', 'Bio', 'trimesh', 'torch', 'matplotlib.pyplot'))
structure = topoppi.load_structure(sys.argv[1])
interface = topoppi.get_interface(structure, 'A', 'B')
assert interface.contacts()
assert not any(n in sys.modules for n in (
    'trimesh', 'topoppi.mesh.surface', 'topoppi.pipeline', 'topoppi.optimization.optcuts',
    'matplotlib.pyplot', 'tkinter', 'torch', 'topoppi_affinity', 'topoppi_hotspot'))
"""
    environment = dict(os.environ, PYTHONPATH=str(Path(topoppi.__file__).resolve().parents[1]),
                       PYTHONDONTWRITEBYTECODE="1")
    subprocess.run([sys.executable, "-c", script, str(FIXTURE.resolve())],
                   env=environment, check=True, capture_output=True, text=True)


def _brute_force_contacts(interface, cutoff):
    left, right = interface.partners
    distances = np.linalg.norm(left.coords[:, None, :] - right.coords[None, :, :], axis=2)
    result = {}
    for i, j in np.argwhere(distances <= cutoff):
        values = result.setdefault(str(left.residue_keys[i]), {})
        label = str(right.residue_keys[j])
        values[label] = values.get(label, 0) + 1
    return result


def test_loaded_structure_survives_source_removal_and_preserves_contact_rule(tmp_path):
    source = tmp_path / "synthetic.pdb"
    source.write_bytes(FIXTURE.read_bytes())
    structure = load_structure(source, source="synthetic", complex_id="toy")
    source.unlink()

    interface = get_interface(structure, "B", "A", contact_distance=4.5, analysis_id="toy-analysis")
    assert [p.chains for p in interface.partners] == [("B",), ("A",)]
    assert interface.source["path"] == str(source)
    assert interface.contacts() == _brute_force_contacts(interface, 4.5)
    reverse = {right: {left: count for left, values in interface.contacts().items()
                       for partner, count in values.items() if partner == right}
               for right in {right for values in interface.contacts().values() for right in values}}
    assert interface.contacts(side=1) == reverse

    geometry = compute_geometry(interface)
    assert isinstance(geometry, GeometryRepresentation)
    assert GeometryRepresentation is PreparedInterface
    assert geometry.interface is interface
    assert geometry.analysis_id == "toy-analysis"
    assert all(partner.patches for partner in geometry.partners)
    assert all(not hasattr(partner, "mesh") for partner in interface.partners)

    restored = load_interface(save_interface(geometry, tmp_path / "native.npz"))
    assert restored.metadata() == geometry.metadata()
    assert restored.interface.contact_distance_angstrom == 4.5
    assert restored.contacts() == geometry.contacts()
    assert restored.contacts(cutoff=6.0) == _brute_force_contacts(interface, 6.0)
    for before, after in zip(geometry.partners, restored.partners, strict=True):
        np.testing.assert_array_equal(before.coords, after.coords)
        np.testing.assert_array_equal(before.mesh.vertices, after.mesh.vertices)
        np.testing.assert_array_equal(before.mesh.faces, after.mesh.faces)


def test_contact_distance_and_surface_face_cutoff_are_independent(tmp_path):
    interface = get_interface(load_structure(FIXTURE), "A", "B", contact_distance=0.01)
    assert interface.contacts() == {}
    geometry = compute_geometry(interface, sides=(0,))
    assert geometry.partners[0].patches
    assert geometry.partners[1].mesh is None
    assert geometry.topology_config.distance_cutoff == 4.0
    assert geometry.contacts() == {}


def test_default_native_record_keeps_existing_version_and_metadata(tmp_path):
    geometry = compute_geometry(get_interface(load_structure(FIXTURE), "A", "B"), sides=(0,))
    destination = save_interface(geometry, tmp_path / "native.npz")
    with np.load(destination, allow_pickle=False) as arrays:
        record = json.loads(str(arrays["document_json"]))
    assert (record["format"], record["version"]) == ("topoppi-native-interface", 1)
    assert "contact_distance_angstrom" not in record["interface"]
    restored = load_interface(destination)
    assert restored.metadata() == geometry.metadata()
    assert restored.contacts() == geometry.contacts()


def test_partner_selection_retains_complete_groups_and_rejects_overlap(tmp_path):
    lines = FIXTURE.read_text().splitlines()
    extra = [line[:21] + "C" + line[22:26] + "A" + line[27:]
             for line in lines if line.startswith("ATOM") and line[21] == "B"]
    source = tmp_path / "group.pdb"
    source.write_text("\n".join([line for line in lines if line not in {"TER", "END"}] + extra + ["END"]))
    structure = load_structure(source)
    interface = get_interface(structure, "C", "B,A")
    assert [partner.chains for partner in interface.partners] == [("C",), ("A", "B")]
    assert [len(partner.atoms) for partner in interface.partners] == [17, 34]
    assert any(key.endswith("1A") for key in interface.partners[0].residue_keys)
    with pytest.raises(ValueError, match="overlap"):
        get_interface(structure, "A,B", "B,C")


def test_shared_native_area_and_boundary_integrals_have_analytic_values():
    import trimesh

    mesh = trimesh.Trimesh(vertices=[[0., 0., 0.], [2., 0., 0.], [0., 2., 0.]],
                           faces=[[0, 1, 2]], process=False)
    np.testing.assert_allclose(vertex_area_masses(mesh), [2 / 3, 2 / 3, 2 / 3])
    boundary = residue_dual_boundaries(mesh, np.asarray([0, 1, 0]))
    expected = np.sqrt(5) / 3 + np.sqrt(2) / 3
    assert boundary == pytest.approx({(0, 1): expected})
