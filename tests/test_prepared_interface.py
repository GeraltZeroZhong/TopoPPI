"""Shared native geometry must keep chain groups, identity and source domains intact."""
import copy
from pathlib import Path

import numpy as np

from topoppi.interface import prepare_interface
from topoppi.interface_io import load_interface, save_interface
from topoppi.visualization.surface_rendering import surface_geometry

FIXTURE = Path(__file__).parent / "fixtures/tiny_complex.pdb"


def test_group_preparation_roundtrip_preserves_sources_without_uv(tmp_path):
    lines = FIXTURE.read_text().splitlines()
    extra = [line[:21] + "C" + line[22:26] + "A" + line[27:]
             for line in lines if line.startswith("ATOM") and line[21] == "B"]
    path = tmp_path / "group.pdb"
    path.write_text("\n".join([line for line in lines if line not in {"TER", "END"}] + extra + ["END"]))
    prepared = prepare_interface(path, "B,A", "C", source="cleaned", complex_id="group")
    assert prepared.partners[0].chains == ("A", "B")
    assert any(key.endswith("1A") for key in prepared.partners[1].residue_keys)
    original = copy.deepcopy(prepared.partners[0].patches)
    branch = prepared.partners[0].layout_patches()
    branch[0].vertices += 100
    branch[0].metadata["uv"] = np.zeros((len(branch[0].vertices), 2))
    save_interface(prepared, tmp_path / "native.npz")
    restored = load_interface(tmp_path / "native.npz")
    assert restored.metadata() == prepared.metadata()
    assert restored.contacts() == prepared.contacts()
    assert restored.visualizer().interaction_residue_source == "geometric"
    assert restored.visualizer().interaction_partner_map == prepared.contacts()
    assert restored.visualizer()._residue_label("A:GLY:1") == "A:Gly1"
    assert restored.visualizer()._residue_label("B:GLY:1") == "B:Gly1"
    assert restored.visualizer(1)._label_for_residue("chain_b", "C:GLY:1A").startswith(("A:", "B:"))
    for before, after in zip(original, restored.partners[0].patches, strict=True):
        np.testing.assert_array_equal(before.vertices, after.vertices)
        np.testing.assert_array_equal(before.faces, after.faces)
        for name in ("source_face_ids", "source_vertex_ids"):
            np.testing.assert_array_equal(before.metadata[name], after.metadata[name])
        assert "uv" not in after.metadata
    from topoppi.atlas.footprints import mesh_vertex_residue_labels

    patch = restored.partners[0].patches[0]
    geometry = surface_geometry(patch, mesh_vertex_residue_labels(patch, restored.partners[0].residue_keys))
    assert len(geometry["cells"]) == 3 * len(patch.faces)
    assert len(geometry["seam_segments"]) == 0


