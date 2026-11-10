"""RDE omission handling, exact mapping, and information-matched controls."""
import numpy as np
import pytest

from topoppi_affinity.features import AA, residue_properties
from topoppi_affinity.structure_features import (
    align_structure_embedding,
    matched_structure_features,
    structure_pool,
)


@pytest.fixture()
def graph():
    x = np.zeros((6, 45), dtype=np.float32)
    for index, aa in enumerate("ACDDEF"):
        x[index, AA.index(aa)] = 1.
        x[index, 21:26] = list(residue_properties(aa).values())
    x[:, 26] = [2., 5., 20., 3., 7., 20.]
    return {"x": x, "side": np.array([0, 0, 0, 1, 1, 1]),
            "chain": np.array(["A", "A", "A", "B", "B", "B"]),
            "residue_number": np.array([1, 2, 3, 1, 2, 3]), "insertion": np.array([""] * 6),
            "residue_key": np.array(["A:1", "A:2", "A:3", "B:1", "B:2", "B:3"]),
            "aa_index": np.array([AA.index(aa) for aa in "ACDDEF"])}


def test_exact_key_alignment_with_missing_values(graph):
    original = np.arange(18, dtype=np.float32).reshape(6, 3)
    original[2] = np.nan
    order = np.array([4, 2, 3, 0, 1, 5])
    bundle = {"embedding": original[order], "residue_key": graph["residue_key"][order],
              "encoded_mask": np.array([True, True, False, True, True, True])[order]}
    aligned, mask = align_structure_embedding(graph, bundle)
    np.testing.assert_array_equal(aligned, original)
    assert not mask[2]


def test_absolute_valid_only_and_partner_invariance(graph):
    embedding = np.arange(18, dtype=np.float32).reshape(6, 3)
    embedding[1] = np.nan
    valid = np.array([True, False, True, True, True, True])
    vector, names, coverage = structure_pool(graph, embedding, valid, target="pkd")
    assert np.isfinite(vector).all()
    assert len(vector) == len(names) == 12
    assert coverage["side0_contact_encoded"] == 1
    np.testing.assert_allclose(vector[:3], embedding[[0, 2]].mean(axis=0) + embedding[3:].mean(axis=0))
    graph["side"] = 1 - graph["side"]
    np.testing.assert_array_equal(vector, structure_pool(graph, embedding, valid, target="pkd")[0])
    valid[:3] = False
    with pytest.raises(ValueError, match="no encoded side"):
        structure_pool(graph, embedding, valid, target="pkd")


def test_mutation_requires_encoded_sites_and_is_conditionally_odd(graph):
    from topoppi_affinity.graph_features import mutant_node_features

    embedding = np.arange(18, dtype=np.float32).reshape(6, 3)
    valid = np.array([True, True, False, True, True, True])
    embedding[2] = np.nan
    vector, names, _ = structure_pool(graph, embedding, valid, target="ddg", mutations="CA2W,EB2A")
    assert np.isfinite(vector).all() and len(vector) == len(names) == 26 * 3
    reverse = graph.copy()
    reverse["x"], reverse["aa_index"], _ = mutant_node_features(graph, "CA2W,EB2A")
    np.testing.assert_allclose(vector, -structure_pool(reverse, embedding, valid, target="ddg",
                                                       mutations="WA2C,AB2E")[0])
    np.testing.assert_array_equal(structure_pool(graph, embedding, valid, target="ddg", mutations="")[0], 0.)
    with pytest.raises(ValueError, match="RDE missing mutation residues: A:3"):
        structure_pool(graph, embedding, valid, target="ddg", mutations="DA3A")


def test_four_controls_share_dimensions_and_nonablated_values():
    base = np.array([[1., 2., 3.], [4., 5., 6.]])
    rde = np.array([[7., 8.], [9., 10.]])
    values, names = matched_structure_features(base, np.array(["plm.a", "chemical.b", "surface.c"]),
                                               rde, np.array(["rde.a", "rde.b"]))
    assert len(names) == 5
    for matrix in values.values():
        assert matrix.shape == (2, 5)
        np.testing.assert_array_equal(matrix[:, :2], base[:, :2])
    np.testing.assert_array_equal(values["rde_control_contact"][:, 2:], 0.)
    np.testing.assert_array_equal(values["rde_context"][:, 2], 0.)
    np.testing.assert_array_equal(values["rde_control_surface"][:, :3], base)
    np.testing.assert_array_equal(values["rde_control_surface"][:, 3:], 0.)
    np.testing.assert_array_equal(values["rde_context_surface"][:, 3:], rde)
