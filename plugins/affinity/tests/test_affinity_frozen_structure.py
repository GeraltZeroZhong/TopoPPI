"""Actual graph alignment and complete core coverage for frozen RDE transfer."""

import numpy as np

from topoppi_affinity.frozen_structure import align_author_residues, covering_patches


def test_covering_patches_are_local_complete_deterministic_and_bounded():
    positions = np.random.default_rng(31).normal(size=(300, 3))
    patches = covering_patches(positions)
    covered = np.zeros(len(positions), int)
    for core, patch in patches:
        assert len(patch) == 128
        assert len(core) == 32
        assert set(core).issubset(patch)
        assert np.all(patch[:-1] < patch[1:])
        covered[core] += 1
        distances = np.linalg.norm(positions[:, None] - positions[core][None], axis=-1).min(axis=1)
        outside = np.setdiff1d(np.arange(len(positions)), patch)
        assert distances[patch].max() <= distances[outside].min() + 1e-12
    assert covered.min() >= 1
    for (core, patch), (again_core, again_patch) in zip(patches, covering_patches(positions), strict=True):
        np.testing.assert_array_equal(core, again_core)
        np.testing.assert_array_equal(patch, again_patch)


def test_residue_mapping_distinguishes_chain_insertion_and_identity():
    graph = {
        'chain': np.asarray(['B', 'A', 'A', 'C']), 'residue_number': np.asarray([10, 10, 10, 12]),
        'insertion': np.asarray(['', 'A', '', '']), 'aa_index': np.asarray([1, 0, 2, 3]),
        'residue_key': np.asarray(['B:10:CYS', 'A:10A:ALA', 'A:10:ASP', 'C:12:GLU']),
    }
    mapping = {('A', 10, ' '): 0, ('A', 10, 'A'): 1, ('B', 10, ' '): 2}
    indices, errors = align_author_residues(graph, mapping, np.asarray([4, 0, 1]))
    np.testing.assert_array_equal(indices, [2, 1, -1, -1])
    assert errors[0]['reason'] == 'identity_mismatch:D:F'
    assert errors[1]['reason'] == 'omitted_by_author_parser'


