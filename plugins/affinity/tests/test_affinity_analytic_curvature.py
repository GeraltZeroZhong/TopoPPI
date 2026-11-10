"""Analytic derivatives, physical scaling, and isolated graph-column changes."""
import json
from pathlib import Path

import numpy as np

from topoppi_affinity.analytic_curvature import (
    GRAPH_PROTOCOL,
    gaussian_derivatives,
    implicit_curvatures,
    replace_graph_curvature,
)
from topoppi_affinity.graph_features import extract_complex_graph


def test_gaussian_sphere_rotation_and_physical_scaling():
    points = np.array([[2., 0., 0.], [0., -2., 0.], [2., 2., 2.] / np.sqrt(3.)])
    _, gradient, hessian = gaussian_derivatives(points, np.zeros((1, 3)))
    gaussian, mean = implicit_curvatures(gradient, hessian)
    np.testing.assert_allclose(gaussian, .25, atol=1e-12)
    np.testing.assert_allclose(mean, .5, atol=1e-12)
    _, gradient, hessian = gaussian_derivatives(3. * points, np.zeros((1, 3)), sigma=3.)
    larger, larger_mean = implicit_curvatures(gradient, hessian)
    np.testing.assert_allclose(larger, gaussian / 9., atol=1e-12)
    np.testing.assert_allclose(larger_mean, mean / 3., atol=1e-12)


def test_sparse_mixture_derivatives_match_finite_differences_and_chunking():
    atoms = np.array([[0., 0., 0.], [1., -.5, .3], [-.5, .4, .7], [20., 0., 0.]])
    point = np.array([[1.2, .7, -.3]])
    density, gradient, hessian = gaussian_derivatives(point, atoms, chunk_size=1)
    step = 1e-5
    for axis in range(3):
        shift = np.eye(3)[axis][None] * step
        plus = gaussian_derivatives(point + shift, atoms)
        minus = gaussian_derivatives(point - shift, atoms)
        np.testing.assert_allclose((plus[0] - minus[0]) / (2. * step), gradient[:, axis], rtol=1e-8)
        np.testing.assert_allclose((plus[1] - minus[1]) / (2. * step), hessian[:, :, axis], rtol=1e-8)
    repeated = gaussian_derivatives(np.repeat(point, 3, axis=0), atoms, chunk_size=2)
    for original, batch in zip((density, gradient, hessian), repeated, strict=True):
        np.testing.assert_allclose(batch, np.repeat(original, 3, axis=0))


def test_graph_refinement_changes_only_two_curvature_columns():
    structure = Path(__file__).parent / "fixtures/tiny_complex.pdb"
    graph = extract_complex_graph(structure, "A", "B", max_surface_points=64)
    refined = replace_graph_curvature(graph, structure)
    direct = extract_complex_graph(structure, "A", "B", max_surface_points=64,
                                   curvature_method="analytic_gaussian")
    columns = np.ones(45, dtype=bool)
    columns[[34, 35]] = False
    np.testing.assert_array_equal(graph["x"][:, columns], refined["x"][:, columns])
    np.testing.assert_allclose(refined["x"], direct["x"], rtol=1e-6)
    for name in graph:
        if name not in {"x", "metadata_json"}:
            np.testing.assert_array_equal(graph[name], refined[name])
    assert json.loads(str(refined["metadata_json"]))["protocol"] == GRAPH_PROTOCOL
    assert json.loads(str(graph["metadata_json"]))["protocol"] != GRAPH_PROTOCOL


def test_pooled_patch_matches_full_recomputation_for_both_targets():
    import torch
    from torch_geometric.data import Data

    from topoppi_affinity.analytic_curvature import update_pooled_curvature
    from topoppi_affinity.embedding_baselines import pooled_features
    from topoppi_affinity.graph_features import mutant_node_features

    structure = Path(__file__).parent / "fixtures/tiny_complex.pdb"
    old = extract_complex_graph(structure, "A", "B", max_surface_points=64)
    new = replace_graph_curvature(old, structure)
    rng = np.random.default_rng(8)
    esm = rng.normal(size=(len(old["x"]), 4)).astype(np.float32)
    mutant_esm = esm + .1
    for target in ("pkd", "ddg"):
        vectors, names = [], None
        for graph in (old, new):
            mutant, _, mask = mutant_node_features(graph, "GA1W")
            data = Data(x=torch.from_numpy(graph["x"]), x_mut=torch.from_numpy(mutant),
                        esm=torch.from_numpy(esm), esm_mut=torch.from_numpy(mutant_esm),
                        side=torch.from_numpy(graph["side"]), contact_mask=torch.from_numpy(graph["x"][:, 26] <= 10.),
                        mutation_mask=torch.from_numpy(mask))
            vector, names = pooled_features(data, target=target, use_surface=True)
            vectors.append(vector)
        patched = update_pooled_curvature(vectors[0], names, new, target=target, mutations="GA1W")
        np.testing.assert_array_equal(patched, vectors[1])
