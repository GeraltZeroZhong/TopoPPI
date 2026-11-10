"""Gradient, thermodynamic symmetry, and actual no-surface input-path tests."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import torch
from torch_geometric.data import Batch, Data

from topoppi_affinity.graph_dataset import GraphDataset
from topoppi_affinity.graph_features import extract_complex_graph, save_graph
from topoppi_affinity.graph_models import AffinityGraphNetwork, GraphModelConfig

OPERATORS = ("gine", "gatv2", "appnp", "mlp")


@pytest.fixture(scope="module", autouse=True)
def bounded_cpu_threads():
    previous = torch.get_num_threads()
    torch.set_num_threads(1)
    yield
    torch.set_num_threads(previous)


@pytest.fixture()
def paired_batch():
    generator = torch.Generator().manual_seed(17)
    graphs = []
    for n_nodes in (6, 5):
        x = torch.rand((n_nodes, 45), generator=generator)
        x[:, :21] = 0.
        x[:, 0] = 1.
        x_mut = x.clone()
        mutation_mask = torch.zeros(n_nodes, dtype=torch.bool)
        mutation_mask[[1, 4]] = True
        x_mut[mutation_mask, 0] = 0.
        x_mut[mutation_mask, 1] = 1.
        x_mut[mutation_mask, 21:26] += .5
        esm = torch.rand((n_nodes, 8), generator=generator)
        esm_mut = esm.clone()
        esm_mut[mutation_mask, :4] += .7
        forward = torch.stack([torch.arange(n_nodes - 1), torch.arange(1, n_nodes)])
        edge_index = torch.cat([forward, forward.flip(0)], dim=1)
        edge_attr = torch.rand((edge_index.shape[1], 8), generator=generator)
        side = torch.zeros(n_nodes, dtype=torch.long)
        side[3:] = 1
        contact_mask = torch.ones(n_nodes, dtype=torch.bool)
        contact_mask[0] = False
        graphs.append(Data(x=x, x_mut=x_mut, esm=esm, esm_mut=esm_mut,
                           mutation_mask=mutation_mask, edge_index=edge_index, edge_attr=edge_attr,
                           side=side, contact_mask=contact_mask, y=torch.tensor([1.3])))
    return Batch.from_data_list(graphs)


def _network(operator, target="ddg", use_plm=True):
    torch.manual_seed(29)
    return AffinityGraphNetwork(GraphModelConfig(architecture=operator, target=target,
                                                hidden=16, layers=2, dropout=.2,
                                                plm_dimension=8, use_plm=use_plm))


@pytest.mark.parametrize("operator", OPERATORS)
@pytest.mark.parametrize("target", ("pkd", "ddg"))
def test_original_operator_paths_backpropagate_finite_gradients(operator, target, paired_batch):
    model = _network(operator, target).train()
    prediction = model(paired_batch)
    assert prediction.shape == (2,)
    loss = (prediction - paired_batch.y).square().mean()
    assert torch.isfinite(loss)
    loss.backward()
    for module in (model.chemistry, model.language, model.layers,
                   model.absolute_head if target == "pkd" else model.mutation_head):
        gradients = [parameter.grad for parameter in module.parameters() if parameter.grad is not None]
        assert gradients
        assert all(torch.isfinite(gradient).all() for gradient in gradients)
        assert any(torch.count_nonzero(gradient) > 0 for gradient in gradients)


@pytest.mark.parametrize("operator", OPERATORS)
def test_absolute_partner_exchange_is_invariant(operator, paired_batch):
    model = _network(operator, "pkd").eval()
    exchanged = paired_batch.clone()
    exchanged.side = 1 - exchanged.side
    with torch.inference_mode():
        torch.testing.assert_close(model(paired_batch), model(exchanged), atol=1e-7, rtol=1e-6)


@pytest.mark.parametrize("operator", OPERATORS)
def test_mutation_zero_reverse_symmetry_and_target_center(operator, paired_batch):
    model = _network(operator).eval()
    noop = paired_batch.clone()
    noop.x_mut = noop.x.clone()
    noop.esm_mut = noop.esm.clone()
    reverse = paired_batch.clone()
    reverse.x, reverse.x_mut = paired_batch.x_mut.clone(), paired_batch.x.clone()
    reverse.esm, reverse.esm_mut = paired_batch.esm_mut.clone(), paired_batch.esm.clone()
    with torch.inference_mode():
        prediction = model(paired_batch)
        assert prediction.abs().max() > 1e-8
        torch.testing.assert_close(model(noop), torch.zeros(2), atol=0., rtol=0.)
        torch.testing.assert_close(model(reverse), -prediction, atol=1e-7, rtol=1e-6)
        # Fitting a mean to ddG labels must never introduce a nonzero identity
        # mutation or destroy antisymmetry when converting back to kcal/mol.
        model.target_center.fill_(7.5)
        torch.testing.assert_close(model(paired_batch), prediction, atol=0., rtol=0.)
        torch.testing.assert_close(model(noop), torch.zeros(2), atol=0., rtol=0.)
        model.target_scale.fill_(2.3)
        torch.testing.assert_close(model(paired_batch), prediction * 2.3)


@pytest.fixture()
def no_surface_batches(tmp_path):
    graph = extract_complex_graph(Path(__file__).parent / "fixtures" / "tiny_complex.pdb",
                                  "A", "B", max_surface_points=16)
    frame = pd.DataFrame([{"source": "fixture", "complex_id": "tiny:A:B", "sample_id": "mutant",
                           "mutations": "GA1W", "ddg": 1.}])
    batches = []
    for index, factor in enumerate((1., 31.)):
        altered = {name: value.copy() for name, value in graph.items()}
        altered["x"][:, graph["surface_feature_mask"]] *= factor
        altered["edge_attr"][:, graph["edge_surface_feature_mask"]] *= factor
        for key in ("surface_normal", "footprint_normal", "footprint_offset", "point_pos", "point_normal",
                    "point_weight", "point_interface_weight"):
            altered[key] *= factor
        path = tmp_path / f"graph_{index}.npz"
        save_graph(altered, path)
        manifest = pd.DataFrame([{"source": "fixture", "complex_id": "tiny:A:B", "graph_path": str(path), "status": "ok"}])
        dataset = GraphDataset(frame, manifest, None, use_surface=False, use_plm=False)
        batches.append(Batch.from_data_list([dataset[0]]))
    return batches


@pytest.mark.parametrize("operator", OPERATORS)
def test_no_surface_dataset_has_no_model_side_channel(operator, no_surface_batches):
    baseline, changed_surface = no_surface_batches
    model = _network(operator, use_plm=False).eval()
    model.config.use_surface = False
    assert baseline.surface_area.count_nonzero() == baseline.node_weight.count_nonzero() == 0
    np.testing.assert_array_equal(baseline.x, changed_surface.x)
    np.testing.assert_array_equal(baseline.edge_attr, changed_surface.edge_attr)
    with torch.inference_mode():
        torch.testing.assert_close(model(baseline), model(changed_surface), atol=0., rtol=0.)


@pytest.mark.parametrize("operator", OPERATORS)
def test_training_pairs_share_dropout_masks(operator, paired_batch):
    model = _network(operator).train()
    noop = paired_batch.clone()
    noop.x_mut, noop.esm_mut = noop.x.clone(), noop.esm.clone()
    torch.testing.assert_close(model(noop), torch.zeros(2), atol=0., rtol=0.)
    reverse = paired_batch.clone()
    reverse.x, reverse.x_mut = paired_batch.x_mut.clone(), paired_batch.x.clone()
    reverse.esm, reverse.esm_mut = paired_batch.esm_mut.clone(), paired_batch.esm.clone()
    torch.manual_seed(51)
    forward = model(paired_batch)
    torch.manual_seed(51)
    backward = model(reverse)
    torch.testing.assert_close(forward, -backward, atol=1e-7, rtol=1e-6)
    # Sharing masks within a pair must still draw a new mask for the next pair.
    assert not torch.equal(backward, model(reverse))
