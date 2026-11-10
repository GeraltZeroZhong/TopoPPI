"""Primary-fold isolation, measured WT deduplication, and state algebra."""

import numpy as np
import pandas as pd
import pytest
import torch
from torch_geometric.data import Batch, Data

from topoppi_affinity.graph_training import development_partition
from topoppi_affinity.joint_training import (
    JointAffinityNetwork,
    build_joint_views,
    state_difference_to_ddg,
    wt_anchor_view,
)


def _synthetic_frames():
    mutation = pd.DataFrame({
        'sample_id': [f'm{i}' for i in range(24)], 'complex_id': [f'c{i // 4}' for i in range(24)],
        'family_group': [f'f{i // 4}' for i in range(24)],
        'initial_interaction_group': [f'i{i // 4}' for i in range(24)],
        'mutations': [f'AA{i % 4 + 1}V' for i in range(24)], 'n_mutations': 1,
        'sequence_a': [f'seq{i // 4}' for i in range(24)], 'sequence_b': 'other',
        'reference': 'paper', 'method': 'assay', 'temperature_K': 298., 'kd_wt_M': 1e-8,
        'pkd_wt': 8., 'pkd': 6., 'ddg': 2.,
    })
    absolute = mutation.iloc[::4].copy().reset_index(drop=True)
    absolute['sample_id'] = ['a' + str(i) for i in range(len(absolute))]
    absolute['mutations'] = ''
    return absolute, mutation


def test_wt_anchor_count_is_measurements_not_mutation_rows():
    _, mutation = _synthetic_frames()
    mutation.loc[0, ['temperature_K', 'kd_wt_M', 'pkd_wt']] = [303., 1e-9, 9.]
    anchors = wt_anchor_view(mutation)
    assert len(anchors) == 7
    assert anchors.anchor_repetitions.sum() == len(mutation)
    assert set(anchors.pkd) == {8., 9.}
    assert anchors.mutations.eq('').all()
    assert anchors.ddg.isna().all()
    np.testing.assert_allclose(anchors.pkd, -np.log10(anchors.kd_M))


@pytest.mark.parametrize('task', ['absolute_family', 'mutation_family', 'mutation_within'])
def test_auxiliary_families_and_mt_labels_cannot_cross_primary_validation(task):
    absolute, mutation = _synthetic_frames()
    views = build_joint_views(absolute, mutation, task)
    primary = absolute if task.startswith('absolute') else mutation
    train, validation = development_partition(primary, task, 0)
    assert views.primary_train.sample_id.tolist() == primary.loc[train, 'sample_id'].tolist()
    assert views.primary_validation.sample_id.tolist() == primary.loc[validation, 'sample_id'].tolist()
    held = set(views.primary_validation.family_group)
    for auxiliary in (views.auxiliary_absolute, views.auxiliary_mutation):
        assert held.isdisjoint(auxiliary.family_group)
    assert views.auxiliary_absolute.mutations.eq('').all()
    before = views.auxiliary_absolute[['sample_id', 'pkd']].copy()
    if not task.startswith('absolute'):
        mutation.loc[validation, ['pkd', 'pkd_wt', 'ddg']] = 10000.
    else:
        mutation.loc[mutation.family_group.isin(held), ['pkd', 'pkd_wt', 'ddg']] = 10000.
    after = build_joint_views(absolute, mutation, task).auxiliary_absolute[['sample_id', 'pkd']]
    pd.testing.assert_frame_equal(before, after)


def _pair(a, b):
    n = len(a)
    return Batch.from_data_list([Data(
        x=a.clone(), x_mut=b.clone(), edge_index=torch.tensor([[0, 1, 2, 3], [1, 0, 3, 2]]),
        edge_attr=torch.zeros((4, 8)), side=torch.tensor([0, 0, 1, 1]),
        contact_mask=torch.ones(n, dtype=torch.bool), mutation_mask=torch.ones(n, dtype=torch.bool),
        temperature_K=torch.tensor([298.]), esm=torch.zeros((n, 1152)), esm_mut=torch.zeros((n, 1152)),
    )])


def test_scalar_state_sign_reverse_identity_cycle_and_pair_mask_independence():
    torch.manual_seed(4)
    model = JointAffinityNetwork('state', use_plm=False, hidden=8, layers=1).train()
    a = torch.rand(4, 45)
    b, c = a.clone(), a.clone()
    b[:, :26], c[:, :26] = torch.rand(4, 26), torch.rand(4, 26)
    ab, bc, ac = (_pair(a, b), _pair(b, c), _pair(a, c))
    torch.testing.assert_close(model(ac, 'ddg'), model(ab, 'ddg') + model(bc, 'ddg'), atol=2e-7, rtol=1e-5)
    torch.testing.assert_close(model(ab, 'ddg'), -model(_pair(b, a), 'ddg'), atol=2e-7, rtol=1e-5)
    torch.testing.assert_close(model(_pair(a, a), 'ddg'), torch.zeros(1))
    before = model(ab, 'ddg')
    ab.mutation_mask.zero_()
    torch.testing.assert_close(before, model(ab, 'ddg'))
    assert state_difference_to_ddg(torch.tensor([8.]), torch.tensor([7.]), torch.tensor([298.])).item() > 0
    with pytest.raises(ValueError, match='temperature'):
        state_difference_to_ddg(torch.tensor([8.]), torch.tensor([7.]), torch.tensor([float('nan')]))


