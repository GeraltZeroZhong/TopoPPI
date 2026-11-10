"""Grouped partitions, graph normalization and prediction summaries."""
from __future__ import annotations

import pandas as pd
from sklearn.model_selection import GroupKFold

from topoppi_affinity.evaluation import regression_metrics


def development_partition(frame, task, fold, seed=11):
    """Construct three folds grouped by mutation, family or interaction."""
    if 'within' in task:
        units = frame[['complex_id', 'mutations']].drop_duplicates().reset_index(drop=True)
        units['unit'] = units.complex_id + '|' + units.mutations
        groups = units.unit.to_numpy()
        row_units = frame.complex_id + '|' + frame.mutations
    else:
        units = frame[['complex_id', 'family_group', 'initial_interaction_group']].drop_duplicates('complex_id')
        column = 'family_group' if 'family' in task else 'initial_interaction_group'
        groups = units[column].to_numpy()
        row_units = frame.complex_id
        units = units.assign(unit=units.complex_id)
    folds = GroupKFold(n_splits=3, shuffle=True, random_state=seed) if 'within' in task else GroupKFold(n_splits=3)
    train_indices, valid_indices = list(folds.split(units, groups=groups))[fold]
    train = row_units.isin(set(units.iloc[train_indices].unit)).to_numpy()
    valid = row_units.isin(set(units.iloc[valid_indices].unit)).to_numpy()
    return train, valid


def fit_normalization(model, dataset):
    import torch

    from topoppi_affinity.graph_models import physical_transform

    subset = dataset.frame.drop_duplicates(['source', 'complex_id']).index
    total = torch.zeros(model.config.node_dimension, dtype=torch.float64)
    squares = torch.zeros_like(total)
    count = 0
    original_plm = dataset.use_plm
    dataset.use_plm = False
    for i in subset:
        x = physical_transform(dataset[int(i)].x).double()
        total += x.sum(0)
        squares += (x * x).sum(0)
        count += len(x)
    dataset.use_plm = original_plm
    mean = total / count
    std = (squares / count - mean.square()).clamp_min(0.).sqrt()
    std[std < .05] = 1.
    model.feature_center.copy_(mean.float())
    model.feature_scale.copy_(std.float())
    labels = dataset.frame[model.config.target].astype(float).to_numpy()
    model.target_center.fill_(float(labels.mean()) if model.config.target == 'pkd' else 0.)
    model.target_scale.fill_(max(float(labels.std()), 1.))


def predict(model, dataset, device, batch_size=8):
    import torch
    from torch_geometric.loader import DataLoader

    model.eval()
    rows = []
    with torch.inference_mode():
        for data in DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=0):
            data = data.to(device)
            prediction = model(data).cpu().numpy()
            truth = data.y.cpu().numpy()
            rows.extend({'sample_id': sample, 'complex_id': complex_id, 'observed': float(y), 'predicted': float(p)}
                for sample, complex_id, y, p in zip(data.sample_id, data.complex_id, truth, prediction, strict=True))
    frame = pd.DataFrame(rows)
    return frame, regression_metrics(frame.observed, frame.predicted, frame.complex_id)
