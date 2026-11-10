"""Joint absolute-affinity and mutation model components."""
from __future__ import annotations

import json
from dataclasses import dataclass

import numpy as np
import pandas as pd
import torch
from torch_geometric.loader import DataLoader

from topoppi_affinity.evaluation import regression_metrics
from topoppi_affinity.graph_dataset import GraphDataset
from topoppi_affinity.graph_models import AffinityGraphNetwork, GraphModelConfig, physical_transform
from topoppi_affinity.graph_training import development_partition
from topoppi_affinity.models import complex_weights
from topoppi_affinity.thermodynamics import R_KCAL

WT_MEASUREMENT_KEYS = ('sequence_a', 'sequence_b', 'reference', 'method', 'temperature_K', 'kd_wt_M')


def wt_anchor_view(mutations):
    """One WT measurement per sequence/experiment/T/Kd key, never MT labels."""
    records = []
    for number, (_, group) in enumerate(mutations.sort_values('sample_id').groupby(
            list(WT_MEASUREMENT_KEYS), sort=True, dropna=False)):
        record = group.iloc[0].to_dict()
        record.update(sample_id=f'wt_anchor:{number:04d}', pkd=float(record['pkd_wt']),
                      mutations='', mutations_pdb='', n_mutations=0, kd_M=float(record['kd_wt_M']), ddg=np.nan,
                      dg_kcal_mol=-R_KCAL * float(record['temperature_K']) * np.log(10.) * float(record['pkd_wt']),
                      label_origin='deduplicated_experimental_WT_Kd', anchor_members=json.dumps(group.sample_id.tolist()),
                      anchor_repetitions=len(group), anchor_kind='SKEMPI_WT')
        records.append(record)
    if not records:
        return mutations.iloc[:0].assign(anchor_members='', anchor_repetitions=0, anchor_kind='SKEMPI_WT')
    return pd.DataFrame(records)


@dataclass
class JointViews:
    primary_train: pd.DataFrame
    primary_validation: pd.DataFrame
    auxiliary_absolute: pd.DataFrame
    auxiliary_mutation: pd.DataFrame
    excluded_auxiliary: pd.DataFrame
    metadata: dict


def build_joint_views(absolute, mutation, task, fold=0, seed=11):
    """Preserve primary folds; remove every held-out family from auxiliary data."""
    primary = absolute if task.startswith('absolute') else mutation
    train, validation = development_partition(primary, task, fold, seed)
    families = set(primary.loc[validation, 'family_group'])
    allowed_absolute = absolute.loc[~absolute.family_group.isin(families)].copy()
    allowed_mutation = mutation.loc[~mutation.family_group.isin(families)].copy()
    anchors = wt_anchor_view(allowed_mutation)
    if task.startswith('absolute'):
        auxiliary_absolute, auxiliary_mutation = anchors, allowed_mutation
        excluded = mutation.loc[mutation.family_group.isin(families), ['sample_id', 'family_group']].copy()
    else:
        # A within-family WT anchor is also excluded, even when its ddg row is
        # in primary training. This deliberately applies the strict aux rule.
        auxiliary_absolute = pd.concat([allowed_absolute, anchors], ignore_index=True)
        auxiliary_mutation = mutation.iloc[:0].copy()
        excluded = pd.concat([
            absolute.loc[absolute.family_group.isin(families), ['sample_id', 'family_group']],
            mutation.loc[mutation.family_group.isin(families), ['sample_id', 'family_group']],
        ], ignore_index=True)
    excluded['exclusion_reason'] = 'primary_validation_family'
    metadata = {
        'task': task, 'fold': fold, 'seed': seed, 'validation_families': sorted(families),
        'wt_measurement_keys': list(WT_MEASUREMENT_KEYS), 'n_wt_anchors': len(anchors),
        'n_wt_source_rows': len(allowed_mutation), 'mt_absolute_supervision': False,
        'auxiliary_exclusion': 'All primary validation families excluded, including within-task WT auxiliary anchors',
        'primary_train_sample_ids': primary.loc[train, 'sample_id'].tolist(),
        'primary_validation_sample_ids': primary.loc[validation, 'sample_id'].tolist(),
        'auxiliary_absolute_sample_ids': auxiliary_absolute.sample_id.tolist(),
        'auxiliary_mutation_sample_ids': auxiliary_mutation.sample_id.tolist(),
        'auxiliary_wt_source_sample_ids': allowed_mutation.sample_id.tolist(),
        'confirmation_evaluated': False,
    }
    return JointViews(primary.loc[train].copy().reset_index(drop=True),
                      primary.loc[validation].copy().reset_index(drop=True),
                      auxiliary_absolute.reset_index(drop=True), auxiliary_mutation.reset_index(drop=True),
                      excluded, metadata)


class JointGraphDataset(GraphDataset):
    def __getitem__(self, index):
        data = super().__getitem__(index)
        value = self.frame.iloc[index].get('temperature_K', '')
        data.temperature_K = torch.tensor([float(value) if value != '' and pd.notna(value) else float('nan')])
        return data


def state_difference_to_ddg(wild_type_pkd, mutant_pkd, temperature_K):
    if not torch.isfinite(temperature_K).all() or not (temperature_K > 0).all():
        raise ValueError('State differences require each observation\'s positive reported temperature')
    return -R_KCAL * np.log(10.) * temperature_K * (mutant_pkd - wild_type_pkd)


class JointAffinityNetwork(AffinityGraphNetwork):
    def __init__(self, mode='independent', *, use_surface=True, use_plm=True, hidden=64, layers=2):
        if mode not in {'independent', 'state'}:
            raise ValueError('Joint mode must be independent or state')
        super().__init__(GraphModelConfig(architecture='gine', hidden=hidden, layers=layers,
            dropout=0., use_surface=use_surface, use_plm=use_plm, target='ddg', max_nodes=None))
        self.mode = mode
        self.register_buffer('pkd_center', torch.tensor(0.))
        self.register_buffer('pkd_scale', torch.tensor(1.))

    def state_score(self, data, *, mutant=False):
        return self.absolute_readout(self.encode(data, mutant=mutant), data) * self.pkd_scale + self.pkd_center

    def forward(self, data, target='pkd'):
        if target == 'pkd':
            return self.state_score(data)
        if target != 'ddg':
            raise ValueError('Target must be pkd or ddg')
        if self.mode == 'independent':
            return super().forward(data)
        return state_difference_to_ddg(self.state_score(data), self.state_score(data, mutant=True), data.temperature_K)


def _fit_training_normalization(model, datasets):
    total = torch.zeros(model.config.node_dimension, dtype=torch.float64)
    squares, count, seen = torch.zeros_like(total), 0, set()
    by_target = {'pkd': [], 'ddg': []}
    for dataset in datasets:
        by_target[dataset.target].append(dataset.frame)
        use_plm = dataset.use_plm
        dataset.use_plm = False
        for index, row in dataset.frame.iterrows():
            key = (row.source, row.complex_id)
            if key in seen:
                continue
            seen.add(key)
            x = physical_transform(dataset[index].x).double()
            total += x.sum(0)
            squares += x.square().sum(0)
            count += len(x)
        dataset.use_plm = use_plm
    center = total / count
    scale = (squares / count - center.square()).clamp_min(0.).sqrt()
    scale[scale < .05] = 1.
    model.feature_center.copy_(center.float())
    model.feature_scale.copy_(scale.float())
    scales = {}
    for target, frames in by_target.items():
        if not frames:
            scales[target] = 1.
            continue
        frame = pd.concat(frames, ignore_index=True)
        weights, labels = complex_weights(frame.complex_id), frame[target].to_numpy(dtype=float)
        mean = float(np.average(labels, weights=weights))
        standard_deviation = max(float(np.sqrt(np.average((labels - mean)**2, weights=weights))), 1.)
        scales[target] = standard_deviation
        if target == 'pkd':
            model.pkd_center.fill_(mean)
            model.pkd_scale.fill_(standard_deviation)
        else:
            model.target_scale.fill_(standard_deviation)
    return scales


def predict_joint(model, dataset, device):
    model.eval()
    rows = []
    with torch.inference_mode():
        for data in DataLoader(dataset, batch_size=1, shuffle=False, num_workers=0):
            data = data.to(device)
            prediction = model(data, dataset.target).cpu().numpy()
            rows.extend({'sample_id': sample, 'complex_id': group, 'observed': float(observed), 'predicted': float(predicted)}
                        for sample, group, observed, predicted in zip(data.sample_id, data.complex_id,
                            data.y.cpu().numpy(), prediction, strict=True))
    frame = pd.DataFrame(rows)
    return frame, regression_metrics(frame.observed, frame.predicted, frame.complex_id)
