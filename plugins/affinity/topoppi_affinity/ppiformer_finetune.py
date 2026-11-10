"""Same-split fine-tuning of the released structural masked-modeling weights.

The initialization was selected using SKEMPI validation ddG in the original
release. Local train/test isolation does not erase that selection exposure.
This is a fixed-budget adaptation, not the author's three-fold tuning protocol.
"""
from __future__ import annotations

import json
import random
import time
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd

METHOD = 'PPIformer-local-finetune-selection-exposed'
EXPOSURE = 'released_SSL_checkpoint_selected_on_SKEMPI_validation; local_training_targets_only'


def enable_activation_checkpointing(model):
    """Recompute encoder activations to fit the observed 2120-node interfaces."""
    from torch.utils.checkpoint import checkpoint

    for block in model.encoder.layers.blocks:
        for layer in block:
            layer.forward = partial(checkpoint, layer.forward, use_reentrant=False)


def _author_types():
    import ppiformer.definitions as definitions

    definitions.SE3TRANSFORMER_REQUIRED = False
    from ppiformer.data.transforms import DDGLabelPretransform, PDBToPyGPretransform
    from ppiformer.tasks.node import DDGPPIformer
    return DDGPPIformer, DDGLabelPretransform, PDBToPyGPretransform


def load_initialization(repository, device):
    import torch

    cls, _, _ = _author_types()
    checkpoint = torch.load(Path(repository)/'weights/masked_modeling.ckpt', map_location='cpu', weights_only=False)
    parameters = dict(checkpoint['hyper_parameters'])
    parameters.update(kind='masked_marginals', correction=False, test_csv_path=None)
    model = cls(**parameters)
    model.load_state_dict(checkpoint['state_dict'], strict=True)
    return model.to(device), checkpoint['hyper_parameters']['run_cfg']


def prepare_graphs(frame, output):
    """Cache label-free original author interface graphs, with explicit failures."""
    import torch

    from topoppi_affinity.expanded_references import prepare_ppiformer_interface

    _, _, transform_type = _author_types()
    transform = transform_type(k=None)
    output = Path(output)
    paths, failures = {}, {}
    for complex_id, subset in frame.groupby('complex_id', sort=False):
        work = output/complex_id.replace(':', '_')
        path = work/'graph.pt'
        try:
            if not path.exists():
                work.mkdir(parents=True, exist_ok=True)
                interface = prepare_ppiformer_interface(subset.iloc[0].to_dict(), work)
                graph = transform(interface)
                torch.save(graph, path)
            paths[complex_id] = path
        except Exception as error:
            failures[complex_id] = f'{type(error).__name__}: {error}'
    return paths, failures


def labelled_graph(graph, mutations, ddg=np.nan):
    _, label_type, _ = _author_types()
    frame = pd.DataFrame([{'#Pdb': graph.path.stem, 'Mutation(s)': mutations, 'ddG': ddg}])
    value = label_type(df=frame, mut_col='Mutation(s)', strict_wt=True)(graph.clone())
    if value.n_muts != 1:
        raise ValueError('Mutation absent or WT mismatched in original 10A interface')
    return value


def fit_ppiformer(train, heldout, output, repository, graph_root, *, device='cuda', epochs=40, seed=11):
    """Fixed epochs, one mutation per complex per epoch, Adam and accumulation32.

    No early stopping, outer validation or test-label access is performed.
    The 40-epoch budget is fixed before outcomes and is not a best-model search.
    """
    import torch
    from torch_geometric.data import Batch

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    heldout = heldout.drop(columns=['ddg', 'pkd', 'dg_kcal_mol'], errors='ignore').copy()
    if set(train.sample_id) & set(heldout.sample_id):
        raise ValueError('Local PPIformer train/test IDs overlap')
    torch.set_num_threads(2)
    torch.manual_seed(seed)
    rng = random.Random(seed)
    started = time.monotonic()
    model, run_cfg = load_initialization(repository, device)
    enable_activation_checkpointing(model)
    identity_columns = ['sample_id', 'complex_id', 'pdb_id', 'structure_path', 'chains_a', 'chains_b', 'mutations']
    paths, failures = prepare_graphs(pd.concat([train[identity_columns], heldout[identity_columns]]), graph_root)
    graphs = {key: torch.load(path, map_location='cpu', weights_only=False) for key, path in paths.items()}
    training_groups, coverage = {}, []
    for row in train.to_dict('records'):
        state = dict(sample_id=row['sample_id'], status='ok', error='')
        try:
            if row['complex_id'] in failures:
                raise ValueError(failures[row['complex_id']])
            labelled_graph(graphs[row['complex_id']], row['mutations'], row['ddg'])
            training_groups.setdefault(row['complex_id'], []).append(row)
        except Exception as error:
            state.update(status='error', error=f'{type(error).__name__}: {error}')
        coverage.append(state)
    pd.DataFrame(coverage).to_csv(output/'training_coverage.csv', index=False)
    if not training_groups:
        raise ValueError('No usable local PPIformer training complexes')
    optimizer = torch.optim.Adam(model.parameters(), lr=.0003, weight_decay=0.)
    used_ids, losses = set(), []
    model.train()
    for epoch in range(epochs):
        keys = list(training_groups)
        rng.shuffle(keys)
        optimizer.zero_grad()
        epoch_losses = []
        for i, key in enumerate(keys):
            row = rng.choice(training_groups[key])
            graph = labelled_graph(graphs[key], row['mutations'], row['ddg'])
            batch = Batch.from_data_list([graph], follow_batch=['mut_ddg']).to(device)
            pred = model(batch)
            loss = torch.nn.functional.mse_loss(pred, batch.mut_ddg)
            if not torch.isfinite(loss):
                raise ValueError('Nonfinite training loss; no result checkpoint written')
            group_size = min(32, len(keys) - (i//32)*32)
            (loss/group_size).backward()
            if (i+1) % 32 == 0 or i+1 == len(keys):
                optimizer.step()
                optimizer.zero_grad()
            epoch_losses.append(float(loss.detach().cpu()))
            used_ids.add(row['sample_id'])
        losses.append({'epoch': epoch+1, 'train_mse': float(np.mean(epoch_losses))})
        print(f'PPIformer epoch {epoch+1}/{epochs}: {len(keys)} complexes', flush=True)
    model.eval()
    rows = []
    for row in heldout.to_dict('records'):
        result = dict(sample_id=row['sample_id'], complex_id=row['complex_id'], predicted=np.nan,
                      method=METHOD, target='ddg', status='error', error='', training_overlap=EXPOSURE,
                      comparison_scope='local_retraining_with_author_selection_exposure')
        try:
            if row['complex_id'] in failures:
                raise ValueError(failures[row['complex_id']])
            graph = labelled_graph(graphs[row['complex_id']], row['mutations'])
            batch = Batch.from_data_list([graph], follow_batch=['mut_ddg']).to(device)
            with torch.inference_mode():
                value = float(model(batch)[0].cpu())
            if not np.isfinite(value):
                raise ValueError('Nonfinite fine-tuned prediction')
            result.update(predicted=value, status='ok')
        except Exception as error:
            result['error'] = f'{type(error).__name__}: {error}'
        rows.append(result)
    result = pd.DataFrame(rows)
    result.to_csv(output/'predictions.csv', index=False)
    train[['sample_id']].to_csv(output/'training_ids.csv', index=False)
    pd.DataFrame({'sample_id': sorted(used_ids)}).to_csv(output/'sampled_training_ids.csv', index=False)
    heldout[['sample_id']].to_csv(output/'test_ids.csv', index=False)
    pd.DataFrame(losses).to_csv(output/'training_loss.csv', index=False)
    torch.save({'state_dict': model.state_dict(), 'initialization': str(Path(repository)/'weights/masked_modeling.ckpt'),
                'epochs': epochs, 'seed': seed}, output/'model.pt')
    (output/'initialization_run_cfg.yaml').write_text(run_cfg)
    metadata = dict(method=METHOD, epochs=epochs, optimizer='Adam', learning_rate=.0003,
                    mutation_samples_per_complex_per_epoch=1, gradient_accumulation=32,
                    usable_training_complexes=len(training_groups), sampled_training_rows=len(used_ids),
                    activation_checkpointing=True, training_size_exclusions=0,
                    declared_test=len(heldout), successful_test=int(result.status.eq('ok').sum()),
                    heldout_targets_read=False, training_exposure=EXPOSURE,
                    selection='Fixed final epoch; no validation or test outcome selection',
                    seconds=time.monotonic()-started)
    (output/'fit.json').write_text(json.dumps(metadata, indent=2)+'\n')
    return result
