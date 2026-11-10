"""Equal-readout feature ablations and label-free footprint permutations."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from topoppi_affinity.benchmark_models import LoadedFeatures
from topoppi_affinity.embedding_baselines import zero_surface_slots


def nested_features(sequence, native, *, view):
    """All views have the same columns, bandwidth rule and learning budget.

    Contact retains the legacy native non-surface readout, including additional
    sequence-local terms. Only surface-versus-contact isolates surface slots.
    Sequence-versus-contact therefore tests the added interface readout as a
    whole, not a unique causal contribution of contact geometry.
    """
    if not np.array_equal(sequence.sample_ids, native.sample_ids):
        raise ValueError('Nested feature blocks must have identical row order')
    if view == 'sequence':
        block = np.zeros_like(native.X)
    elif view == 'contact':
        block = zero_surface_slots(native.X, native.feature_names)
    elif view == 'surface':
        block = native.X
    else:
        raise ValueError(f'Unknown nested view: {view}')
    return LoadedFeatures(sequence.sample_ids, np.concatenate([sequence.X, block], axis=1),
                          np.concatenate([sequence.feature_names, native.feature_names]))


def prepare_nested_features(benchmark_root, sequence_root, output, *, kind, permutation_root=None):
    benchmark_root, sequence_root, output = Path(benchmark_root), Path(sequence_root), Path(output)
    target = 'pkd' if kind == 'absolute' else 'ddg'
    seq_view = 'global' if kind == 'absolute' else 'context'
    native_view = 'plm_surface' if kind == 'absolute' else 'plm_analytic'
    ids = pd.read_csv(benchmark_root/'manifests'/f'{kind}.csv', usecols=['sample_id']).sample_id.to_numpy(str)
    seq = LoadedFeatures.read(sequence_root/f'{target}_sequence_{seq_view}.npz')
    raw = LoadedFeatures.read(benchmark_root/'features'/kind/f'{target}_{native_view}.npz')
    seq = LoadedFeatures(ids, seq.aligned(ids), seq.feature_names)
    native = LoadedFeatures(ids, raw.aligned(ids), raw.feature_names)
    output.mkdir(parents=True, exist_ok=True)
    paths = {}
    for view in ['sequence', 'contact', 'surface']:
        value = nested_features(seq, native, view=view)
        path = output/f'{view}.npz'
        np.savez_compressed(path, X=value.X, sample_ids=ids, feature_names=value.feature_names)
        paths[view] = str(path)
    if permutation_root is not None:
        for seed in [11, 29]:
            perm = LoadedFeatures.read(Path(permutation_root)/f'{target}_permuted_{seed}.npz')
            perm = LoadedFeatures(ids, perm.aligned(ids), perm.feature_names)
            np.testing.assert_allclose(zero_surface_slots(native.X, native.feature_names),
                                       zero_surface_slots(perm.X, perm.feature_names), equal_nan=True, rtol=0, atol=0)
            value = nested_features(seq, perm, view='surface')
            path = output/f'permuted_{seed}.npz'
            np.savez_compressed(path, X=value.X, sample_ids=ids, feature_names=value.feature_names)
            paths[f'permuted_{seed}'] = str(path)
    report = {'kind': kind, 'native_source': str(benchmark_root/'features'/kind/f'{target}_{native_view}.npz'),
              'pure_sequence_view': seq_view, 'n_declared': len(ids), 'dimension': len(value.feature_names),
              'views': paths, 'no_labels_loaded': True,
              'interpretation': 'Surface-contact differs only in surface slots; contact-sequence adds the entire native interface readout'}
    (output/'protocol.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


def prepare_permutations(benchmark_root, output, *, kind, sample_ids=None):
    """Use analytic mutation graphs and original absolute graphs, as prespecified.

    This consumes node pools only. Permuted graphs are never fed to an edge
    network. The source graphs and saved benchmark predictions remain unchanged.
    """
    import torch

    from topoppi_affinity.embedding_baselines import pooled_features
    from topoppi_affinity.footprint_controls import permute_footprints
    from topoppi_affinity.graph_dataset import GraphDataset

    benchmark_root, output = Path(benchmark_root), Path(output)
    target = 'pkd' if kind == 'absolute' else 'ddg'
    frame = pd.read_csv(benchmark_root/'inputs'/f'{kind}.csv', keep_default_na=False)
    if sample_ids is not None:
        frame = frame[frame.sample_id.isin(sample_ids)].copy()
    graph_dir = 'graphs' if kind == 'absolute' else 'graphs_analytic'
    dataset = GraphDataset(frame.assign(**{target: 0.}), benchmark_root/graph_dir/'manifest.csv',
                           benchmark_root/'sequences', max_nodes=None, target=target)
    native_view = 'plm_surface' if kind == 'absolute' else 'plm_analytic'
    native = LoadedFeatures.read(benchmark_root/'features'/kind/f'{target}_{native_view}.npz')
    aligned = native.aligned(frame.sample_id)
    values, ids, coverage = {11: [], 29: []}, [], []
    torch.set_num_threads(1)
    for index, row in enumerate(frame.itertuples()):
        try:
            data = dataset[index]
            original, names = pooled_features(data, target=target, use_surface=True)
            # Detect a curvature/readout mismatch before attributing it to shuffle.
            np.testing.assert_allclose(original, aligned[index], rtol=2e-6, atol=2e-6)
            pending = {}
            for seed in values:
                vector, names = pooled_features(permute_footprints(data, seed), target=target, use_surface=True)
                np.testing.assert_array_equal(zero_surface_slots(vector, names), zero_surface_slots(original, names))
                pending[seed] = vector
            for seed in values:
                values[seed].append(pending[seed])
            ids.append(row.sample_id)
            coverage.append({'sample_id': row.sample_id, 'status': 'ok', 'error': ''})
        except (ValueError, KeyError, OSError) as error:
            coverage.append({'sample_id': row.sample_id, 'status': 'error', 'error': str(error)})
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(coverage).to_csv(output/f'{target}_coverage.csv', index=False)
    if not ids:
        raise ValueError('No successful footprint permutations')
    for seed, matrix in values.items():
        np.savez_compressed(output/f'{target}_permuted_{seed}.npz', X=np.asarray(matrix),
                            sample_ids=np.asarray(ids, dtype=str), feature_names=names)
    report = {'kind': kind, 'declared': len(frame), 'successful': len(ids), 'seeds': list(values),
              'curvature': native_view, 'control': 'Joint surface-column shuffle within partner/contact strata',
              'contact_slots_unchanged': True, 'labels_loaded': False, 'edge_network_evaluated': False}
    (output/f'{target}_protocol.json').write_text(json.dumps(report, indent=2)+'\n')
    return report
