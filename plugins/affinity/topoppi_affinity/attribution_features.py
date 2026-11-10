"""Explicit absolute-affinity blocks and residue-weight correspondence controls."""
from __future__ import annotations

import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from topoppi_affinity.benchmark_models import LoadedFeatures

SEEDS = (11, 29, 47, 71, 101)
VIEWS = {
    'S': ('S',), 'B': ('S', 'I', 'C'), 'BA': ('S', 'I', 'C', 'A'),
    'BAK': ('S', 'I', 'C', 'A', 'K'), 'BAW': ('S', 'I', 'C', 'A', 'W'),
    'BAKW': ('S', 'I', 'C', 'A', 'K', 'W'),
}
WEIGHT_CONTROLS = ('uniform', 'sasa', *(f'perm_{s}' for s in SEEDS), *(f'aa_perm_{s}' for s in SEEDS))


def save_features(path, ids, matrix, names):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, sample_ids=np.asarray(ids, dtype=str), X=np.asarray(matrix, dtype=np.float32),
                        feature_names=np.asarray(names, dtype=str))


def split_blocks(sequence, native, legacy):
    """Discard duplicated native global PLM; select 3 area versus 12 shape columns."""
    ids = sequence.sample_ids
    native_x, legacy_x = native.aligned(ids), legacy.aligned(ids)
    names = native.feature_names.astype(str)
    stats = np.array([name.startswith(('surface.mean.', 'surface.sum.', 'surface.std.')) for name in names])
    area = np.array([int(name.rsplit('.', 1)[1]) % 15 < 3 if selected else False
                     for name, selected in zip(names, stats, strict=True)])
    masks = {'I': np.char.startswith(names, 'plm.contact.'), 'A': area, 'K': stats & ~area,
             'W': np.char.startswith(names, 'surface.area_weighted_plm.')}
    blocks = {'S': sequence}
    for key, mask in masks.items():
        blocks[key] = LoadedFeatures(ids, native_x[:, mask], names[mask])
    mask = np.char.startswith(legacy.feature_names.astype(str), 'contact.')
    blocks['C'] = LoadedFeatures(ids, legacy_x[:, mask], legacy.feature_names[mask])
    if any(not block.X.shape[1] for block in blocks.values()):
        raise ValueError('Required attribution block has no columns')
    return blocks


def prepare_blocks(root, benchmark, prep):
    root, benchmark, prep = Path(root), Path(benchmark), Path(prep)
    inputs = pd.read_csv(benchmark/'inputs/absolute.csv', keep_default_na=False)
    sequence = LoadedFeatures.read(prep/'sequence/absolute/pkd_sequence_global.npz')
    native = LoadedFeatures.read(benchmark/'features/absolute/pkd_plm_surface.npz')
    legacy = LoadedFeatures.read(benchmark/'features/absolute/pkd_legacy_v1.npz')
    blocks = split_blocks(sequence, native, legacy)
    for name, block in blocks.items():
        save_features(root/'features/absolute'/f'{name}.npz', block.sample_ids, block.X, block.feature_names)
    # PPB's actual released sequence state remains its sequence input.
    ppb = pd.read_csv(prep/'ppb/structure_manifest.csv', keep_default_na=False)
    actual = LoadedFeatures.read(prep/'ppb/features/pkd_sequence_global.npz')
    for name, block in blocks.items():
        matrix = actual.aligned(ppb.sample_id) if name == 'S' else block.aligned(ppb.cached_sample_id)
        names = actual.feature_names if name == 'S' else block.feature_names
        save_features(root/'features/ppb'/f'{name}.npz', ppb.sample_id, matrix, names)
    dimensions = {key: block.X.shape[1] for key, block in blocks.items()}
    report = {'absolute_rows': len(inputs), 'ppb_rows': len(ppb), 'block_dimensions': dimensions,
              'duplicate_native_global_plm_included': False, 'labels_used': False,
              'A': 'Per-partner mean/sum/std of surface area, footprint area, footprint fraction',
              'K': 'Remaining 12 surface descriptors: patch, curvature, normal dispersion and extent',
              'C': 'Every contact.* column in the frozen legacy matrix; no legacy surface columns',
              'W': 'Native residue footprint-area-weighted ESMC, normalized by max(area sum, 1)'}
    (root/'features/protocol.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


def weighted_pool(embedding, side, weights):
    embedding, weights = np.asarray(embedding, dtype=np.float64), np.asarray(weights, dtype=np.float64)
    vectors = [(embedding[side == partner] * weights[side == partner, None]).sum(axis=0)
               / max(float(weights[side == partner].sum()), 1.) for partner in (0, 1)]
    return np.concatenate((vectors[0]+vectors[1], np.abs(vectors[0]-vectors[1]))).astype(np.float32)


def permuted_weights(weights, side, contact, aa, seed, *, condition_aa):
    rng = np.random.default_rng(seed)
    result = weights.copy()
    for partner in (0, 1):
        for interface in (False, True):
            mask = (side == partner) & (contact == interface)
            for identity in np.unique(aa[mask]) if condition_aa else [None]:
                indices = np.flatnonzero(mask if identity is None else mask & (aa == identity))
                result[indices] = weights[rng.permutation(indices)]
    return result


def buried_sasa_weights(row, graph):
    """Shrake-Rupley isolated-partner minus complex SASA on the same selected atoms."""
    from Bio.PDB import Chain, Model, Residue
    from Bio.PDB.SASA import ShrakeRupley
    from topoppi.atlas.footprints import source_atom_residue_labels
    from topoppi.io.io_loader import PDBLoader

    from topoppi_affinity.features import _chains, _group_data

    loader = PDBLoader(row['structure_path'])
    groups = sorted((_chains(row['chains_a']), _chains(row['chains_b'])))
    atoms = [_group_data(loader, group)[1] for group in groups]
    model = Model.Model(0)
    node_lookup = {key: i for i, key in enumerate(graph['residue_key'].astype(str))}
    for partner, source in enumerate(atoms):
        chain = Chain.Chain(str(partner))
        model.add(chain)
        residues = {}
        labels = source_atom_residue_labels(source)
        for atom, label in zip(source, labels, strict=True):
            node = node_lookup[str(label)]
            if node not in residues:
                residue = Residue.Residue((' ', node+1, ' '), atom.get_parent().resname, ' ')
                chain.add(residue)
                residues[node] = residue
            residues[node].add(atom.copy())
    sr = ShrakeRupley(probe_radius=1.4, n_points=100)
    sr.compute(model, level='R')
    bound = {residue.id[1]-1: residue.sasa for residue in model.get_residues()}
    weights = np.zeros(len(graph['x']))
    for chain in model:
        sr.compute(chain, level='R')
        for residue in chain:
            node = residue.id[1]-1
            weights[node] = max(0., residue.sasa-bound[node])
    return weights


def _weight_job(row, graph_manifest, sequences, folder):
    import torch

    from topoppi_affinity.graph_dataset import GraphDataset
    from topoppi_affinity.graph_features import load_graph

    torch.set_num_threads(1)
    sample = row['sample_id']
    path = Path(folder)/(sample.replace(':', '_')+'.npz')
    if path.exists():
        with np.load(path) as old:
            return json.loads(str(old['report']))
    started = time.monotonic()
    data = GraphDataset(pd.DataFrame([row]).assign(pkd=0.), graph_manifest, sequences,
                        max_nodes=None, target='pkd')[0]
    embedding, side = data.esm.numpy(), data.side.numpy()
    weights, contact = data.x[:, 31].numpy(), data.contact_mask.numpy()
    aa = data.x[:, :21].numpy().argmax(axis=1)
    values, coverage = {}, []
    for name in WEIGHT_CONTROLS:
        error = ''
        try:
            if name == 'uniform':
                changed = contact.astype(float)
            elif name == 'sasa':
                graphs = pd.read_csv(graph_manifest, keep_default_na=False)
                graph_path = graphs.loc[graphs.source.eq(row['source']) & graphs.complex_id.eq(row['complex_id']),
                                        'graph_path'].iloc[0]
                changed = buried_sasa_weights(row, load_graph(graph_path))
            else:
                changed = permuted_weights(weights, side, contact, aa, int(name.rsplit('_', 1)[1]),
                                           condition_aa=name.startswith('aa_'))
            values[name] = weighted_pool(embedding, side, changed)
            altered = ~np.isclose(changed, weights, rtol=0, atol=1e-8)
            coverage.append(dict(sample_id=sample, control=name, status='ok', error='',
                changed_residue_fraction=float(altered.mean()),
                original_area_on_changed_fraction=float(weights[altered].sum()/max(weights.sum(), 1.)),
                total_weight=float(changed.sum())))
        except (ValueError, KeyError, OSError, IndexError) as exc:
            values[name] = np.full(2*embedding.shape[1], np.nan, dtype=np.float32)
            error = f'{type(exc).__name__}: {exc}'
            coverage.append(dict(sample_id=sample, control=name, status='error', error=error))
    report = dict(sample_id=sample, coverage=coverage, seconds=time.monotonic()-started)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **values, report=np.asarray(json.dumps(report)))
    return report


def prepare_weights(root, benchmark, *, workers=4, sample_ids=None):
    root, benchmark = Path(root), Path(benchmark)
    frame = pd.read_csv(benchmark/'inputs/absolute.csv', keep_default_na=False)
    if sample_ids is not None:
        frame = frame[frame.sample_id.isin(sample_ids)]
    folder = root/'weight_cache'
    folder.mkdir(parents=True, exist_ok=True)
    coverage = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_weight_job, row, str(benchmark/'graphs/manifest.csv'), str(benchmark/'sequences'), folder)
                   for row in frame.to_dict('records')]
        for i, future in enumerate(as_completed(futures), 1):
            result = future.result()
            coverage.extend(result['coverage'])
            if i % 100 == 0 or i == len(futures):
                print(f'Weight inputs {i}/{len(futures)}', flush=True)
    table = pd.DataFrame(coverage)
    table.to_csv(root/'weight_coverage.csv', index=False)
    for control in WEIGHT_CONTROLS:
        matrix = []
        for sample in frame.sample_id:
            with np.load(folder/(sample.replace(':', '_')+'.npz')) as source:
                matrix.append(source[control])
        names = [f'weight.{control}.{i}' for i in range(len(matrix[0]))]
        save_features(root/'features/absolute'/f'W_{control}.npz', frame.sample_id, matrix, names)
    report = {'rows': len(frame), 'controls': list(WEIGHT_CONTROLS),
              'counts': table.groupby(['control', 'status']).size().to_dict(),
              'labels_used': False, 'native_summaries_A_K_unchanged': True}
    report['counts'] = {f'{a}:{b}': int(n) for (a, b), n in report['counts'].items()}
    (root/'weight_protocol.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


def load_view(root, cohort, view):
    blocks = VIEWS[view] if view in VIEWS else (*VIEWS['BAK'], 'W_'+view.removeprefix('control_'))
    values = [LoadedFeatures.read(Path(root)/'features'/cohort/f'{block}.npz') for block in blocks]
    ids = values[0].sample_ids
    return LoadedFeatures(ids, np.concatenate([block.aligned(ids) for block in values], axis=1),
                          np.concatenate([block.feature_names for block in values]))
