"""Label-free, resumable native inputs for the fixed affinity benchmark recipes."""
from __future__ import annotations

import json
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from topoppi_affinity.features import parse_mutations
from topoppi_affinity.graph_features import GRAPH_PROTOCOL, extract_complex_graph, load_graph, save_graph
from topoppi_affinity.sequence_embeddings import chain_records, extract_embeddings


def read_inputs(root):
    return pd.concat([pd.read_csv(Path(root) / 'inputs' / f'{task}.csv', keep_default_na=False)
                      for task in ('absolute', 'mutation')], ignore_index=True)


def _graph_task(task):
    row, path, cache, previous = task
    started = time.monotonic()
    record = {key: row[key] for key in ('source', 'complex_id', 'structure_path', 'chains_a', 'chains_b')}
    record['graph_path'] = str(path)
    try:
        cached = False
        candidate = path if path.exists() else previous
        if candidate is not None and candidate.exists():
            graph = load_graph(candidate)
            metadata = json.loads(str(graph['metadata_json']))
            cached = metadata['protocol'] == GRAPH_PROTOCOL and metadata['max_surface_points_per_side'] == 2048
            if cached:
                record['graph_path'] = str(candidate)
        if not cached:
            graph = extract_complex_graph(row['structure_path'], row['chains_a'], row['chains_b'], cache, 2048)
            save_graph(graph, path)
        record.update(status='ok', error='', cached=cached, n_nodes=len(graph['x']))
    except (ValueError, KeyError, OSError, IndexError) as error:
        record.update(status='error', error=f'{type(error).__name__}: {error}', cached=False)
    record['seconds'] = time.monotonic() - started
    return record


def prepare_graphs(root, reuse_root, *, workers=2):
    root, reuse_root = Path(root), Path(reuse_root)
    folder = root / 'graphs'
    folder.mkdir(parents=True, exist_ok=True)
    paths = [folder / 'manifest.csv', *sorted((reuse_root / 'graphs').glob('*manifest.csv'))]
    sources = [pd.read_csv(path, keep_default_na=False) for path in paths if path.exists()]
    cached = pd.concat(sources, ignore_index=True) if sources else pd.DataFrame()
    lookup = {}
    if len(cached):
        for row in cached[cached.status.eq('ok')].to_dict('records'):
            key = tuple(str(row.get(name, '')) for name in
                        ('source', 'complex_id', 'structure_path', 'chains_a', 'chains_b'))
            lookup.setdefault(key, Path(row['graph_path']))
    frame = read_inputs(root).drop_duplicates(['source', 'complex_id'])
    tasks = []
    for row in frame.to_dict('records'):
        name = row['complex_id'].replace(':', '_').replace(',', '+') + '.npz'
        key = tuple(str(row[name]) for name in ('source', 'complex_id', 'structure_path', 'chains_a', 'chains_b'))
        path = folder / row['source'] / name
        tasks.append((row, path, root / 'native_cache', lookup.get(key)))
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        for future in as_completed([executor.submit(_graph_task, task) for task in tasks]):
            row = future.result()
            rows.append(row)
            pd.DataFrame(rows).to_csv(folder / 'manifest.csv', index=False)
            print(f'Graph {len(rows)}/{len(tasks)}: {row["complex_id"]} {row["status"]}', flush=True)
    return pd.DataFrame(rows)


def prepare_sequences(root, reuse_root):
    """Keep per-row mapping failures instead of replacing a missing mutant by WT."""
    root, reuse_root = Path(root), Path(reuse_root)
    folder = root / 'sequences'
    folder.mkdir(parents=True, exist_ok=True)
    for name in ('sequences.csv', 'wild_type_chains.csv', 'mutant_chains.csv'):
        if not (folder / name).exists() and (reuse_root / 'sequences' / name).exists():
            shutil.copyfile(reuse_root / 'sequences' / name, folder / name)
    (folder / 'embeddings').mkdir(exist_ok=True)
    existing = pd.read_csv(folder / 'sequences.csv', keep_default_na=False) if (folder / 'sequences.csv').exists() else pd.DataFrame()
    registry = dict(zip(existing.sequence, existing.sequence_id, strict=True)) if len(existing) else {}
    wt_path, mt_path = folder / 'wild_type_chains.csv', folder / 'mutant_chains.csv'
    wild = pd.read_csv(wt_path, keep_default_na=False).to_dict('records') if wt_path.exists() else []
    mutants = pd.read_csv(mt_path, keep_default_na=False).to_dict('records') if mt_path.exists() else []
    coverage, requested = [], set()

    def register(sequence):
        if sequence not in registry:
            registry[sequence] = f'sequence_{len(registry):06d}'
        requested.add(registry[sequence])
        return registry[sequence]

    for (source, complex_id), group in read_inputs(root).groupby(['source', 'complex_id'], sort=True):
        row = group.iloc[0]
        try:
            chains = chain_records(row.structure_path, row.chains_a, row.chains_b)
        except (ValueError, KeyError, OSError, IndexError) as error:
            coverage.extend({'sample_id': sample, 'status': 'error', 'error': str(error)} for sample in group.sample_id)
            continue
        for chain, records in chains.items():
            sequence = ''.join(value[0] for value in records.values())
            wild.append({'source': source, 'complex_id': complex_id, 'chain': chain,
                         'sequence_id': register(sequence), 'sequence': sequence, 'residue_keys': json.dumps(list(records))})
        for row in group.itertuples():
            try:
                mutation_list = parse_mutations(row.mutations) if row.mutations else []
                pending = []
                for chain in sorted({mutation['chain'] for mutation in mutation_list}):
                    records = list(chains[chain].values())
                    positions = {(value[1], value[2]): index for index, value in enumerate(records)}
                    sequence = [value[0] for value in records]
                    for mutation in mutation_list:
                        if mutation['chain'] != chain:
                            continue
                        index = positions[(int(mutation['position']), mutation['insertion'])]
                        if sequence[index] != mutation['wt']:
                            raise ValueError(f'{row.sample_id}: structural sequence WT mismatch')
                        sequence[index] = mutation['mt']
                    pending.append({'sample_id': row.sample_id, 'complex_id': complex_id,
                                    'chain': chain, 'sequence_id': register(''.join(sequence))})
                mutants.extend(pending)
                coverage.append({'sample_id': row.sample_id, 'status': 'ok', 'error': ''})
            except (ValueError, KeyError, IndexError) as error:
                coverage.append({'sample_id': row.sample_id, 'status': 'error', 'error': str(error)})
    pd.DataFrame(wild).drop_duplicates(['source', 'complex_id', 'chain'], keep='last').to_csv(wt_path, index=False)
    pd.DataFrame(mutants, columns=['sample_id', 'complex_id', 'chain', 'sequence_id']).drop_duplicates(
        ['sample_id', 'chain'], keep='last').to_csv(mt_path, index=False)
    pd.DataFrame([{'sequence_id': value, 'sequence': key, 'length': len(key)} for key, value in registry.items()]).to_csv(
        folder / 'sequences.csv', index=False)
    pd.DataFrame(coverage).to_csv(folder / 'coverage.csv', index=False)
    (folder / 'requested_sequence_ids.json').write_text(json.dumps(sorted(requested)) + '\n')
    # Copy only
    # requested cached chains; keep the full registry to preserve old IDs.
    for identifier in requested:
        old = reuse_root / 'sequences/embeddings' / f'{identifier}.npz'
        destination = folder / 'embeddings' / old.name
        if old.exists() and not destination.exists():
            shutil.copyfile(old, destination)
    return pd.DataFrame(coverage)


def prepare_embeddings(root, reuse_root, *, device='cuda'):
    coverage = prepare_sequences(root, reuse_root)
    folder = Path(root) / 'sequences'
    requested = json.loads((folder / 'requested_sequence_ids.json').read_text())
    report = extract_embeddings(folder, device=device, model_name='esmc_600m', sequence_ids=requested)
    report['n_mapping_success'] = int(coverage.status.eq('ok').sum())
    return report


def pool_features(root, task):
    """Save only successful vectors, plus the full declared coverage denominator."""
    import torch
    from threadpoolctl import threadpool_limits

    from topoppi_affinity.bilinear_features import bilinear_features
    from topoppi_affinity.embedding_baselines import pooled_features, zero_surface_slots
    from topoppi_affinity.graph_dataset import GraphDataset

    root = Path(root)
    target = 'pkd' if task == 'absolute' else 'ddg'
    frame = pd.read_csv(root / 'inputs' / f'{task}.csv', keep_default_na=False)
    mapping = pd.read_csv(root / 'sequences/coverage.csv', keep_default_na=False).set_index('sample_id')
    dataset = GraphDataset(frame.assign(**{target: 0.}), root / 'graphs/manifest.csv', root / 'sequences',
                           max_nodes=None, target=target)
    folder = root / 'features' / task
    folder.mkdir(parents=True, exist_ok=True)
    vectors, ids, pairs, pair_ids, coverage = [], [], [], [], []
    names, pair_names = None, None
    projection = np.random.default_rng(11).normal(size=(1152, 32)).astype(np.float32) / np.sqrt(32.)
    torch.set_num_threads(1)
    with threadpool_limits(limits=1):
        for index, row in enumerate(frame.itertuples()):
            record = {'sample_id': row.sample_id}
            try:
                if row.sample_id not in mapping.index or mapping.loc[row.sample_id, 'status'] != 'ok':
                    raise ValueError('Sequence mapping failed; mutant representation is not available')
                data = dataset[index]
                vector, current_names = pooled_features(data, target=target, use_surface=True)
                if not np.isfinite(vector).all():
                    raise ValueError('Nonfinite pooled input')
                if names is not None and not np.array_equal(names, current_names):
                    raise ValueError('Pooled feature columns changed')
                names = current_names
                vectors.append(vector)
                ids.append(row.sample_id)
                record.update(status='ok', error='')
                if task == 'absolute':
                    try:
                        pair, pair_names, _ = bilinear_features(data, projection, target=target)
                        if not np.isfinite(pair).all():
                            raise ValueError('Nonfinite bilinear input')
                        pairs.append(pair)
                        pair_ids.append(row.sample_id)
                        record.update(bilinear_status='ok', bilinear_error='')
                    except (ValueError, KeyError, IndexError) as error:
                        record.update(bilinear_status='error', bilinear_error=str(error))
            except (ValueError, KeyError, OSError, IndexError) as error:
                record.update(status='error', error=f'{type(error).__name__}: {error}')
            coverage.append(record)
            if (index + 1) % 100 == 0:
                print(f'Pool {task}: {index + 1}/{len(frame)}', flush=True)
    pd.DataFrame(coverage).to_csv(folder / 'pooling_coverage.csv', index=False)
    if not vectors:
        raise ValueError(f'No successful {task} pooled vectors; see pooling_coverage.csv')
    matrix, sample_ids = np.asarray(vectors), np.asarray(ids, dtype=str)
    for representation, values in [('plm_surface', matrix), ('plm_contact', zero_surface_slots(matrix, names))]:
        np.savez_compressed(folder / f'{target}_{representation}.npz', X=values,
                            sample_ids=sample_ids, feature_names=names)
    if pairs:
        np.savez_compressed(folder / 'pkd_bilinear_pair_surface.npz', X=np.asarray(pairs),
                            sample_ids=np.asarray(pair_ids, dtype=str), feature_names=pair_names)
    return pd.DataFrame(coverage)


def analytic_features(root):
    """Refine mutation graphs separately; keep each original contact slot intact."""
    from topoppi_affinity.analytic_curvature import replace_graph_curvature, update_pooled_curvature
    from topoppi_affinity.embedding_baselines import zero_surface_slots

    root = Path(root)
    folder = root / 'features/mutation'
    frame = pd.read_csv(root / 'inputs/mutation.csv', keep_default_na=False).set_index('sample_id')
    graphs = pd.read_csv(root / 'graphs/manifest.csv', keep_default_na=False)
    lookup = {(row.source, row.complex_id): row for row in graphs.itertuples() if row.status == 'ok'}
    with np.load(folder / 'ddg_plm_surface.npz', allow_pickle=False) as values:
        matrix, ids, names = values['X'], values['sample_ids'].astype(str), values['feature_names']
    indices = {sample: index for index, sample in enumerate(ids)}
    selected = frame.loc[ids].reset_index().sort_values(['source', 'complex_id', 'sample_id'])
    vectors, success, coverage, graph_rows = [], [], [], []
    for key, group in selected.groupby(['source', 'complex_id'], sort=True):
        record = lookup[key]
        output = root / 'graphs_analytic' / record.source / Path(record.graph_path).name
        try:
            if output.exists():
                graph = load_graph(output)
            else:
                graph = replace_graph_curvature(load_graph(record.graph_path), record.structure_path)
                save_graph(graph, output)
            graph_rows.append({'source': key[0], 'complex_id': key[1], 'graph_path': str(output), 'status': 'ok'})
        except (ValueError, KeyError, OSError, IndexError) as error:
            coverage.extend({'sample_id': sample, 'status': 'error', 'error': str(error)} for sample in group.sample_id)
            continue
        for row in group.itertuples():
            try:
                vector = update_pooled_curvature(matrix[indices[row.sample_id]], names, graph,
                                                 target='ddg', mutations=row.mutations)
                if not np.isfinite(vector).all():
                    raise ValueError('Nonfinite analytic vector')
                vectors.append(vector)
                success.append(row.sample_id)
                coverage.append({'sample_id': row.sample_id, 'status': 'ok', 'error': ''})
            except (ValueError, KeyError, IndexError) as error:
                coverage.append({'sample_id': row.sample_id, 'status': 'error', 'error': str(error)})
    pd.DataFrame(coverage).to_csv(folder / 'analytic_coverage.csv', index=False)
    if not vectors:
        raise ValueError('No successful analytic mutation vectors; see analytic_coverage.csv')
    values = np.asarray(vectors)
    original = matrix[[indices[sample] for sample in success]]
    if not np.array_equal(zero_surface_slots(values, names), zero_surface_slots(original, names)):
        raise ValueError('Analytic conversion changed the contact control')
    np.savez_compressed(folder / 'ddg_plm_analytic.npz', X=values, sample_ids=np.asarray(success, dtype=str),
                        feature_names=names, curvature_protocol=np.asarray('analytic_Gaussian_Goldman'))
    pd.DataFrame(graph_rows).to_csv(root / 'graphs_analytic/manifest.csv', index=False)
    return pd.DataFrame(coverage)
