"""The released Proteinea PPB split, with actual provided mutant sequences.

The filtered labels average Kd before taking its logarithm. Source metadata can
represent only the first contributing record. Neither its official test split
nor a new local split erases previous project or author-checkpoint exposure.
"""
from __future__ import annotations

import json
import shutil
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from topoppi_affinity.sequence_baselines import sequence_state_features


def _chains(value):
    return tuple(part.strip() for part in str(value).split(',') if part.strip())


def sequence_pair(first, second):
    return tuple(sorted((tuple(sorted(first)), tuple(sorted(second)))))


def _exact_components(frame):
    """Grouping for uncertainty, not a claim of homology isolation."""
    parents = list(range(len(frame)))

    def root(index):
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    seen = {}
    for index, row in enumerate(frame.itertuples()):
        for sequence in row.sequence_a.split('|') + row.sequence_b.split('|'):
            if sequence in seen:
                parents[root(index)] = root(seen[sequence])
            else:
                seen[sequence] = index
    identifiers = {}
    result = []
    for index in range(len(frame)):
        key = root(index)
        if key not in identifiers:
            identifiers[key] = f'ppb_exact_component_{len(identifiers):05d}'
        result.append(identifiers[key])
    return result


def prepare_ppb_protocol(source, output, *, benchmark_root=None):
    source, output = Path(source), Path(output)
    raw = pd.read_csv(source, keep_default_na=False)
    rows = []
    for index, row in enumerate(raw.to_dict('records')):
        first, second = _chains(row['Ligand Sequences']), _chains(row['Receptor Sequences'])
        ca, cb = _chains(row['Ligand Chains']), _chains(row['Receptor Chains'])
        if not first or not second or len(ca) != len(first) or len(cb) != len(second):
            raise ValueError(f'PPB row {index}: chains and provided sequences disagree')
        kd = float(row['KD(M)'])
        if not np.isfinite(kd) or kd <= 0:
            raise ValueError(f'PPB row {index}: invalid released Kd')
        pdb = str(row['PDB']).lower()
        partner_groups = sorted([','.join(sorted(ca)), ','.join(sorted(cb))])
        rows.append({
            'sample_id': f'ppb_filtered:{index + 1:05d}', 'source': 'PPB_Proteinea_filtered',
            'source_row_1based': index + 1, 'source_datasets': row['Source Data Set'],
            'complex_id': ':'.join([pdb, *partner_groups]), 'source_complex_id': row['Complex ID'], 'pdb_id': pdb,
            'chains_a': ','.join(ca), 'chains_b': ','.join(cb),
            'sequence_a': '|'.join(first), 'sequence_b': '|'.join(second),
            'mutations': row['Mutations'], 'kd_M': kd, 'pkd': -np.log10(kd),
            'partition': {'train': 'train', 'val': 'validation', 'test': 'test'}[row['split']],
            'reference': row['Affinity PubMed ID'], 'temperature_reported': row['Temperature(K)'],
            'label_protocol': 'released arithmetic mean Kd then -log10; first-record metadata',
        })
    frame = pd.DataFrame(rows)
    frame['family_group'] = _exact_components(frame)
    pairs, ids = {}, []
    for row in frame.itertuples():
        pair = sequence_pair(row.sequence_a.split('|'), row.sequence_b.split('|'))
        if pair not in pairs:
            pairs[pair] = f'ppb_exact_pair_{len(pairs):05d}'
        ids.append(pairs[pair])
    frame['interaction_group'] = ids
    train_sequences = {seq for row in frame[frame.partition.eq('train')].itertuples()
                       for side in (row.sequence_a, row.sequence_b) for seq in side.split('|')}
    validation_sequences = {seq for row in frame[frame.partition.eq('validation')].itertuples()
                            for side in (row.sequence_a, row.sequence_b) for seq in side.split('|')}
    frame['exact_chain_seen_in_official_train'] = [
        bool(set(row.sequence_a.split('|') + row.sequence_b.split('|')) & train_sequences)
        for row in frame.itertuples()]
    frame['exact_chain_seen_in_official_train_or_validation'] = [
        bool(set(row.sequence_a.split('|') + row.sequence_b.split('|')) & (train_sequences | validation_sequences))
        for row in frame.itertuples()]
    overlap = {}
    for first, second in [('train', 'validation'), ('train', 'test'), ('validation', 'test')]:
        a, b = frame[frame.partition.eq(first)], frame[frame.partition.eq(second)]
        seq_a = {s for row in a.itertuples() for side in (row.sequence_a, row.sequence_b) for s in side.split('|')}
        seq_b = {s for row in b.itertuples() for side in (row.sequence_a, row.sequence_b) for s in side.split('|')}
        overlap[f'{first}/{second}'] = {
            'exact_sequence_pairs': len(set(a.interaction_group) & set(b.interaction_group)),
            'identical_individual_chain_sequences': len(seq_a & seq_b),
        }
    output.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output / 'manifest.csv', index=False)
    frame[['sample_id', 'partition']].to_csv(output / 'partitions.csv', index=False)
    matched = pd.DataFrame()
    if benchmark_root is not None:
        matched = match_cached_structures(frame, Path(benchmark_root))
        matched.to_csv(output / 'structure_coverage.csv', index=False)
        selected = frame.merge(matched[matched.status.eq('ok')], on='sample_id', validate='one_to_one')
        selected.to_csv(output / 'structure_manifest.csv', index=False)
    report = {
        'source': str(source), 'source_url': 'https://huggingface.co/datasets/proteinea/ppb_affinity',
        'partitions': frame.partition.value_counts().to_dict(), 'label_protocol': rows[0]['label_protocol'],
        'sequence_protocol': 'Use released full recovered sequences as supplied; mutations already applied',
        'grouping_for_uncertainty': 'Exact individual-chain connected components; not CATH or sequence-homology groups',
        'complex_grouping': 'PDB and unordered partner chain sets; source record including mutation/reference retained separately',
        'exact_overlap': overlap, 'training_started': False,
        'structural_subset': 'Only cached same-PDB/same-partner and exact provided-sequence matches; all exclusions retained',
        'structure_matched': int(matched.status.eq('ok').sum()) if len(matched) else 0,
        'structure_failures': matched.loc[matched.status.ne('ok'), 'error'].value_counts().to_dict() if len(matched) else {},
        'structure_partition_counts': (frame[frame.sample_id.isin(matched.loc[matched.status.eq('ok'), 'sample_id'])]
                                       .partition.value_counts().to_dict()) if len(matched) else {},
        'independent_confirmation': False,
    }
    (output / 'preparation.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def match_cached_structures(frame, benchmark_root):
    """Reuse absolute cached inputs only when their entire molecular state matches.

    A sequence match alone cannot license reuse of an unrelated experimental
    structure. Mutation annotations alone cannot license use of wild-type graphs.
    """
    cached = pd.read_csv(benchmark_root / 'manifests/absolute.csv', keep_default_na=False)
    lookup = {}
    def key_for(row):
        # Chain-to-sequence assignment matters within a multi-chain partner.
        sides = [tuple(sorted(zip(_chains(getattr(row, f'chains_{side}')),
                                  getattr(row, f'sequence_{side}').split('|'), strict=True)))
                 for side in ('a', 'b')]
        return row.pdb_id.lower(), tuple(sorted(sides))

    for row in cached.itertuples():
        lookup.setdefault(key_for(row), []).append(row.sample_id)
    result = []
    for row in frame.itertuples():
        matches = lookup.get(key_for(row), [])
        result.append({'sample_id': row.sample_id, 'cached_sample_id': matches[0] if len(matches) == 1 else '',
                       'status': 'ok' if len(matches) == 1 else 'excluded',
                       'error': '' if len(matches) == 1 else 'No unique exact cached molecular-state match'})
    return pd.DataFrame(result)


def prepare_ppb_sequences(manifest, output, *, reuse_sequences=None):
    """Register actual released sequences and reuse embeddings by exact string."""
    frame = pd.read_csv(manifest, usecols=['sample_id', 'sequence_a', 'sequence_b'], keep_default_na=False)
    output = Path(output)
    (output / 'embeddings').mkdir(parents=True, exist_ok=True)
    reuse = {}
    if reuse_sequences is not None:
        reuse_sequences = Path(reuse_sequences)
        old = pd.read_csv(reuse_sequences / 'sequences.csv', keep_default_na=False)
        reuse = dict(zip(old.sequence, old.sequence_id, strict=True))
    registry, mapping, reused = {}, [], 0
    for row in frame.itertuples():
        for side in ('a', 'b'):
            for chain_index, sequence in enumerate(getattr(row, f'sequence_{side}').split('|')):
                if sequence not in registry:
                    sid = f'ppb_sequence_{len(registry):06d}'
                    registry[sequence] = sid
                    if sequence in reuse:
                        cached = reuse_sequences / 'embeddings' / f'{reuse[sequence]}.npz'
                        dest = output / 'embeddings' / f'{sid}.npz'
                        if cached.exists():
                            if not dest.exists():
                                shutil.copy2(cached, dest)
                            reused += 1
                mapping.append({'sample_id': row.sample_id, 'side': side, 'chain_index': chain_index,
                                'sequence_id': registry[sequence]})
    pd.DataFrame([{'sequence_id': sid, 'sequence': seq, 'length': len(seq)}
                  for seq, sid in registry.items()]).to_csv(output / 'sequences.csv', index=False)
    pd.DataFrame(mapping).to_csv(output / 'sample_chains.csv', index=False)
    report = {'unique_sequences': len(registry), 'reusable_v1_embeddings': reused,
              'sequence_scope': 'Released recovered actual sequences, including jointly applied mutations',
              'not_mutated_twice': True}
    (output / 'registration.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def pool_ppb_sequences(manifest, sequence_dir, output, *, limit=None):
    frame = pd.read_csv(manifest, usecols=['sample_id'], keep_default_na=False)
    if limit is not None:
        frame = frame.iloc[:limit]
    sequence_dir, output = Path(sequence_dir), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    registry = pd.read_csv(sequence_dir / 'sequences.csv', keep_default_na=False).set_index('sequence_id')
    mapping = pd.read_csv(sequence_dir / 'sample_chains.csv', keep_default_na=False)
    lookup = {(sid, side): group.sort_values('chain_index').sequence_id.tolist()
              for (sid, side), group in mapping.groupby(['sample_id', 'side'])}

    @lru_cache(maxsize=16)
    def read(sid):
        sequence = registry.loc[sid, 'sequence']
        with np.load(sequence_dir / 'embeddings' / f'{sid}.npz', allow_pickle=False) as saved:
            if str(saved['sequence']) != sequence or str(saved['model_name']) != 'esmc_600m':
                raise ValueError('Cached PPB sequence or model mismatch')
            values = saved['embedding'].astype(np.float32)
        if values.shape != (len(sequence), 1152) or not np.isfinite(values).all():
            raise ValueError('Invalid PPB residue embedding')
        return sequence, values

    vectors, ids, coverage, names = [], [], [], None
    for row in frame.itertuples():
        try:
            vector, names = sequence_state_features([[read(sid) for sid in lookup[(row.sample_id, side)]]
                                                     for side in ('a', 'b')])
            vectors.append(vector)
            ids.append(row.sample_id)
            coverage.append({'sample_id': row.sample_id, 'status': 'ok', 'error': ''})
        except (KeyError, ValueError, OSError) as error:
            coverage.append({'sample_id': row.sample_id, 'status': 'error', 'error': str(error)})
    pd.DataFrame(coverage).to_csv(output / 'coverage.csv', index=False)
    if not vectors:
        raise ValueError('No usable PPB embeddings; see coverage.csv')
    np.savez_compressed(output / 'pkd_sequence_global.npz', X=np.asarray(vectors),
                        sample_ids=np.asarray(ids, dtype=str), feature_names=names)
    return {'declared': len(frame), 'success': len(ids), 'failed': len(frame)-len(ids)}
