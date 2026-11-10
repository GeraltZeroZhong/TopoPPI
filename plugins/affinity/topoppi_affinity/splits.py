"""Reproducible experimental-system and partner-family splits for affinity data."""
from __future__ import annotations

import json
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


class DisjointSet:
    def __init__(self, values):
        self.parent = {value: value for value in values}

    def find(self, value):
        if self.parent[value] != value:
            self.parent[value] = self.find(self.parent[value])
        return self.parent[value]

    def union(self, left, right):
        left, right = self.find(left), self.find(right)
        if left != right:
            small, large = sorted((left, right))
            self.parent[large] = small


def group_split(groups, seed: int = 42, fractions=(0.7, 0.15, 0.15)) -> dict[str, str]:
    """Assign groups, without labels, to train/validation/test; large groups first.

    Counts are balanced approximately: no component is broken to meet a quota.
    At least three groups are required for an interpretable three-way split.
    """
    sizes = pd.Series(list(groups), dtype=str).value_counts()
    if len(sizes) < 3:
        raise ValueError('A three-way split requires at least three independent groups.')
    rng = np.random.default_rng(seed)
    keys = list(sizes.index.sort_values())
    rng.shuffle(keys)
    keys.sort(key=lambda value: -sizes[value])
    targets = np.asarray(fractions) * sizes.sum()
    counts = np.zeros(3)
    names = ('train', 'validation', 'test')
    assignment = {}
    for index, key in enumerate(keys):
        empty = np.flatnonzero(counts == 0)
        remaining = len(keys) - index
        if len(empty) == remaining:
            destination = int(empty[0])
        else:
            destination = int(np.argmax((targets - counts) / targets))
        assignment[str(key)] = names[destination]
        counts[destination] += sizes[key]
    return assignment


def within_complex_split(frame: pd.DataFrame, seed: int = 42) -> pd.Series:
    """Hold out mutation identities within each known complex.

    All temperatures/references for one mutation stay together. Complexes with
    fewer than five distinct variants remain train-only, explicitly reported.
    """
    result = pd.Series('train', index=frame.index, dtype=str)
    rng = np.random.default_rng(seed)
    for _, subset in frame.groupby('complex_id', sort=True):
        variants = sorted(subset['mutations'].fillna('').unique())
        if len(variants) < 5:
            continue
        rng.shuffle(variants)
        count = max(1, round(0.15 * len(variants)))
        assignments = {value: 'test' for value in variants[:count]}
        assignments.update({value: 'validation' for value in variants[count:2 * count]})
        result.loc[subset.index] = subset['mutations'].map(assignments).fillna('train')
    return result


def mmseqs_sequence_groups(frame: pd.DataFrame, output_dir: Path,
                           mmseqs: str = 'mmseqs', identity: float = 0.30,
                           coverage: float = 0.80, threads: int = 8) -> pd.Series:
    """Join systems sharing either partner family, including swapped partners.

    Use all discovered MMseqs sequence edges, rather than only representatives.
    Each chain is searched separately: concatenated antibody chains are not
    mistaken for one sequence. The split is an operational MMseqs 30%/80%
    homology test, not a guarantee that remote homology has been eliminated.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    sequences = sorted({sequence for column in ('sequence_a', 'sequence_b')
                        for group in frame[column] for sequence in str(group).split('|') if sequence})
    sequence_to_id = {sequence: f'seq{i:06d}' for i, sequence in enumerate(sequences)}
    fasta = output_dir / 'partner_sequences.fasta'
    fasta.write_text(''.join(f'>{identifier}\n{sequence}\n' for sequence, identifier in sequence_to_id.items()))
    hits = output_dir / 'sequence_hits.tsv'
    # Keep MMseqs scratch on the system temporary filesystem, which must
    # support the symbolic links used by MMseqs.
    scratch = tempfile.TemporaryDirectory(prefix='topoppi_affinity_mmseqs_')
    command = [str(mmseqs), 'easy-search', str(fasta), str(fasta), str(hits),
               scratch.name, '--min-seq-id', str(identity),
               '-c', str(coverage), '--cov-mode', '0', '--seq-id-mode', '0',
               '--alignment-mode', '3', '--max-seqs', str(len(sequences)),
               '-s', '7.5', '-e', '0.001', '--threads', str(threads),
               '--format-output', 'query,target,fident,qcov,tcov,evalue']
    try:
        with (output_dir / 'mmseqs.log').open('w') as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    finally:
        scratch.cleanup()
    ds = DisjointSet(sequence_to_id.values())
    with hits.open() as handle:
        for line in handle:
            left, right, similarity, qcov, tcov, _ = line.rstrip().split('\t')
            if float(similarity) >= identity and min(float(qcov), float(tcov)) >= coverage:
                ds.union(left, right)
    row_representatives = {}
    # Joining all partners of a complex creates the strict unseen-partner
    # connected component; it can be much broader than an unseen-PDB split.
    for index, row in frame.iterrows():
        members = [sequence_to_id[sequence] for column in ('sequence_a', 'sequence_b')
                   for sequence in row[column].split('|')]
        for member in members[1:]:
            ds.union(members[0], member)
        row_representatives[index] = members[0]
    # Different constructs from the same deposited PDB must also remain together.
    for _, subset in frame.groupby('pdb_id'):
        members = [row_representatives[index] for index in subset.index]
        for member in members[1:]:
            ds.union(members[0], member)
    roots = sorted({ds.find(value) for value in row_representatives.values()})
    root_to_group = {value: f'family_{i:04d}' for i, value in enumerate(roots)}
    result = pd.Series({index: root_to_group[ds.find(value)] for index, value in row_representatives.items()})
    version = subprocess.check_output([str(mmseqs), 'version'], text=True).strip()
    (output_dir / 'method.json').write_text(json.dumps({'command': command, 'version': version,
        'identity': identity, 'coverage_both_chains': coverage, 'evalue': 0.001,
        'sequence_count': len(sequences), 'family_count': len(roots)}, indent=2) + '\n')
    return result



def exact_interaction_groups(frame: pd.DataFrame) -> pd.Series:
    """Group same PDB and exact unordered full partner sequence sets."""
    ds = DisjointSet(frame.index)
    seen_pdb, seen_partners = {}, {}
    for index, row in frame.iterrows():
        partners = tuple(sorted(tuple(sorted(row[column].split('|')))
                                for column in ('sequence_a', 'sequence_b')))
        for value, seen in [(row['pdb_id'], seen_pdb), (partners, seen_partners)]:
            if value in seen:
                ds.union(index, seen[value])
            else:
                seen[value] = index
    roots = sorted({ds.find(index) for index in frame.index})
    names = {root: f'interaction_{i:04d}' for i, root in enumerate(roots)}
    return pd.Series({index: names[ds.find(index)] for index in frame.index})

def audit_splits(frame: pd.DataFrame) -> dict:
    """Fail on actual leakage that would invalidate held-out estimates."""
    report = {}
    for split, group in [('split_complex', 'pdb_id'), ('split_sequence', 'sequence_group')]:
        leaking = frame.groupby(group)[split].nunique()
        if (leaking > 1).any():
            raise ValueError(f'{split}: groups occur across partitions: {leaking[leaking > 1].index.tolist()}')
        report[split] = {str(key): int(value) for key, value in frame[split].value_counts().items()}
    exact_sequences = defaultdict(set)
    if {'sequence_a', 'sequence_b'}.issubset(frame.columns):
        for _, row in frame.iterrows():
            for column in ('sequence_a', 'sequence_b'):
                for sequence in row[column].split('|'):
                    exact_sequences[sequence].add(row['split_sequence'])
        if any(len(partitions) > 1 for partitions in exact_sequences.values()):
            raise ValueError('Identical partner sequence crosses sequence partitions.')
    report['exact_partner_sequence_leaks'] = 0
    if 'interaction_group' in frame and (frame.groupby('interaction_group')['split_complex'].nunique() > 1).any():
        raise ValueError('Exact interaction group crosses complex partitions.')
    mutations = frame[frame['n_mutations'] > 0]
    if len(mutations) and (mutations.groupby(['complex_id', 'mutations'])['split_within'].nunique() > 1).any():
        raise ValueError('Same experimental variant crosses within-complex partitions.')
    report['split_within'] = {str(key): int(value) for key, value in mutations['split_within'].value_counts().items()}
    return report


def add_splits(frame: pd.DataFrame, output_dir: Path, mmseqs: str = 'mmseqs', seed: int = 42) -> pd.DataFrame:
    frame = frame.reset_index(drop=True).copy()
    frame['sequence_group'] = mmseqs_sequence_groups(frame, output_dir, mmseqs)
    # Balance unique complexes, not repeated labels from large mutation series.
    systems = frame.drop_duplicates('complex_id')
    frame['interaction_group'] = exact_interaction_groups(frame)
    systems = frame.drop_duplicates('complex_id')
    complex_mapping = group_split(systems['interaction_group'], seed)
    sequence_mapping = group_split(systems['sequence_group'], seed)
    frame['split_complex'] = frame['interaction_group'].map(complex_mapping)
    frame['split_sequence'] = frame['sequence_group'].map(sequence_mapping)
    frame['split_within'] = within_complex_split(frame, seed)
    report = audit_splits(frame)
    (output_dir / 'split_audit.json').write_text(json.dumps(report, indent=2) + '\n')
    return frame
