"""Group-preserving partitions and prior-use metadata for supplied observations."""
from __future__ import annotations

import json
from collections import defaultdict

import numpy as np
import pandas as pd

from .data import parse_mutations, read_structure
from .splits import DisjointSet

TASKS = {
    'absolute_interaction': ('absolute', 'fold_interaction'),
    'absolute_family': ('absolute', 'fold_family'),
    'mutation_variant': ('mutation', 'fold_variant'),
    'mutation_interaction': ('mutation', 'fold_interaction'),
    'mutation_family': ('mutation', 'fold_family'),
}
IDENTITY_COLUMNS = ['complex_id', 'sequence_a', 'sequence_b']
INPUT_COLUMNS = [
    'sample_id', 'source', 'complex_id', 'pdb_id', 'structure_path',
    'source_structure_path', 'chains_a', 'chains_b', 'mutations', 'n_mutations',
    'mutation_numbering', 'structure_variant', 'sequence_a', 'sequence_b',
    'mutant_sequence_a', 'mutant_sequence_b', 'temperature_K', 'temperature_status',
    'family_group', 'interaction_group', 'initial_interaction_group', 'variant_group',
    'fold_interaction', 'fold_family', 'fold_variant', 'complex_seen_in_training',
    'exposure_roles', 'exposure_manifests', 'ever_used_row', 'ever_used_variant',
    'ever_used_interaction', 'ever_used_family', 'fresh_row', 'fresh_label',
    'fresh_variant', 'fresh_interaction', 'fresh_family',
]


def attach_global_groups(frame: pd.DataFrame, groups: pd.DataFrame) -> pd.DataFrame:
    """Match actual source structure sequences, not ambiguous complex IDs alone."""
    columns = IDENTITY_COLUMNS + ['family_group', 'initial_interaction_group']
    table = groups[columns].drop_duplicates()
    if table.duplicated(IDENTITY_COLUMNS).any():
        raise ValueError('Global group table has conflicting entries for the same structure sequences.')
    result = frame.merge(table, on=IDENTITY_COLUMNS, how='left', validate='many_to_one', sort=False)
    missing = result['family_group'].isna() | result['initial_interaction_group'].isna()
    if missing.any():
        raise ValueError(f'Eligible systems missing from global groups: {result.loc[missing, "sample_id"].tolist()}')
    result['interaction_group'] = result['initial_interaction_group']
    for column in ['family_group', 'interaction_group']:
        if result.groupby('pdb_id')[column].nunique().gt(1).any():
            raise ValueError(f'Same deposited PDB crosses {column} in the supplied global groups.')
    sequence_families = defaultdict(set)
    for row in result.drop_duplicates(IDENTITY_COLUMNS).itertuples():
        for sequence in (row.sequence_a + '|' + row.sequence_b).split('|'):
            sequence_families[sequence].add(row.family_group)
    if any(len(families) > 1 for families in sequence_families.values()):
        raise ValueError('An identical partner sequence crosses global family groups.')
    return result


def canonical_state_pair(wt: dict[str, str], mt: dict[str, str], sides: list[list[str]]) -> tuple:
    """Unordered partner/chain sets, with the entire WT/MT edge undirected.

    Chain multiplicity and partner membership are retained. WT and MT strings
    remain paired per physical chain, preventing a chain-specific change from
    being mistaken for a change in another sequence.
    """
    forward = tuple(sorted(tuple(sorted((wt[chain], mt[chain]) for chain in side)) for side in sides))
    reverse = tuple(sorted(tuple(sorted((mt[chain], wt[chain]) for chain in side)) for side in sides))
    return min(forward, reverse)


def physical_variant_groups(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Resolve actual PDB positions and group repeated physical mutations.

    This reads residue identities only. It never creates mutant coordinates or
    looks at Kd, temperature, experimental reference, or an effect label.
    """
    result = frame.copy()
    keys, mutant_a, mutant_b = [], [], []
    positions = {}
    for row in result.itertuples():
        sequences, identities = read_structure(str(row.structure_path))
        if row.structure_path not in positions:
            order = defaultdict(dict)
            for chain, number, insertion in identities:
                order[chain][(number, insertion)] = len(order[chain])
            positions[row.structure_path] = order
        sides = [row.chains_a.split(','), row.chains_b.split(',')]
        chains = sides[0] + sides[1]
        wt = {chain: sequences[chain] for chain in chains}
        actual_sides = ['|'.join(wt[chain] for chain in side) for side in sides]
        if actual_sides != [row.sequence_a, row.sequence_b]:
            raise ValueError(f'{row.sample_id}: structure sequences changed since eligible preparation.')
        mt = {chain: list(sequence) for chain, sequence in wt.items()}
        for original, chain, number, insertion, replacement in parse_mutations(row.mutations):
            if chain not in mt or (number, insertion) not in positions[row.structure_path][chain]:
                raise ValueError(f'{row.sample_id}: mutation not mapped to an observed partner residue.')
            index = positions[row.structure_path][chain][(number, insertion)]
            if index >= len(mt[chain]) or mt[chain][index] != original:
                raise ValueError(f'{row.sample_id}: mutation WT identity does not match structure sequence.')
            mt[chain][index] = replacement
        mt = {chain: ''.join(sequence) for chain, sequence in mt.items()}
        keys.append(canonical_state_pair(wt, mt, sides))
        mutant_a.append('|'.join(mt[chain] for chain in sides[0]))
        mutant_b.append('|'.join(mt[chain] for chain in sides[1]))
    names = {key: f'variant_{index:06d}' for index, key in enumerate(sorted(set(keys)))}
    result['variant_group'] = [names[key] for key in keys]
    result['mutant_sequence_a'], result['mutant_sequence_b'] = mutant_a, mutant_b
    identities = pd.DataFrame([
        {'variant_group': name, 'canonical_state_pair': json.dumps(key, separators=(',', ':'))}
        for key, name in names.items()
    ])
    return result, identities


def balanced_group_folds(frame: pd.DataFrame, group: str, *, n_folds: int = 5,
                         seed: int = 20260923, joint_systems: bool = True) -> dict[str, int]:
    """Label-blind greedy group allocation, balancing both task system counts.

    No connected component is split to satisfy quotas. Mutation-variant folds
    instead balance record counts, keeping every physical variant together.
    """
    if n_folds < 2:
        raise ValueError('At least two folds are required.')
    units = frame.drop_duplicates(['source', 'complex_id']) if joint_systems else frame
    sizes = units.groupby([group, 'benchmark_target']).size().unstack(fill_value=0).sort_index()
    if len(sizes) < n_folds:
        raise ValueError(f'{group} has fewer independent groups than requested folds.')
    totals = sizes.sum(axis=0).to_numpy(dtype=float)
    normalized = sizes.to_numpy(dtype=float) / totals
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(sizes)).tolist()
    order.sort(key=lambda index: -float(np.square(normalized[index]).sum()))
    counts = np.zeros((n_folds, normalized.shape[1]))
    assignment = {}
    for index in order:
        weight = normalized[index]
        costs = (np.square(counts + weight) - np.square(counts)).sum(axis=1)
        choices = np.flatnonzero(np.isclose(costs, costs.min(), rtol=0, atol=1e-15))
        fold = int(rng.choice(choices))
        assignment[str(sizes.index[index])] = fold
        counts[fold] += weight
    return assignment


def merge_variant_interactions(frame: pd.DataFrame, *, n_folds: int = 5,
                               seed: int = 20260923) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Keep identical physical edges together in unseen-interaction evaluation.

    Reciprocal mutations can connect different WT sequence pairs and therefore
    different legacy exact-interaction groups. Merge those components jointly
    with absolute rows; keep the supplied family components unchanged.
    """
    result = frame.copy()
    original_folds = balanced_group_folds(result, 'interaction_group', n_folds=n_folds, seed=seed)
    mutation = result[result.benchmark_target.eq('mutation')]
    sets = mutation.groupby('variant_group').interaction_group.nunique()
    bridge_groups = set(sets[sets.gt(1)].index)
    bridge = mutation[mutation.variant_group.isin(bridge_groups)][[
        'sample_id', 'source', 'complex_id', 'pdb_id', 'variant_group', 'family_group',
        'initial_interaction_group',
    ]].copy()
    bridge['fold_before_variant_union'] = bridge.initial_interaction_group.map(original_folds)
    bridge['crossed_fold_before_union'] = bridge.groupby('variant_group').fold_before_variant_union.transform('nunique').gt(1)
    ds = DisjointSet(sorted(result.interaction_group.unique()))
    for _, subset in mutation[mutation.variant_group.isin(bridge_groups)].groupby('variant_group'):
        if subset.family_group.nunique() != 1:
            raise ValueError('A canonical variant crosses global families; supplied global grouping is inconsistent.')
        groups = sorted(subset.interaction_group.unique())
        for group in groups[1:]:
            ds.union(groups[0], group)
    roots = sorted({ds.find(group) for group in result.interaction_group})
    names = {root: f'benchmark_interaction_{index:04d}' for index, root in enumerate(roots)}
    result['interaction_group'] = result.interaction_group.map(lambda group: names[ds.find(group)])
    bridge['interaction_group'] = bridge.sample_id.map(result.set_index('sample_id').interaction_group)
    return result, bridge

def add_exposure(frame: pd.DataFrame, exposure_frames: dict[str, tuple[str, pd.DataFrame]]) -> pd.DataFrame:
    """Mark prior use jointly across targets using caller-supplied exposure tables."""
    result = frame.copy()
    sample_lookup = result.set_index('sample_id')
    roles, manifests = defaultdict(set), defaultdict(set)
    for filename, (role, used) in exposure_frames.items():
        unknown = set(used.sample_id) - set(sample_lookup.index)
        if unknown:
            raise ValueError(f'{filename}: historically used samples absent from eligible pool: {sorted(unknown)}')
        for row in used.itertuples():
            expected = sample_lookup.loc[row.sample_id]
            if (row.source, row.complex_id) != (expected.source, expected.complex_id):
                raise ValueError(f'{filename}: changed source/system for {row.sample_id}.')
            roles[row.sample_id].add(role)
            manifests[row.sample_id].add(filename)
    result['exposure_roles'] = result.sample_id.map(lambda sample: '|'.join(sorted(roles.get(sample, ()))))
    result['exposure_manifests'] = result.sample_id.map(lambda sample: '|'.join(sorted(manifests.get(sample, ()))))
    result['ever_used_row'] = result.sample_id.isin(roles)
    used = result[result.ever_used_row]
    for suffix, column in [('interaction', 'interaction_group'), ('family', 'family_group')]:
        result[f'ever_used_{suffix}'] = result[column].isin(used[column])
        result[f'fresh_{suffix}'] = ~result[f'ever_used_{suffix}']
    result['fresh_row'] = ~result.ever_used_row
    result['fresh_label'] = result.fresh_row
    # There is no mutation edge for the absolute task: variant fields remain
    # empty there rather than inventing a second kind of mutation identity.
    mutation = result.benchmark_target.eq('mutation')
    exposed_variants = set(used.loc[used.benchmark_target.eq('mutation'), 'variant_group'])
    result['ever_used_variant'] = pd.Series(pd.NA, index=result.index, dtype='boolean')
    result.loc[mutation, 'ever_used_variant'] = result.loc[mutation, 'variant_group'].isin(exposed_variants)
    result['fresh_variant'] = ~result.ever_used_variant
    return result


def make_splits(frame: pd.DataFrame, *, n_folds: int = 5, seed: int = 20260923
                ) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    result = frame.copy()
    for group, column in [('interaction_group', 'fold_interaction'), ('family_group', 'fold_family')]:
        mapping = balanced_group_folds(result, group, n_folds=n_folds, seed=seed)
        result[column] = result[group].map(mapping).astype(int)
    mutation = result.benchmark_target.eq('mutation')
    mapping = balanced_group_folds(result[mutation], 'variant_group', n_folds=n_folds,
                                   seed=seed, joint_systems=False)
    result['fold_variant'] = pd.Series(pd.NA, index=result.index, dtype='Int64')
    result.loc[mutation, 'fold_variant'] = result.loc[mutation, 'variant_group'].map(mapping)
    result['complex_seen_in_training'] = pd.Series(pd.NA, index=result.index, dtype='boolean')
    mutation_frame = result[mutation]
    seen = mutation_frame.groupby(['source', 'complex_id']).fold_variant.transform('nunique').gt(1)
    result.loc[mutation, 'complex_seen_in_training'] = seen
    splits = {}
    shared = ['sample_id', 'source', 'complex_id', 'pdb_id', 'n_mutations', 'family_group',
              'interaction_group', 'initial_interaction_group', 'variant_group', 'exposure_roles',
              'exposure_manifests', 'fresh_row', 'fresh_label', 'fresh_variant',
              'fresh_interaction', 'fresh_family']
    for task, (target, column) in TASKS.items():
        subset = result[result.benchmark_target.eq(target)]
        output = subset[shared].copy()
        output['fold'] = subset[column].astype(int)
        if task == 'mutation_variant':
            output['complex_seen_in_training'] = subset.complex_seen_in_training
            output['evaluation_eligible'] = subset.complex_seen_in_training
            freshness = subset.fresh_variant
        else:
            output['complex_seen_in_training'] = False
            output['evaluation_eligible'] = True
            freshness = subset['fresh_family' if task.endswith('_family') else 'fresh_interaction']
        output['fresh_only'] = output.evaluation_eligible & subset.fresh_row & freshness
        splits[task] = output.reset_index(drop=True)
    return result, splits
