"""Partner-preserving controls for native-neighbor feature blocks."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .features import NEIGHBOR_FEATURES, neighbor_means
from .partner_features import QUERY_COLUMNS

CONTROL_KINDS = ('native_source_shuffle', 'native_append_shuffle')


def native_raw_columns():
    return ['topology.neighbor_count'] + [f'topology.mean.{source}' for source in NEIGHBOR_FEATURES]


def partner_permutation(frame, seed):
    """A joint row permutation of all WT nodes, separately on each partner."""
    permutation = np.arange(len(frame))
    for (complex_id, side), positions in frame.groupby(QUERY_COLUMNS, sort=True).indices.items():
        complex_bytes, side_bytes = str(complex_id).encode('utf-8'), str(side).encode('utf-8')
        # Length-delimited byte entropy gives each actual query its own stream.
        # Adding another complex to a batch cannot change this realization.
        entropy = [int(seed), len(complex_bytes), *complex_bytes, len(side_bytes), *side_bytes]
        rng = np.random.default_rng(np.random.SeedSequence(entropy))
        ordered = positions[np.argsort(frame.iloc[positions].residue_key.to_numpy())]
        permutation[ordered] = rng.permutation(ordered)
    return permutation


def _recompute_native_relative(frame):
    grouped = frame.groupby(QUERY_COLUMNS, sort=False)
    for column in native_raw_columns():
        values = grouped[column]
        mean, std = values.transform('mean'), values.transform('std', ddof=0)
        frame[f'relative.partner.zscore.{column}'] = (frame[column] - mean) / std.where(std > 0, 1.)
        frame[f'relative.partner.percentile.{column}'] = (values.rank(method='average') - .5) / values.transform('size')
        frame[f'context.partner.mean.{column}'] = mean
        frame[f'context.partner.deviation.{column}'] = frame[column] - mean
    return frame


def control_features(prepared, graphs, kind, seed):
    """Change only the appended native block, preserving focal/spatial inputs.

    Source shuffle preserves native counts but need not preserve readout
    marginals. Append shuffle preserves the full-WT joint readout distribution,
    not necessarily its measured or training subset distribution. Both controls
    use exactly the same full-WT partner permutations for a given seed.
    """
    if kind not in CONTROL_KINDS:
        raise ValueError(f'Unknown neighbor feature control: {kind}')
    result = prepared.copy().reset_index(drop=True)
    permutation = partner_permutation(result, seed)
    columns = native_raw_columns()
    if kind == 'native_append_shuffle':
        result.loc[:, columns] = result.iloc[permutation][columns].to_numpy(float)
    else:
        shuffled_source = result.iloc[permutation][NEIGHBOR_FEATURES].to_numpy(float)
        for cid, positions in result.groupby('complex_id', sort=True).indices.items():
            part = result.iloc[positions]
            graph = graphs[str(cid)]
            graph_keys = np.asarray(graph['residue_keys']).astype(str)
            lookup = pd.Series(np.arange(len(part)), index=part.residue_key)
            graph_to_local = lookup.reindex(graph_keys).to_numpy()
            original_edges = np.asarray(graph['native_neighbor_edges'], int).reshape(-1, 2)
            mapped = graph_to_local[original_edges]
            if not np.isfinite(mapped).all():
                raise ValueError('Native edges reference nodes outside the full WT interface context.')
            edges = mapped.astype(int)
            if len(edges) and not (part.side.to_numpy()[edges[:, 0]] == part.side.to_numpy()[edges[:, 1]]).all():
                raise ValueError('Native edges must remain within a partner query.')
            means, counts = neighbor_means(shuffled_source[positions], edges)
            if not np.array_equal(counts, part['topology.neighbor_count'].to_numpy()):
                raise ValueError('Cached native counts disagree with the fixed graph.')
            result.loc[positions, columns[1:]] = means
    return _recompute_native_relative(result)
