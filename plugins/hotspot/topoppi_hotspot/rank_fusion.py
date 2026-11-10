"""Combine ranking scores through within-partner candidate percentiles."""
from __future__ import annotations

import numpy as np
import pandas as pd
from Bio.Data.PDBData import protein_letters_3to1


def standard_nonalanine_candidates(frame):
    names = frame.residue_key.str.split(':').str[1]
    return frame.wt.isin(list('CDEFGHIKLMNPQRSTVWY')).to_numpy() & names.isin(protein_letters_3to1).to_numpy()


def partner_rank_percentiles(scores, queries, eligible):
    """Midrank percentiles over each full non-Ala candidate query; others NaN."""
    result = np.full(len(scores), np.nan)
    for query in np.unique(queries):
        indices = np.flatnonzero((queries == query) & eligible)
        if len(indices):
            rank = pd.Series(np.asarray(scores)[indices]).rank(method='average').to_numpy()
            result[indices] = (rank - .5) / len(indices)
    return result


def aggregate_ranks(score_columns, queries, eligible):
    return np.mean([partner_rank_percentiles(scores, queries, eligible) for scores in score_columns], axis=0)
