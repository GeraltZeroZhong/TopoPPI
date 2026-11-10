import json

import numpy as np
import pandas as pd

from topoppi_affinity.paired_embeddings import PairedResidueFeatures


def test_pair_to_native_residue_mapping_averages_opposite_chains(tmp_path):
    sequences = tmp_path / 'sequences'
    sequences.mkdir()
    output = tmp_path / 'paired'
    embeddings = output / 'embeddings'
    embeddings.mkdir(parents=True)
    rows = []
    for chain, keys in [('A', ['A:10:', 'A:12:']), ('B', ['B:4:']), ('C', ['C:7:'])]:
        rows.append({'source': 'skempi', 'complex_id': 'c', 'chain': chain, 'residue_keys': json.dumps(keys)})
    pd.DataFrame(rows).to_csv(sequences / 'wild_type_chains.csv', index=False)
    pd.DataFrame([
        {'source': 'skempi', 'complex_id': 'c', 'state': 'WT', 'chain_a': 'A', 'chain_b': 'B',
         'a_is_first': True, 'pair_id': 'p1'},
        {'source': 'skempi', 'complex_id': 'c', 'state': 'WT', 'chain_a': 'A', 'chain_b': 'C',
         'a_is_first': False, 'pair_id': 'p2'},
    ]).to_csv(output / 'assignments.csv', index=False)
    np.savez(embeddings / 'p1.npz', embedding_a=np.full((2, 1280), 2.), embedding_b=np.full((1, 1280), 4.))
    np.savez(embeddings / 'p2.npz', embedding_a=np.full((1, 1280), 8.), embedding_b=np.full((2, 1280), 6.))
    graph = {'x': np.zeros((4, 45)), 'chain': np.array(['C', 'A', 'B', 'A']),
             'residue_key': np.array(['C:7:', 'A:12:', 'B:4:', 'A:10:'])}
    actual = PairedResidueFeatures(tmp_path, output).features('skempi', 'c', 'WT', graph)
    np.testing.assert_array_equal(actual[:, 0], [8., 4., 4., 4.])


