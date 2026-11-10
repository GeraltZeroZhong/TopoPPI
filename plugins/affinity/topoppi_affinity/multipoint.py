"""Joint multi-mutant predictions versus explicit sums of singleton predictions."""
from __future__ import annotations

from pathlib import Path

import pandas as pd


def prepare_singletons(manifest, output):
    frame = pd.read_csv(manifest, keep_default_na=False)
    records = []
    for row in frame[frame.n_mutations > 1].to_dict('records'):
        for index, mutation in enumerate(row['mutations'].split(',')):
            record = dict(row)
            record.update(sample_id=f"singleton:{row['sample_id']}:{index}", parent_sample_id=row['sample_id'],
                          mutations=mutation, mutations_pdb='', n_mutations=1,
                          purpose='unlabeled_singleton_inference_only')
            for column in ['ddg', 'pkd', 'pkd_wt', 'kd_M', 'kd_wt_M', 'kd_mut_M', 'dg_kcal_mol']:
                if column in record:
                    record[column] = ''
            records.append(record)
    result = pd.DataFrame(records)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)
    return {'n_parents': result.parent_sample_id.nunique(), 'n_singleton_inputs': len(result),
            'labels_used': False, 'purpose': 'prediction-time additive control, never added as training labels'}
