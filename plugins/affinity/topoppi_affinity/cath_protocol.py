"""Map supplied eligible mutations to released CATH-ddG labels and partitions."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


def prepare_cath_protocol(manifest, raw, author, output):
    frame = pd.read_csv(manifest, keep_default_na=False)
    original = pd.read_csv(raw, sep=';', keep_default_na=False)
    released = pd.read_csv(author, keep_default_na=False)
    # The prefix is the author's raw SKEMPI zero-based row number, not a PDB ID.
    indices = released['#Pdb'].str.split('_').str[0].astype(int)
    if indices.duplicated().any():
        raise ValueError('Author row indices are not unique')
    corresponding = original.iloc[indices].reset_index(drop=True)
    if not (corresponding['#Pdb'].eq(released.complex_PPI).all()
            and corresponding['Mutation(s)_cleaned'].eq(released['Mutation(s)_cleaned']).all()):
        raise ValueError('Author row numbering does not match this raw SKEMPI release')
    released['sample_id'] = [f'skempi2:{i+1:05d}' for i in indices]
    released['partition'] = released.cath_fold.map({'train': 'train', 'val': 'test', '': 'unassigned'})
    if released.partition.isna().any():
        raise ValueError('Unknown official CATH partition')
    source_columns = ['sample_id', 'partition', 'complex_PPI', 'ddG']
    mapped = frame.merge(released[source_columns], on='sample_id', how='left', validate='one_to_one')
    if mapped.partition.isna().any():
        raise ValueError('Eligible rows absent from the author release')
    mapped['ddg_at_reported_temperature'] = mapped.ddg
    mapped['ddg'] = pd.to_numeric(mapped.pop('ddG'), errors='coerce')
    if not np.isfinite(mapped.ddg).all():
        raise ValueError('Nonfinite author target in the eligible intersection')
    mapped['label_protocol'] = 'author_CATH_ddG; original experimental-temperature target retained separately'
    mapped['evaluation_scope'] = 'official_CATH_matched_eligible'
    inventory = released[['sample_id', '#Pdb', 'complex_PPI', 'partition']].copy()
    inventory['eligible'] = inventory.sample_id.isin(frame.sample_id)
    inventory['reason'] = np.where(inventory.eligible, '', 'outside the fixed eligible mutation pool')
    train = mapped[mapped.partition.eq('train')]
    test = mapped[mapped.partition.eq('test')]
    overlap = {key: len(set(train[key]) & set(test[key]))
               for key in ['variant_group', 'interaction_group', 'family_group']}
    def chains(part):
        return {seq for row in part.itertuples() for side in [row.sequence_a, row.sequence_b]
                for seq in side.split('|')}
    overlap['identical_observed_chain_sequences'] = len(chains(train) & chains(test))
    if any(overlap.values()):
        raise ValueError(f'Official matched holdout crosses available isolation groups: {overlap}')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    mapped.to_csv(output / 'manifest.csv', index=False)
    mapped[['sample_id', 'partition']].to_csv(output / 'partitions.csv', index=False)
    inventory.to_csv(output / 'author_coverage.csv', index=False)
    references = []
    # Published predictions, not newly executed local energy calculations.
    for column, method in [('ddG_FoldX', 'FoldX-author-table'), ('Flex-ddG', 'flex-ddG-author-table')]:
        values = released[['sample_id', column]].copy()
        values['predicted'] = pd.to_numeric(values[column], errors='coerce')
        values['method'] = method
        values['status'] = np.where(np.isfinite(values.predicted), 'ok', 'error')
        values['error'] = np.where(values.status.eq('ok'), '', 'Missing published scalar prediction')
        values['target'] = 'ddg'
        values['provenance'] = 'author supplied CATH table; not recomputed locally'
        values['comparison_scope'] = 'published_reference'
        references.append(values.drop(columns=column))
    pd.concat(references, ignore_index=True).to_csv(output / 'reference_predictions.csv', index=False)
    report = {
        'official_counts': released.partition.value_counts().to_dict(),
        'matched_counts': mapped.partition.value_counts().to_dict(),
        'official_test_complexes': int(released[released.partition.eq('test')].complex_PPI.nunique()),
        'matched_test_complexes': int(test.complex_id.nunique()),
        'available_overlap_checks': overlap,
        'domain_assignment_audit': 'not independently repeated; full author domain annotation absent',
        'partition_semantics': 'Author val is the published heldout TEST, never a tuning validation set',
        'target_semantics': 'Use released CATH-ddG targets; preserve SKEMPI targets calculated at the reported temperature',
        'target_mean_absolute_difference_from_reported_temperature': float(np.abs(mapped.ddg - mapped.ddg_at_reported_temperature).mean()),
        'source': str(Path(author).resolve()), 'source_url': 'https://github.com/ak422/CATH-ddG',
        'no_new_independent_confirmation': True, 'models_trained': False,
    }
    (output / 'preparation.json').write_text(json.dumps(report, indent=2) + '\n')
    return report
