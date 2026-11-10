"""BindPred-style ESM2/CatBoost with explicit sequence pooling and local splits.

The released head is separate from the locally retrained baseline. The author
provides ordered 1280+1280 features but no complete multichain extractor; our
residue/chain pooling convention is therefore an explicit adaptation.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from topoppi_affinity.benchmark_models import LoadedFeatures
from topoppi_affinity.sequence_embeddings import window_starts

ESM2_NAME = 'esm2_t33_650M_UR50D'
METHOD = 'BindPred-ESM2-retrained'
PARAMETERS = dict(iterations=2000, learning_rate=.08, depth=4, loss_function='RMSE',
                  random_seed=11, thread_count=2, task_type='CPU', verbose=False,
                  allow_writing_files=False)
POOLING = ('Final-layer residue mean excluding BOS/EOS; mean over overlapping '
           '1000-residue windows at each residue; residue-count weighted chains; '
           'ordered manifest partner A then B; no contacts or surface inputs')


def register_sequences(manifests, output):
    """Deduplicate actual supplied sequences, including PPB's already mutated chains."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    frame = pd.concat([pd.read_csv(p, usecols=['sample_id', 'sequence_a', 'sequence_b'],
                                 keep_default_na=False) for p in manifests], ignore_index=True)
    frame = frame.drop_duplicates('sample_id')
    sequences, mapping = {}, []
    for row in frame.itertuples():
        for side in ('a', 'b'):
            for order, sequence in enumerate(getattr(row, f'sequence_{side}').split('|')):
                if not sequence:
                    raise ValueError(f'{row.sample_id}: empty partner sequence')
                if sequence not in sequences:
                    sequences[sequence] = f'esm2_{len(sequences):06d}'
                mapping.append(dict(sample_id=row.sample_id, side=side, chain_order=order,
                                    sequence_id=sequences[sequence], length=len(sequence)))
    pd.DataFrame([dict(sequence_id=sid, sequence=seq, length=len(seq))
                  for seq, sid in sequences.items()]).to_csv(output/'sequences.csv', index=False)
    pd.DataFrame(mapping).to_csv(output/'mapping.csv', index=False)
    return {'samples': len(frame), 'unique_sequences': len(sequences)}


def embed_sequences(root, checkpoint, *, device='cuda', limit=None):
    """Frozen ESM2-650M; cache means only, without fitting or reading affinity labels."""
    import torch
    from transformers import AutoModel, AutoTokenizer

    root = Path(root)
    table = pd.read_csv(root/'sequences.csv', keep_default_na=False)
    if limit is not None:
        table = table.iloc[:limit]
    output = root/'means'
    output.mkdir(exist_ok=True)
    pending = table[~table.sequence_id.map(lambda sid: (output/f'{sid}.npz').exists())]
    started = time.monotonic()
    if len(pending):
        torch.set_num_threads(2)
        tokenizer = AutoTokenizer.from_pretrained(str(checkpoint), local_files_only=True)
        model = AutoModel.from_pretrained(str(checkpoint), local_files_only=True,
                   torch_dtype=torch.float16 if device.startswith('cuda') else torch.float32,
                   low_cpu_mem_usage=True, add_pooling_layer=False).to(device).eval()
        if model.config.hidden_size != 1280 or model.config.num_hidden_layers != 33:
            raise ValueError('BindPred baseline requires ESM2-650M final-layer 1280-dimensional features')
        for index, row in enumerate(pending.itertuples(), 1):
            accumulated = np.zeros((len(row.sequence), 1280), dtype=np.float32)
            count = np.zeros((len(row.sequence), 1), dtype=np.float32)
            for start in window_starts(len(row.sequence)):
                sequence = row.sequence[start:start+1000]
                tokens = tokenizer(sequence, return_tensors='pt', add_special_tokens=True).to(device)
                with torch.inference_mode():
                    value = model(**tokens).last_hidden_state[0, 1:len(sequence)+1].float().cpu().numpy()
                accumulated[start:start+len(sequence)] += value
                count[start:start+len(sequence)] += 1
            mean = (accumulated/count).mean(axis=0)
            np.savez_compressed(output/f'{row.sequence_id}.npz', mean=mean, sequence=row.sequence,
                                model_name=ESM2_NAME)
            if index % 100 == 0:
                print(f'BindPred ESM2: {index}/{len(pending)} new chains, {time.monotonic()-started:.1f}s', flush=True)
    result = {'requested': len(table), 'newly_encoded': len(pending), 'checkpoint': str(checkpoint),
              'pooling': POOLING, 'seconds': time.monotonic()-started, 'labels_read': False}
    (root/'embedding_provenance.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


def pool_features(manifest, root, output):
    root, output = Path(root), Path(output)
    frame = pd.read_csv(manifest, usecols=['sample_id'], keep_default_na=False)
    registry = pd.read_csv(root/'sequences.csv', keep_default_na=False).set_index('sequence_id')
    mapping = pd.read_csv(root/'mapping.csv', keep_default_na=False)
    groups = {key: part for key, part in mapping.groupby(['sample_id', 'side'], sort=False)}
    values, rows, cache = [], [], {}
    for row in frame.itertuples():
        result = dict(sample_id=row.sample_id, status='ok', error='')
        try:
            partners = []
            for side in ('a', 'b'):
                chains = groups[(row.sample_id, side)].sort_values('chain_order')
                means = []
                for chain in chains.itertuples():
                    if chain.sequence_id not in cache:
                        with np.load(root/'means'/f'{chain.sequence_id}.npz', allow_pickle=False) as data:
                            if str(data['model_name']) != ESM2_NAME or str(data['sequence']) != registry.loc[chain.sequence_id, 'sequence']:
                                raise ValueError('Cached ESM2 model/sequence differs from registry')
                            cache[chain.sequence_id] = data['mean'].astype(np.float32)
                    means.append(cache[chain.sequence_id])
                partners.append(np.average(means, axis=0, weights=chains.length))
            value = np.concatenate(partners).astype(np.float32)
            if value.shape != (2560,) or not np.isfinite(value).all():
                raise ValueError('Invalid pooled ESM2 pair')
        except Exception as error:
            value = np.full(2560, np.nan, dtype=np.float32)
            result.update(status='error', error=f'{type(error).__name__}: {error}')
        values.append(value)
        rows.append(result)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, X=np.stack(values), sample_ids=frame.sample_id.to_numpy(str),
                        feature_names=np.asarray([f'esm2.{side}.{i}' for side in ('a', 'b') for i in range(1280)]))
    pd.DataFrame(rows).to_csv(output.with_suffix('.coverage.csv'), index=False)
    report = dict(declared=len(frame), success=sum(r['status'] == 'ok' for r in rows),
                  pooling=POOLING, exact_author_multichain_preprocessing_claimed=False)
    output.with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


def fit_bindpred(train, heldout, feature_path, output, *, iterations=2000, thread_count=2):
    """Retrain the author's tree recipe using only local training pKd labels."""
    from catboost import CatBoostRegressor

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    heldout = heldout[['sample_id', 'complex_id']].copy()
    features = LoadedFeatures.read(feature_path)
    if features.X.shape[1] != 2560:
        raise ValueError('Expected ordered ESM2-650M partner features, dimension 2560')
    x, xv = features.aligned(train.sample_id), features.aligned(heldout.sample_id)
    good, valid = np.isfinite(x).all(axis=1), np.isfinite(xv).all(axis=1)
    params = {**PARAMETERS, 'iterations': iterations, 'thread_count': thread_count}
    model = CatBoostRegressor(**params)
    started = time.monotonic()
    # Original train.py target is log10(Kd_M) = -pKd, with no sample reweighting.
    model.fit(x[good], -train.loc[good, 'pkd'].to_numpy(float))
    predicted = np.full(len(heldout), np.nan)
    if valid.any():
        predicted[valid] = -model.predict(xv[valid])
    result = heldout.copy()
    result['predicted'], result['method'], result['target'] = predicted, METHOD, 'pkd'
    finite = valid & np.isfinite(predicted)
    result['status'] = np.where(finite, 'ok', 'error')
    result['error'] = np.where(finite, '', 'Missing/nonfinite ESM2 feature or prediction')
    result['comparison_scope'] = 'matched_retraining_explicit_pooling_adaptation'
    result['training_overlap'] = 'local_split_training_only; frozen_ESM2_pretraining_not_audited'
    result.to_csv(output/'predictions.csv', index=False)
    train.loc[good, ['sample_id']].to_csv(output/'training_ids.csv', index=False)
    model.save_model(str(output/'model.cbm'))
    (output/'fit.json').write_text(json.dumps(dict(method=METHOD, parameters=params, pooling=POOLING,
        n_train_declared=len(train), n_train_success=int(good.sum()), n_test_declared=len(heldout),
        n_test_success=int(finite.sum()), heldout_targets_read=False, seconds=time.monotonic()-started,
        author_code='https://huggingface.co/hbp5181/BindPred',
        adaptations='explicit residue/chain pooling; local grouped split; CPU backend; fixed seed'), indent=2)+'\n')
    return result
