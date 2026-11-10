"""Frozen author PPLM representations, mapped to native TopoPPI residues.

PPLM (Liu et al., Nature Communications 2026) learns paired-sequence context.
Only the unsupervised language-model checkpoint is used here. For multi-chain
partners, each cross-partner chain pair is encoded separately; representations
of a chain are averaged over the opposite partner's chains. No synthetic linker
or affinity-finetuned checkpoint is introduced.
"""
from __future__ import annotations

import gc
import json
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import torch


def load_pplm(repository, checkpoint, device='cuda'):
    """mmap avoids materializing the checkpoint's unused optimizer tensors."""
    sys.path.insert(0, str(Path(repository).resolve()))
    from pplm import PPLM, Alphabet

    payload = torch.load(checkpoint, map_location='cpu', weights_only=False, mmap=True)
    params = payload['param']
    alphabet = Alphabet.from_architecture()
    with torch.device('meta'):
        model = PPLM(num_layers=params['encoder_layers'], embed_dim=params['encoder_embed_dim'],
                     attention_heads=params['encoder_attention_heads'], token_dropout=False, alphabet=alphabet)
    model.load_state_dict(payload['model'], strict=True, assign=True)
    dtype = torch.bfloat16 if str(device).startswith('cuda') else torch.float32
    model = model.to(device=device, dtype=dtype).eval().requires_grad_(False)
    del payload
    gc.collect()
    return model, alphabet.get_batch_converter()


def encode_pair(model, converter, sequence_a, sequence_b, *, return_attention=False):
    device = next(model.parameters()).device
    _, _, a = converter([('a', sequence_a)])
    _, _, b = converter([('b', sequence_b)])
    tokens = torch.cat([a, b], dim=-1).to(device)
    split, length = a.shape[1], tokens.shape[1]
    mask = torch.ones((length, length), device=device, dtype=next(model.parameters()).dtype)
    mask[:split, :split] = 0
    mask[split:, split:] = 0
    with torch.inference_mode():
        result = model(tokens, mask, repr_layers=[model.num_layers], need_head_weights=return_attention)
        embedding = result['representations'][model.num_layers][0]
    return (embedding[1:len(sequence_a) + 1].float().cpu().numpy(),
            embedding[split + 1:split + len(sequence_b) + 1].float().cpu().numpy())


@lru_cache(maxsize=32)
def _pair_arrays(path):
    with np.load(path, allow_pickle=False) as source:
        return source['embedding_a'].astype(np.float32), source['embedding_b'].astype(np.float32)


class PairedResidueFeatures:
    def __init__(self, root, output):
        self.output = Path(output)
        assignments = pd.read_csv(self.output / 'assignments.csv', keep_default_na=False)
        self.assignments = dict(tuple(assignments.groupby(['source', 'complex_id', 'state'], sort=False)))
        chains = pd.read_csv(Path(root) / 'sequences/wild_type_chains.csv', keep_default_na=False)
        self.positions = {(r.source, r.complex_id, r.chain): {key: i for i, key in enumerate(json.loads(r.residue_keys))}
                          for r in chains.itertuples()}

    def features(self, source, complex_id, state, graph):
        rows = self.assignments[(source, complex_id, state)]
        values = np.zeros((len(graph['x']), 1280), dtype=np.float32)
        count = np.zeros(len(values), dtype=np.float32)
        for row in rows.itertuples():
            a, b = _pair_arrays(str(self.output / 'embeddings' / (row.pair_id + '.npz')))
            for chain, embedding in ((row.chain_a, a if row.a_is_first else b), (row.chain_b, b if row.a_is_first else a)):
                selected = np.flatnonzero(graph['chain'] == chain)
                lookup = self.positions[(source, complex_id, chain)]
                positions = [lookup[str(graph['residue_key'][index])] for index in selected]
                values[selected] += embedding[positions]
                count[selected] += 1
        if (count == 0).any():
            raise ValueError(f'{source}/{complex_id}: PPLM chain pair does not cover every graph residue')
        return values / count[:, None]
