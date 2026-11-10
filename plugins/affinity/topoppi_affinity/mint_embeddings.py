"""Load and apply the MINT sequence encoder."""
from __future__ import annotations

import gc
import sys
from pathlib import Path

import torch


def load_mint(repository, checkpoint, device):
    sys.path.insert(0, str(Path(repository).resolve()))
    from mint.helpers.extract import CollateFn, load_config
    from mint.model.esm2 import ESM2

    config = load_config(Path(repository) / 'data/esm2_t33_650M_UR50D.json')
    payload = torch.load(checkpoint, map_location='cpu', weights_only=False, mmap=True)
    weights = {key.removeprefix('model.'): value for key, value in payload['state_dict'].items()}
    with torch.device('meta'):
        model = ESM2(num_layers=config.encoder_layers, embed_dim=config.encoder_embed_dim,
            attention_heads=config.encoder_attention_heads, token_dropout=config.token_dropout, use_multimer=True)
    model.load_state_dict(weights, strict=True, assign=True)
    dtype = torch.bfloat16 if str(device).startswith('cuda') else torch.float32
    model = model.to(device=device, dtype=dtype).eval().requires_grad_(False)
    del weights, payload
    gc.collect()
    return model, CollateFn(truncation_seq_length=None)


def encode_multimer(model, converter, sequences):
    tokens, chains = converter([sequences])
    device = next(model.parameters()).device
    with torch.inference_mode():
        encoded = model(tokens.to(device), chains.to(device), repr_layers=[model.num_layers])['representations'][model.num_layers][0]
    embeddings, start = [], 0
    for sequence in sequences:
        embeddings.append(encoded[start+1:start+len(sequence)+1].float().cpu().numpy())
        start += len(sequence) + 2
    return embeddings
