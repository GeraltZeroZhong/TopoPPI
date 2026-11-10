"""Create residue-aligned embeddings using an explicitly supplied local checkpoint."""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd


def _write_chain_registry(interfaces, output_dir):
    from Bio.Data.PDBData import protein_letters_3to1_extended

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    sequences, chains, known, complex_ids = [], [], {}, set()
    for complex_id, interface in interfaces:
        complex_ids.add(complex_id)
        by_chain = {}
        for partner in interface.partners:
            for atom, key in zip(partner.atoms, partner.residue_keys, strict=True):
                residue = atom.parent
                if atom.name != "CA" or residue.resname not in protein_letters_3to1_extended:
                    continue
                chain_id = str(residue.parent.id)
                by_chain.setdefault(chain_id, []).append((str(key), residue.resname))
        for chain_id, residues in sorted(by_chain.items()):
            sequence = ''.join(protein_letters_3to1_extended[name] for _, name in residues)
            keys = [key for key, _ in residues]
            if sequence not in known:
                known[sequence] = f'hotspot_wt_{len(known):04d}'
                sequences.append({'sequence_id': known[sequence], 'sequence': sequence})
            chains.append({'complex_id': complex_id, 'chain': chain_id,
                           'sequence_id': known[sequence], 'residue_keys': json.dumps(keys)})
    pd.DataFrame(sequences, columns=['sequence_id', 'sequence']).to_csv(output / 'sequences.csv', index=False)
    pd.DataFrame(chains, columns=['complex_id', 'chain', 'sequence_id', 'residue_keys']).to_csv(
        output / 'chain_sequences.csv', index=False)
    return {'n_complexes': len(complex_ids), 'n_chains': len(chains),
            'n_unique_sequences': len(sequences), 'labels_used': False}


def prepare_interface_embedding_registry(prepared, output_dir, *, complex_id=None):
    """Reuse Core-selected atoms for full observed WT chains; no file parsing."""
    identity = complex_id or prepared.source.get("complex_id") or prepared.analysis_id
    return _write_chain_registry([(identity, prepared)], output_dir)


def prepare_chain_embedding_registry(complexes, output_dir):
    """Create a label-free registry for complete observed PDB/mmCIF chains.

    The legacy table interface needs complex_id, structure_path, chains_a and
    chains_b. Core controls first-model, heavy-atom and alternate-conformer
    selection. Only residues with a resolved CA enter the sequence registry.
    """
    from topoppi import get_interface, load_structure

    def interfaces():
        columns = ['complex_id', 'structure_path', 'chains_a', 'chains_b']
        for row in complexes[columns].drop_duplicates().itertuples(index=False):
            yield row.complex_id, get_interface(load_structure(row.structure_path), row.chains_a, row.chains_b)

    return _write_chain_registry(interfaces(), output_dir)


def embed_sequences(scan_dir, checkpoint, *, device='cuda'):
    import importlib.metadata

    import torch
    from esm.models.esmc import ESMC
    from esm.tokenization import get_esmc_model_tokenizers

    root = Path(scan_dir)
    table = pd.read_csv(root / 'sequences.csv', keep_default_na=False)
    destination = root / 'sequence_embeddings'
    destination.mkdir(exist_ok=True)
    # Whole observed chains, no interface crops or missing-residue imputation.
    if table.sequence.str.len().max() > 1000:
        raise ValueError('Whole-chain embedding supports observed chains of at most 1000 residues')
    torch.set_num_threads(2)
    started = time.monotonic()
    model = ESMC(d_model=1152, n_heads=18, n_layers=36,
                 tokenizer=get_esmc_model_tokenizers(), use_flash_attn=device.startswith('cuda')).eval()
    model.load_state_dict(torch.load(checkpoint, map_location='cpu', weights_only=True))
    model = model.to(device)
    if device.startswith('cuda'):
        model = model.to(torch.bfloat16)
    for index, row in enumerate(table.itertuples(index=False)):
        with torch.inference_mode():
            tokens = model._tokenize([row.sequence])
            value = model(sequence_tokens=tokens).embeddings[0, 1:len(row.sequence) + 1].float().cpu().numpy()
        if value.shape != (len(row.sequence), 1152) or not np.isfinite(value).all():
            raise ValueError(f'Invalid embeddings for {row.sequence_id}')
        np.savez_compressed(destination / f'{row.sequence_id}.npz', sequence=row.sequence,
                            embedding=value, model_name='esmc_600m')
        print(f'Hotspot WT sequence {index + 1}/{len(table)}: {len(row.sequence)} residues', flush=True)
    protocol = {'encoder': 'frozen ESMC 600M architecture', 'checkpoint': str(Path(checkpoint).resolve()),
                'source_url': 'https://huggingface.co/EvolutionaryScale/esmc-600m-2024-12',
                'esm_version': importlib.metadata.version('esm'), 'torch_version': torch.__version__,
                'device': device, 'dtype': 'bfloat16' if device.startswith('cuda') else 'float32',
                'n_sequences': len(table), 'embedding_dimension': 1152,
                'new_forward_passes': len(table), 'affinity_predictions_read': False,
                'readout': 'focal WT residue final layer only; no identity/contact/geometry/partner features',
                'elapsed_seconds': time.monotonic() - started}
    if device.startswith('cuda'):
        protocol['peak_allocated_MiB'] = torch.cuda.max_memory_allocated() / 2**20
    (root / 'sequence_embedding_protocol.json').write_text(json.dumps(protocol, indent=2) + '\n')
    return protocol
