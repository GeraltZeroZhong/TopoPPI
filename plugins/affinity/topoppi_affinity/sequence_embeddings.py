"""Frozen protein-language-model residue representations for structure-matched WT/MT chains.

ESM-2: Lin et al., Science 2023, facebookresearch/esm. No affinity supervision
is used. Sequences contain the observed structural residues in author order;
long chains use overlapping windows, with representations averaged in overlaps.
"""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Data.PDBData import protein_letters_3to1_extended
from topoppi.atlas.footprints import source_atom_residue_labels
from topoppi.io.io_loader import PDBLoader

from topoppi_affinity.features import _chains, parse_mutations
from topoppi_affinity.runtime import check_cancel


def chain_records(structure_path, chains_a, chains_b, *, prepared=None):
    loader = PDBLoader(str(structure_path)) if prepared is None else None
    result = {}
    for chain in sorted(set(_chains(chains_a) + _chains(chains_b))):
        if prepared is None:
            _, atoms = loader.get_chain_data(chain)
        else:
            atoms = [atom for partner in prepared.partners for atom in partner.atoms
                     if str(atom.get_parent().get_parent().id) == chain]
        labels = source_atom_residue_labels(atoms)
        records = {}
        for label, atom in zip(labels, atoms, strict=True):
            if str(label) not in records:
                residue = atom.get_parent()
                records[str(label)] = (protein_letters_3to1_extended.get(residue.resname, 'X'),
                                      int(residue.id[1]), str(residue.id[2]).strip())
        result[chain] = records
    return result


def prepare_sequences(manifests, output_dir, *, prepared=None):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.concat([pd.read_csv(p, keep_default_na=False, dtype={"sample_id": str, "source": str,
        "complex_id": str, "chains_a": str, "chains_b": str}) for p in manifests], ignore_index=True)
    frame = frame.drop_duplicates('sample_id')
    registry = output_dir / 'sequences.csv'
    existing = pd.read_csv(registry, keep_default_na=False) if registry.exists() else pd.DataFrame()
    sequences = dict(zip(existing.sequence, existing.sequence_id, strict=True)) if len(existing) else {}
    wt_path, mt_path = output_dir / 'wild_type_chains.csv', output_dir / 'mutant_chains.csv'
    wt_records = pd.read_csv(wt_path, keep_default_na=False, dtype={"chain": str, "source": str, "complex_id": str}).to_dict('records') if wt_path.exists() else []
    mutant_records = pd.read_csv(mt_path, keep_default_na=False, dtype={"chain": str, "sample_id": str, "complex_id": str}).to_dict('records') if mt_path.exists() else []
    # A requested variant replaces its whole chain mapping, including chains
    # that have returned to WT since this sample ID was last prepared.
    current_ids = set(frame.sample_id)
    mutant_records = [record for record in mutant_records if record['sample_id'] not in current_ids]

    def register(sequence):
        if sequence not in sequences:
            sequences[sequence] = f'sequence_{len(sequences):06d}'
        return sequences[sequence]

    for index, ((source, complex_id), subset) in enumerate(frame.groupby(['source', 'complex_id'], sort=True)):
        row = subset.iloc[0]
        chains = chain_records(row.structure_path, row.chains_a, row.chains_b, prepared=prepared)
        for chain, records in chains.items():
            sequence = ''.join(v[0] for v in records.values())
            wt_records.append({'source': source, 'complex_id': complex_id, 'chain': chain, 'sequence_id': register(sequence),
                               'sequence': sequence, 'residue_keys': json.dumps(list(records))})
        for row in subset.itertuples():
            if not row.mutations:
                continue
            mutations = parse_mutations(row.mutations)
            for chain in sorted({m['chain'] for m in mutations}):
                records = list(chains[chain].values())
                lookup = {(v[1], v[2]): i for i, v in enumerate(records)}
                sequence = [v[0] for v in records]
                for mutation in mutations:
                    if mutation['chain'] != chain:
                        continue
                    position = lookup[(int(mutation['position']), mutation['insertion'])]
                    if sequence[position] != mutation['wt']:
                        raise ValueError(f'{row.sample_id}: structural sequence WT mismatch')
                    sequence[position] = mutation['mt']
                mutant_records.append({'sample_id': row.sample_id, 'complex_id': complex_id,
                                       'chain': chain, 'sequence_id': register(''.join(sequence))})
        if (index + 1) % 25 == 0:
            print(f'Sequence mapping {index + 1} complexes; {len(sequences)} unique sequences', file=sys.stderr, flush=True)
    pd.DataFrame(wt_records).drop_duplicates(['source', 'complex_id', 'chain'], keep='last').to_csv(wt_path, index=False)
    pd.DataFrame(mutant_records, columns=['sample_id', 'complex_id', 'chain', 'sequence_id']).drop_duplicates(['sample_id', 'chain'], keep='last').to_csv(mt_path, index=False)
    pd.DataFrame([{'sequence_id': identifier, 'sequence': sequence, 'length': len(sequence)}
                  for sequence, identifier in sequences.items()]).to_csv(output_dir / 'sequences.csv', index=False)
    return len(sequences)


def window_starts(length, window=1000, overlap=128):
    if length <= window:
        return [0]
    starts = list(range(0, length - window + 1, window - overlap))
    if starts[-1] != length - window:
        starts.append(length - window)
    return starts


def extract_embeddings(sequence_dir, checkpoint=None, *, device='cuda', limit=None, model_name='esmc_600m', sequence_ids=None,
                       shared_cache=None, cancel_event=None, progress=None, log=None, cache_dir=None):
    import torch

    log = log or (lambda text: print(text, file=sys.stderr, flush=True))
    check_cancel(cancel_event)
    torch.set_num_threads(2)
    sequence_dir = Path(sequence_dir)
    table = pd.read_csv(sequence_dir / 'sequences.csv', keep_default_na=False)
    if sequence_ids is not None:
        table = table[table.sequence_id.isin(set(sequence_ids))]
    if limit is not None:
        table = table.iloc[:limit]
    destination = sequence_dir / 'embeddings'
    destination.mkdir(exist_ok=True)
    shared_cache = Path(shared_cache) if shared_cache else None
    if shared_cache:
        shared_cache.mkdir(parents=True, exist_ok=True)

    def shared_path(sequence):
        # A sequence digest replaces a complete ESMC forward pass on repeated input.
        return shared_cache / (hashlib.sha256(sequence.encode()).hexdigest() + '.npz') if shared_cache else None

    def matches(path, sequence):
        if path is None or not path.exists():
            return False
        try:
            with np.load(path, allow_pickle=False) as cached:
                return str(cached['sequence']) == sequence and str(cached['model_name']) == model_name
        except (OSError, ValueError, KeyError):
            return False

    missing = []
    for row in table.itertuples():
        check_cancel(cancel_event)
        path = destination / f'{row.sequence_id}.npz'
        shared = shared_path(row.sequence)
        if matches(path, row.sequence):
            if shared and not matches(shared, row.sequence):
                shutil.copyfile(path, shared)
            continue
        if matches(shared, row.sequence):
            shutil.copyfile(shared, path)
            continue
        missing.append(row.sequence_id)
    requested = len(table)
    if not missing:
        log(f'Reused cached embeddings for all {requested} requested sequences')
        if progress:
            progress({'stage': 'embeddings', 'status': 'ok', 'current': requested, 'total': requested,
                      'message': f'Embeddings ready: {requested}/{requested} (cached)'})
        return {'model': model_name, 'requested': len(table), 'newly_embedded': 0,
                'protocol': 'frozen sequence pretraining; verified cached sequences and model names'}
    table = table[table.sequence_id.isin(missing)]
    if model_name == 'esmc_600m':
        from esm.models.esmc import ESMC
        from esm.tokenization import get_esmc_model_tokenizers

        from topoppi_affinity.resources import download_models, esmc_checkpoint

        checkpoint = Path(checkpoint) if checkpoint else esmc_checkpoint(cache_dir=cache_dir)
        if checkpoint is None:
            log('Downloading ESMC-600M from its official source (about 2.3 GB); subsequent runs reuse it')
            download_models(tasks=(), include_esmc=True, progress=progress, cancel_event=cancel_event, cache_dir=cache_dir)
            checkpoint = esmc_checkpoint(cache_dir=cache_dir)
        check_cancel(cancel_event)
        # Same architecture, state, device dtype and attention default as esm 3.2.1's frozen builder.
        with torch.device(device):
            model = ESMC(d_model=1152, n_heads=18, n_layers=36,
                         tokenizer=get_esmc_model_tokenizers(), use_flash_attn=True).eval()
        model.load_state_dict(torch.load(checkpoint, map_location=device, weights_only=True))
        if torch.device(device).type != 'cpu':
            model = model.to(torch.bfloat16)
        dimension = 1152
        model_source = 'https://huggingface.co/biohub/esmc-600m-2024-12'
    else:
        from transformers import AutoModel, AutoTokenizer
        tokenizer = AutoTokenizer.from_pretrained(str(checkpoint), local_files_only=True)
        model = AutoModel.from_pretrained(str(checkpoint), local_files_only=True,
            torch_dtype=torch.float16 if device.startswith('cuda') else torch.float32,
            low_cpu_mem_usage=True).to(device).eval()
        dimension = model.config.hidden_size
        model_source = 'https://github.com/facebookresearch/esm'
    started, completed = time.monotonic(), 0
    for row in table.itertuples():
        check_cancel(cancel_event)
        path = destination / f'{row.sequence_id}.npz'
        if path.exists():
            with np.load(path, allow_pickle=False) as cached:
                if str(cached['sequence']) == row.sequence and str(cached['model_name']) == model_name:
                    continue
        accumulated = np.zeros((len(row.sequence), dimension), dtype=np.float32)
        count = np.zeros((len(row.sequence), 1), dtype=np.float32)
        for start in window_starts(len(row.sequence)):
            check_cancel(cancel_event)
            sequence = row.sequence[start:start + 1000]
            with torch.inference_mode():
                if model_name == 'esmc_600m':
                    tokens = model._tokenize([sequence])
                    value = model(sequence_tokens=tokens).embeddings[0, 1:len(sequence) + 1].float().cpu().numpy()
                else:
                    tokens = tokenizer(sequence, return_tensors='pt', add_special_tokens=True).to(device)
                    value = model(**tokens).last_hidden_state[0, 1:len(sequence) + 1].float().cpu().numpy()
            accumulated[start:start + len(sequence)] += value
            count[start:start + len(sequence)] += 1
        np.savez_compressed(path, embedding=(accumulated / count).astype(np.float16), sequence=row.sequence, model_name=model_name)
        if shared_cache:
            shutil.copyfile(path, shared_path(row.sequence))
        completed += 1
        message = f'Embedded {completed}/{len(table)} new sequences; {time.monotonic() - started:.1f}s'
        log(message)
        if progress:
            progress({'stage': 'embeddings', 'status': 'running', 'current': completed,
                      'total': len(table), 'message': message})
    provenance = {'model': model_name,
        'source_url': model_source, 'checkpoint': str(checkpoint),
        'protocol': 'frozen sequence pretraining; no affinity fine-tuning',
        'structural_sequences': 'Observed residues, author ordering, no unresolved residues imputed',
        'long_chain_window': 1000, 'window_overlap': 128, 'overlap_aggregation': 'mean',
        'newly_embedded': completed, 'requested': requested, 'device': device,
        'elapsed_seconds': time.monotonic() - started}
    (sequence_dir / 'embedding_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    return provenance
