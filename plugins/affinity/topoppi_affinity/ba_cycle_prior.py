"""Frozen author BA-Cycle with seeded decoding and complete partner groups.

The score uses WT-conditioned likelihoods, mutation-last decoding, and the
bound-minus-isolated-partner cycle. It is not supervised BA-DDG, nor an exact
antisymmetric state energy. Raw dimensionless scores and RT-scaled values are
recorded separately.
"""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from topoppi_affinity.graph_features import load_graph, mutant_node_features
from topoppi_affinity.runtime import check_cancel
from topoppi_affinity.thermodynamics import R_KCAL

ALPHABET = 'ACDEFGHIKLMNPQRSTVWYX'


def mpnn_inputs(graph):
    """Author chain/gap numbering on exact graph residues, with explicit masks."""
    n = len(graph['aa_index'])
    order = np.asarray(sorted(range(n), key=lambda i: (str(graph['chain'][i]), int(graph['residue_number'][i]),
                                                      str(graph['insertion'][i]))))
    valid = np.asarray(graph['backbone_mask']).all(1) & (graph['aa_index'] < 20)
    valid &= np.isfinite(graph['backbone']).all(axis=(1, 2))
    chain_ids = list(dict.fromkeys(graph['chain'][order].astype(str)))
    residue_idx, chain_encoding = np.zeros(n, dtype=np.int64), np.zeros(n, dtype=np.int64)
    offset = 0
    for chain_number, chain in enumerate(chain_ids, 1):
        nodes = [i for i in order if str(graph['chain'][i]) == chain and valid[i]]
        sequential, previous = 0, None
        for node in nodes:
            if previous is None:
                sequential = 1
            elif np.linalg.norm(graph['backbone'][node, 1]-graph['backbone'][previous, 1]) <= 4.:
                sequential += 1
            else:
                sequential += max(2, int(graph['residue_number'][node])-int(graph['residue_number'][previous]))
            residue_idx[node], chain_encoding[node] = offset+sequential, chain_number
            previous = node
        offset += sequential + 100
    inverse = np.empty(n, dtype=np.int64)
    inverse[order] = np.arange(n)
    return {'order': order, 'inverse': inverse, 'valid': valid, 'side': graph['side'][order],
            'X': graph['backbone'][order].astype(np.float32), 'S': graph['aa_index'][order].astype(np.int64),
            'mask': valid[order].astype(np.float32), 'residue_idx': residue_idx[order],
            'chain_encoding_all': chain_encoding[order]}


def cycle_score(log_probs, states, wild_type, mutant, site_sides):
    """Return one raw score per MC draw; sum NLL only where identities differ."""
    log_probs = np.asarray(log_probs)
    sites = np.arange(len(wild_type))
    delta_nll = log_probs[:, :, sites, wild_type] - log_probs[:, :, sites, mutant]
    result = np.zeros(len(log_probs))
    for index, state in enumerate(states):
        mask = np.ones(len(sites), bool) if state == -1 else np.asarray(site_sides) == state
        result += (1. if state == -1 else -1.) * delta_nll[:, index, mask].sum(axis=1)
    return result


class BACycleEncoder:
    def __init__(self, repository, checkpoint, *, device='cpu', threads=1):
        self.repository, self.checkpoint = Path(repository).resolve(), Path(checkpoint).resolve()
        self.device = torch.device(device)
        torch.set_num_threads(threads)
        training = self.repository / 'training'
        sys.path.insert(0, str(training))
        spec = importlib.util.spec_from_file_location('ba_cycle_author_mpnn', training / 'protein_mpnn_utils.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        state = torch.load(self.checkpoint, map_location='cpu', weights_only=False)
        self.model = module.ProteinMPNN(num_letters=21, node_features=128, edge_features=128, hidden_dim=128,
            num_encoder_layers=3, num_decoder_layers=3, augment_eps=0., k_neighbors=int(state['num_edges']))
        self.model.load_state_dict(state['model_state_dict'], strict=True)
        self.model.requires_grad_(False).to(self.device).eval()

    def site_log_probs(self, graph, nodes, *, mc_samples=8, seed=42, cancel_event=None, log=None):
        inputs = mpnn_inputs(graph)
        nodes = np.asarray(nodes, dtype=np.int64)
        if not inputs['valid'][nodes].all():
            missing = graph['residue_key'][nodes[~inputs['valid'][nodes]]].astype(str).tolist()
            raise ValueError(f'Mutation sites lack complete canonical N/CA/C/O input: {missing}')
        mutation_flag = np.zeros(len(inputs['order']), dtype=np.float32)
        mutation_flag[inputs['inverse'][nodes]] = 1.
        states = [-1, *sorted(set(graph['side'][nodes].tolist()))]
        results = np.empty((mc_samples, len(states), len(nodes), 21), dtype=np.float32)
        residue_idx = torch.as_tensor(inputs['residue_idx'], device=self.device)[None]
        chain_encoding = torch.as_tensor(inputs['chain_encoding_all'], device=self.device)[None]
        device_ids = [] if self.device.type == 'cpu' else [self.device.index or torch.cuda.current_device()]
        with torch.inference_mode(), torch.random.fork_rng(devices=device_ids):
            for draw in range(mc_samples):
                check_cancel(cancel_event)
                if log:
                    log(f"BA-Cycle draw {draw + 1}/{mc_samples}")
                for state_index, side in enumerate(states):
                    check_cancel(cancel_event)
                    mask = inputs['mask'] * (1. if side == -1 else inputs['side'] == side)
                    # Match author padded isolation, preserving all chains on the selected partner side.
                    x = torch.as_tensor(inputs['X'] * mask[:, None, None], device=self.device)[None]
                    sequence = torch.as_tensor(inputs['S'] * mask.astype(np.int64), device=self.device)[None]
                    mask_tensor = torch.as_tensor(mask, device=self.device)[None]
                    chain_m = torch.as_tensor(mutation_flag * mask, device=self.device)[None]
                    # Common random numbers couple cycle states; each draw remains an author decoding sample.
                    torch.manual_seed(seed + draw)
                    log_probs = self.model.deterministic_forward(x, sequence, mask_tensor, chain_m, residue_idx, chain_encoding)
                    results[draw, state_index] = log_probs[0, inputs['inverse'][nodes]].cpu().numpy()
        return results, np.asarray(states, dtype=np.int64)

    def provenance(self, *, mc_samples, seed):
        source_record = self.repository / '.topoppi-source.json'
        if source_record.is_file():
            commit = json.loads(source_record.read_text(encoding='utf-8'))['commit']
        else:
            commit = subprocess.check_output(['git', '-C', str(self.repository), 'rev-parse', 'HEAD'], text=True).strip()
        return {'method': f'BA-Cycle-partner-MC{mc_samples}', 'target': 'ddg', 'protocol': 'native_pretrained_adapted_input',
            'source_url': 'https://github.com/aim-uofa/BA-DDG', 'paper': 'https://arxiv.org/abs/2410.09543',
            'commit': commit, 'checkpoint': str(self.checkpoint), 'training_overlap': 'possible',
            'license': 'Author README: BSD-2-Clause for non-commercial academic use; commercial use requires contacting author',
            'training_data': 'ProteinMPNN soluble structure/sequence pretraining; no binding-affinity fine-tuning',
            'training_overlap_definition': 'Structural pretraining exposure possible; no SKEMPI supervised BA-DDG weights loaded',
            'checkpoint_scalar': 1., 'mc_samples': mc_samples, 'seed': seed, 'device': str(self.device),
            'backbone_noise_A': 0., 'decoding': 'Author mutation-last chain_M=mutation_flag, seed+draw; shared randomness across cycle states',
            'cache_key': 'source/complex_id/exact mutated-site set; target amino acids do not enter WT-conditioned encoder',
            'partner_policy': 'Isolated complete partner sides retain internal chains; differs from author per-mutated-chain isolation for multichain partners',
            'coverage': 'Graph residues lacking any canonical N/CA/C/O are masked; affected mutation rows fail explicitly',
            'units': 'author_score is dimensionless NLL difference with scalar1; predicted is R_kcal*T*author_score kcal/mol',
            'sign': 'Positive cycle NLL change predicts weaker mutant binding',
            'symmetry': 'WT-context scoring is not automatically antisymmetric when recomputing reversed MT context',
            'labels_used_by_encoder': False}


def run_ba_cycle(root, output, repository, checkpoint, *, device='cpu', mc_samples=8, seed=42,
                 max_complexes=None, max_rows=None, manifest=None, graph_index=None,
                 cancel_event=None, log=None):
    log = log or (lambda text: print(text, file=sys.stderr, flush=True))
    check_cancel(cancel_event)
    root, output = Path(root), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    columns = ['sample_id', 'source', 'complex_id', 'mutations', 'n_mutations', 'temperature_K']
    if manifest is None or graph_index is None:
        raise ValueError('Supply the input manifest and graph_index explicitly')
    manifest = Path(manifest)
    frame = pd.read_csv(manifest, usecols=columns, keep_default_na=False,
                        dtype={"sample_id": str, "source": str, "complex_id": str})
    if max_complexes is not None:
        selected = sorted(frame.complex_id.unique())[:max_complexes]
        frame = frame[frame.complex_id.isin(selected)].copy()
    if max_rows is not None:
        frame = frame.iloc[:max_rows].copy()
    frame.to_csv(output / 'inputs.csv', index=False)
    graphs = pd.read_csv(graph_index, keep_default_na=False, dtype={"source": str, "complex_id": str})
    graphs = graphs.drop_duplicates(['source', 'complex_id']).set_index(['source', 'complex_id'])
    encoder = BACycleEncoder(repository, checkpoint, device=device)
    provenance = encoder.provenance(mc_samples=mc_samples, seed=seed)
    provenance.update(input_manifest=str(manifest), graph_manifest=str(graph_index))
    rows, caches = [], []
    start = time.monotonic()
    for (source, complex_id), group in frame.groupby(['source', 'complex_id'], sort=True):
        check_cancel(cancel_event)
        graph = load_graph(graphs.loc[(source, complex_id), 'graph_path'])
        prepared, site_patterns = [], {}
        for observation in group.itertuples():
            try:
                _, mutant, flag = mutant_node_features(graph, observation.mutations)
                sites = tuple(np.flatnonzero(flag).tolist())
                site_patterns.setdefault(sites, len(site_patterns))
                prepared.append((observation, mutant, sites))
            except (ValueError, KeyError) as exc:
                rows.append({'sample_id': observation.sample_id, 'source': source, 'complex_id': complex_id,
                             'target': 'ddg', 'predicted': np.nan, 'author_score': np.nan, 'status': 'error', 'error': str(exc)})
        case = output / 'log_probs' / source / complex_id.replace(':', '_').replace(',', '+')
        case.mkdir(parents=True, exist_ok=True)
        for sites, pattern in site_patterns.items():
            check_cancel(cancel_event)
            selected_rows = [(row, mutant) for row, mutant, row_sites in prepared if row_sites == sites]
            cache = case / f'sites_{pattern:04d}.npz'
            pattern_start = time.monotonic()
            try:
                nodes = np.asarray(sites)
                if cache.exists():
                    with np.load(cache, allow_pickle=False) as stored:
                        if not np.array_equal(stored['nodes'], nodes) or int(stored['seed']) != seed or int(stored['mc_samples']) != mc_samples:
                            raise ValueError('Cached site mapping/MC configuration differs; use a separate output directory')
                        log_probs, states = stored['log_probs'], stored['states']
                else:
                    log_probs, states = encoder.site_log_probs(graph, nodes, mc_samples=mc_samples, seed=seed,
                        cancel_event=cancel_event, log=log)
                    np.savez_compressed(cache, log_probs=log_probs, nodes=nodes, states=states,
                        residue_key=graph['residue_key'][nodes], seed=seed, mc_samples=mc_samples)
                for observation, mutant in selected_rows:
                    draws = cycle_score(log_probs, states, graph['aa_index'][nodes], mutant[nodes], graph['side'][nodes])
                    score, temperature = float(draws.mean()), float(observation.temperature_K)
                    if not np.isfinite(temperature) or temperature <= 0:
                        raise ValueError('RT conversion requires an explicit positive measurement temperature')
                    rows.append({'sample_id': observation.sample_id, 'source': source, 'complex_id': complex_id,
                        'target': 'ddg', 'method': provenance['method'], 'predicted': R_KCAL*temperature*score,
                        'author_score': score, 'mc_std_author_score': float(draws.std(ddof=1)) if len(draws)>1 else 0.,
                        'temperature_K': temperature, 'n_mutations': observation.n_mutations, 'status': 'ok', 'error': ''})
                caches.append({'source': source, 'complex_id': complex_id, 'pattern': pattern, 'n_sites': len(nodes),
                               'cache_path': str(cache), 'n_rows': len(selected_rows), 'seconds': time.monotonic()-pattern_start})
            except (ValueError, KeyError, RuntimeError) as exc:
                rows.extend({'sample_id': observation.sample_id, 'source': source, 'complex_id': complex_id,
                             'target': 'ddg', 'method': provenance['method'], 'predicted': np.nan, 'author_score': np.nan,
                             'status': 'error', 'error': str(exc)} for observation, _ in selected_rows)
            pd.DataFrame(rows).to_csv(output / 'ba_cycle.csv', index=False)
        log(f'BA-Cycle {len(rows)}/{len(frame)} {complex_id}: {len(site_patterns)} site patterns '
            f'{time.monotonic()-start:.1f}s')
    result = pd.DataFrame(rows).set_index('sample_id').loc[frame.sample_id].reset_index()
    result.to_csv(output / 'ba_cycle.csv', index=False)
    pd.DataFrame(caches).to_csv(output / 'cache_manifest.csv', index=False)
    provenance.update(n_attempted=len(frame), n_success=int(result.status.eq('ok').sum()), n_complexes=frame.complex_id.nunique(),
                      n_successful_site_patterns=len(caches), elapsed_seconds=time.monotonic()-start)
    (output / 'ba_cycle.json').write_text(json.dumps(provenance, indent=2)+'\n')
    return result
