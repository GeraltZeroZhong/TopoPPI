"""Pure released RDE encodings on deterministic, mutation-independent WT patches.

This is frozen structural representation transfer, not the supervised RDE-Net
binding predictor. Only the original RDE.pt encoder and author input parser run.
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from Bio.PDB import PDBParser

from topoppi_affinity.features import AA as TOPO_AA


class FrozenRDEEncoder:
    def __init__(self, repository, *, device='cpu', patch_size=128, core_size=32, threads=2):
        self.repository = Path(repository).resolve()
        self.device, self.patch_size, self.core_size = device, patch_size, core_size
        if core_size > patch_size or core_size < 1:
            raise ValueError('Require 1 <= core_size <= patch_size')
        sys.path.insert(0, str(self.repository))
        from rde.models.rde import CircularSplineRotamerDensityEstimator
        from rde.utils.data import PaddingCollate
        from rde.utils.protein.parsers import parse_biopython_structure
        from rde.utils.transforms._base import _get_CB_positions, _index_select_data
        from rde.utils.transforms.select_atom import SelectAtom

        checkpoints = list((self.repository / 'trained_models').rglob('RDE.pt'))
        if len(checkpoints) != 1:
            raise ValueError('Expected exactly one pure unsupervised RDE.pt checkpoint')
        self.checkpoint_path = checkpoints[0]
        checkpoint = torch.load(self.checkpoint_path, map_location='cpu', weights_only=False)
        torch.set_num_threads(threads)
        self.model = CircularSplineRotamerDensityEstimator(checkpoint['config'].model)
        self.model.load_state_dict(checkpoint['model'], strict=True)
        self.model.requires_grad_(False).to(device).eval()
        self.parse = parse_biopython_structure
        self.select_atoms, self.collate = SelectAtom('backbone+CB'), PaddingCollate()
        self.cb_positions, self.select_indices = _get_CB_positions, _index_select_data
        self.iteration = int(checkpoint['iteration'])
        self.dimension = int(checkpoint['config'].model.encoder.node_feat_dim)

    def provenance(self):
        commit = subprocess.check_output(['git', '-C', str(self.repository), 'rev-parse', 'HEAD'], text=True).strip()
        return {
            'encoder': 'CircularSplineRotamerDensityEstimator.encode', 'checkpoint': str(self.checkpoint_path),
            'source_url': 'https://github.com/luost26/RDE-PPI', 'license': 'Apache-2.0', 'commit': commit,
            'checkpoint_iteration': self.iteration, 'dimension': self.dimension, 'device': str(self.device),
            'pretraining': 'PDB-REDO rotamer maximum likelihood; no supervised binding labels',
            'supervised_ddg_weights_loaded': False, 'rdelinear_calibration_loaded': False,
            'context': 'complete declared partner groups; native WT chi; chi_masked_flag=False',
            'patch_mode': 'deterministic 32-neighbor core cover, author distance-to-core 128-neighbor patch; overlap mean',
            'patch_size': self.patch_size, 'core_size': self.core_size,
            'mutation_independent': True, 'labels_used': False,
            'interpretation': 'Algorithmic representation transfer; not original supervised RDE-Net or a mutant state energy',
            'training_overlap': 'Affinity labels excluded; structural pretraining overlap is possible and not audited',
        }

    def parse_structure(self, row):
        original = PDBParser(QUIET=True).get_structure('rde', row.structure_path)[0]
        chains = sorted(set(str(row.chains_a).split(',') + str(row.chains_b).split(',')))
        selected = [original[chain] for chain in chains]
        data, sequence_map = self.parse(selected)
        if data is None:
            raise ValueError('Author parser returned no supported residues')
        data = self.select_atoms(data)
        data['chi_corrupt'] = data['chi'].clone()
        data['chi_masked_flag'] = torch.zeros_like(data['aa'], dtype=torch.bool)
        return data, sequence_map

    def encode_structure(self, row, graph):
        start = time.monotonic()
        data, sequence_map = self.parse_structure(row)
        positions = self.cb_positions(data['pos_atoms'], data['mask_atoms']).numpy()
        patches = covering_patches(positions, patch_size=self.patch_size, core_size=self.core_size)
        total = np.zeros((len(positions), self.dimension), dtype=np.float64)
        counts, core_counts = np.zeros(len(positions), dtype=np.int32), np.zeros(len(positions), dtype=np.int32)
        with torch.inference_mode():
            for core, indices in patches:
                patch = self.select_indices(data, torch.as_tensor(indices))
                batch = self.collate([patch])
                batch = {key: value.to(self.device) if isinstance(value, torch.Tensor) else value
                         for key, value in batch.items()}
                embedding = self.model.encode(batch)[0, :len(indices)].cpu().numpy()
                if not np.isfinite(embedding).all():
                    raise ValueError('Pure RDE encoder returned nonfinite values')
                total[indices] += embedding
                counts[indices] += 1
                core_counts[core] += 1
        encoded = (total / counts[:, None]).astype(np.float32)
        author_index, errors = align_author_residues(graph, sequence_map, data['aa'].numpy())
        mask = author_index >= 0
        values = np.full((len(author_index), self.dimension), np.nan, dtype=np.float32)
        values[mask] = encoded[author_index[mask]]
        graph_counts, graph_core_counts = np.zeros(len(mask), dtype=np.int32), np.zeros(len(mask), dtype=np.int32)
        graph_counts[mask], graph_core_counts[mask] = counts[author_index[mask]], core_counts[author_index[mask]]
        result = {
            'embedding': values, 'encoded_mask': mask, 'author_index': author_index,
            'patch_count': graph_counts, 'core_count': graph_core_counts, 'residue_key': graph['residue_key'],
        }
        report = {'n_nodes': len(mask), 'n_encoded': int(mask.sum()), 'n_missing': int((~mask).sum()),
                  'n_author_residues': len(positions), 'n_patches': len(patches),
                  'seconds': time.monotonic() - start, 'mapping_errors': errors,
                  'status': 'ok' if mask.all() else 'partial'}
        return result, report


def covering_patches(positions, *, patch_size=128, core_size=32):
    """Cover each WT residue as a core member using geometry only, with no labels."""
    positions = np.asarray(positions, dtype=np.float64)
    if positions.ndim != 2 or positions.shape[1] != 3 or not np.isfinite(positions).all():
        raise ValueError('Patch positions must be finite N by 3 coordinates')
    if not len(positions) or not 1 <= core_size <= patch_size:
        raise ValueError('Nonempty positions and 1 <= core_size <= patch_size are required')
    covered = np.zeros(len(positions), dtype=bool)
    patches = []
    while not covered.all():
        seed = np.flatnonzero(~covered)[0]
        distances = np.linalg.norm(positions - positions[seed], axis=1)
        core = np.argsort(distances, kind='stable')[:core_size]
        # Guarantee progress even for experimentally coincident CA/CB coordinates.
        if seed not in core:
            core[-1] = seed
        distances = np.linalg.norm(positions[:, None] - positions[core][None], axis=-1).min(axis=1)
        remaining = np.flatnonzero(~np.isin(np.arange(len(positions)), core))
        neighbors = remaining[np.argsort(distances[remaining], kind='stable')[:max(0, patch_size-len(core))]]
        indices = np.sort(np.concatenate([core, neighbors]))
        covered[core] = True
        patches.append((core, indices))
    return patches


def align_author_residues(graph, sequence_map, author_aa):
    """Match PDB chain, number, insertion code and identity to exact graph order."""
    indices, errors = np.full(len(graph['aa_index']), -1, dtype=np.int64), []
    for node, (chain, number, insertion, aa) in enumerate(zip(
            graph['chain'], graph['residue_number'], graph['insertion'], graph['aa_index'], strict=True)):
        key = (str(chain), int(number), str(insertion).strip() or ' ')
        if key not in sequence_map:
            errors.append({'node': node, 'residue_key': str(graph['residue_key'][node]), 'reason': 'omitted_by_author_parser'})
            continue
        author_index = sequence_map[key]
        # Both schemas use alphabetical canonical AA order, but compare letters explicitly.
        topo_letter = (TOPO_AA + 'X')[int(aa)]
        author_letter = 'ACDEFGHIKLMNPQRSTVWYX'[int(author_aa[author_index])]
        if topo_letter != author_letter:
            errors.append({'node': node, 'residue_key': str(graph['residue_key'][node]),
                           'reason': f'identity_mismatch:{topo_letter}:{author_letter}'})
            continue
        indices[node] = author_index
    return indices, errors


def extract_structure_embeddings(graph_manifests, output, repository, *, device='cpu',
                                 limit_complexes=None, patch_size=128, core_size=32, threads=2):
    """Build source/complex caches directly from graph manifests, without targets."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    frame = pd.concat([pd.read_csv(path, keep_default_na=False) for path in graph_manifests], ignore_index=True)
    frame = frame[frame.status == 'ok'].drop_duplicates(['source', 'complex_id'])
    if limit_complexes is not None:
        frame = frame.head(limit_complexes)
    encoder = FrozenRDEEncoder(repository, device=device, patch_size=patch_size, core_size=core_size, threads=threads)
    metadata = encoder.provenance()
    metadata['graph_manifests'] = list(map(str, graph_manifests))
    (output / 'provenance.json').write_text(json.dumps(metadata, indent=2) + '\n')
    reports = []
    start = time.monotonic()
    for number, row in enumerate(frame.itertuples(), start=1):
        name = row.complex_id.replace(':', '_').replace(',', '+')
        path = output / row.source / f'{name}.npz'
        path.parent.mkdir(exist_ok=True)
        record = {'source': row.source, 'complex_id': row.complex_id, 'embedding_path': str(path)}
        if path.exists() and path.with_suffix('.json').exists():
            report = json.loads(path.with_suffix('.json').read_text())
            record.update(report, cached=True)
        else:
            try:
                with np.load(row.graph_path, allow_pickle=False) as archive:
                    graph = {key: archive[key] for key in ('aa_index', 'chain', 'residue_number', 'insertion', 'residue_key')}
                embedding, report = encoder.encode_structure(row, graph)
                np.savez_compressed(path, **embedding)
                path.with_suffix('.json').write_text(json.dumps(report, indent=2) + '\n')
                record.update(report, cached=False)
            except (ValueError, KeyError, IndexError) as exc:
                record.update(status='error', error=str(exc), cached=False)
        record.pop('mapping_errors', None)
        reports.append(record)
        pd.DataFrame(reports).to_csv(output / 'manifest.csv', index=False)
        print(f'Frozen RDE {number}/{len(frame)} {row.complex_id} {record["status"]} '
              f'{time.monotonic()-start:.1f}s', flush=True)
    return pd.DataFrame(reports)
