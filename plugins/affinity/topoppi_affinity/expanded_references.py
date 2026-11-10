"""Original PPIformer, Pythia-PPI and SAAMBE-3D reference inference.

Run each method in its recorded Python environment. All manifest rows are
retained, including unsupported multi-mutants and preprocessing failures.
Author supervised checkpoints are references, not same-split retraining.
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from topoppi_affinity.data import parse_mutations
from topoppi_affinity.external_absolute import write_partner_structure
from topoppi_affinity.external_inputs import prepare_numbered_complex, remap_mutations

METHODS = {
    'ppiformer': {'name': 'PPIformer-author', 'single_only': False, 'license': 'MIT',
                 'url': 'https://github.com/anton-bushuiev/PPIformer'},
    'pythia': {'name': 'Pythia-PPI-author', 'single_only': True, 'license': 'Apache-2.0',
              'url': 'https://github.com/Wublab/pythia_ppi'},
    'saambe': {'name': 'SAAMBE-3D-author', 'single_only': True, 'license': 'GPL-3.0',
              'url': 'https://github.com/delphi001/SAAMBE-3D'},
}


def declared_rows(frame, method):
    spec = METHODS[method]
    result = frame[['sample_id', 'complex_id']].copy()
    result['method'], result['target'], result['predicted'] = spec['name'], 'ddg', np.nan
    result['status'], result['error'] = 'error', ''
    result['comparison_scope'] = 'author_supervised_reference'
    result['training_overlap'] = 'possible_SKEMPI_supervised_exposure; not local-fold independent'
    unsupported = frame.mutations.map(lambda value: len(parse_mutations(value)) != 1) if spec['single_only'] else pd.Series(False, index=frame.index)
    result.loc[unsupported, ['status', 'error']] = ['not_applicable', 'Original method supports single substitutions only']
    return result


def _partner_path(record, work):
    chains = record['chains_a'].split(',') + record['chains_b'].split(',')
    return work / (record['pdb_id'].lower()+'_'+'_'.join(sorted(chains))+'.pdb')


def radius_residue_adjacency(atom_df, radius):
    """Same residue contact rule as PPIRef, without its atom-square dense matrix."""
    from scipy.spatial import cKDTree

    codes, labels = pd.factorize(atom_df.res_id, sort=True)
    labels.name = atom_df.res_id.name
    adjacency = np.zeros((len(labels), len(labels)), dtype=bool)
    if radius > 0:
        coordinates = atom_df[['x_coord', 'y_coord', 'z_coord']].to_numpy(dtype=np.float64)
        pairs = cKDTree(coordinates).query_pairs(radius, output_type='ndarray')
        if len(pairs):
            distances = np.linalg.norm(coordinates[pairs[:, 0]]-coordinates[pairs[:, 1]], axis=1)
            pairs = pairs[distances.astype(np.float32) < radius]
            i, j = codes[pairs[:, 0]], codes[pairs[:, 1]]
            adjacency[i, j] = adjacency[j, i] = True
        np.fill_diagonal(adjacency, True)
    return pd.DataFrame(adjacency, index=labels, columns=labels)


def prepare_ppiformer_interface(record, work):
    from ppiref.extraction import PPIExtractor

    source = work/(record['pdb_id'].lower()+'.pdb')
    write_partner_structure(record, source)
    partners = sorted(record['chains_a'].split(',') + record['chains_b'].split(','))
    extractor = PPIExtractor(work/'interface', radius=10., join=True, nest_out_dir=False,
                             max_workers=1, noppi_files=False)
    extractor._get_radius_adjacency = radius_residue_adjacency
    extractor.extract(source, partners=partners)
    path = work/'interface'/(record['pdb_id'].lower()+'_'+'_'.join(partners)+'.pdb')
    if not path.exists():
        raise ValueError('PPIRef did not produce a 10A interface')
    return path


class SaambePredictor:
    def __init__(self, repository, device):
        sys.path.insert(0, str(repository))
        spec = importlib.util.spec_from_file_location('saambe_author', repository/'saambe-3d.py')
        self.author = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.author)

    def prepare(self, record, work):
        path = _partner_path(record, work)
        mapping = prepare_numbered_complex(record, path)
        # The released regression features include resolution, crystallization
        # temperature and pH. PDBIO alone would silently discard these inputs.
        remarks = ''.join(line for line in Path(record['structure_path']).read_text().splitlines(True)
                          if line.startswith('REMARK'))
        path.write_text(remarks + path.read_text())
        structure = self.author.pdy.parsePDB(str(path), header=False)
        return path, mapping, structure

    def predict(self, payload, mutations):
        path, mapping, structure = payload
        mapped = remap_mutations(mutations, mapping)
        wt, chain, position, _, mt = parse_mutations(mapped)[0]
        _, errors = self.author.validate_mt_info(structure, chain, str(position), wt, mt)
        if errors:
            raise ValueError('; '.join(errors))
        # Preserve the original CLI's rounded regression output and sign.
        text = self.author.file_loop(str(path), structure, chain, str(position), wt, mt, 1, False)
        return float(text.split()[0])


class PythiaPredictor:
    def __init__(self, repository, device):
        import torch

        if device != 'cuda':
            raise ValueError('Released Pythia get_neighbor hardcodes CUDA; use --device cuda')
        sys.path.insert(0, str(repository))
        from utils import dataset
        from utils.model import Pythia_PPI, get_torch_model

        self.torch, self.dataset, self.device = torch, dataset, device
        torch.set_num_threads(2)
        encoder = get_torch_model(str(repository/'utils/pythia/pythia-p.pt'), device)
        self.model = Pythia_PPI(encoder)
        self.model.load_state_dict(torch.load(repository/'train_model/ppi_affinity.pth', map_location=device,
                                              weights_only=True))
        self.model.to(device).eval()

    def prepare(self, record, work):
        from Bio.PDB import PDBParser

        path = _partner_path(record, work)
        write_partner_structure(record, path)
        protbb, _ = self.dataset.read_pdb_to_protbb(str(path), return_chain_dict=True)
        seq_index = protbb.seq.detach().cpu().numpy().copy()
        # Match the author's complete-backbone residue iteration, including
        # insertion codes, without changing the residue-number edge features.
        structure = PDBParser(QUIET=True).get_structure('ppi', path)[0]
        keys = [(chain.id, residue.id[1], residue.id[2].strip())
                for chain in structure for residue in chain
                if residue.id[0] == ' ' and all(atom in residue for atom in ('N', 'C', 'O', 'CA'))]
        if len(keys) != len(seq_index):
            raise ValueError('Pythia residue mapping differs from native parser')
        node, edge, _ = self.dataset.get_neighbor(protbb, noise_level=0, mask=False)
        return {key: i for i, key in enumerate(keys)}, seq_index, node, edge

    def predict(self, payload, mutations):
        lookup, seq_index, node, edge = payload
        wt, chain, number, insertion, mt = parse_mutations(mutations)[0]
        position = lookup[(chain, number, insertion)]
        if self.dataset.alphabet.index(wt) != int(seq_index[position, 0]):
            raise ValueError('Pythia WT identity differs from experimental mutation')
        if wt == mt:
            raise ValueError('Original Pythia inference excludes identity substitutions')
        wt_ids, mt_ids, nodes, edges, substitutions = self.dataset.inference_extarct_feature(
            wt, position, seq_index, node, edge, device=self.device)
        i = substitutions.index(mt)
        tensors = [value[i].unsqueeze(0).to(self.device) for value in (wt_ids, mt_ids, nodes, edges)]
        with self.torch.inference_mode():
            binding, _ = self.model(*tensors)
        return float(binding.reshape(-1)[0].cpu())


class PpiformerPredictor:
    def __init__(self, repository, device):
        import ppiformer.definitions as definitions
        import torch

        # This checkpoint uses Equiformer. The shared base environment also has
        # an unrelated SE3Transformer, whose presence triggers optional imports.
        definitions.SE3TRANSFORMER_REQUIRED = False
        from ppiformer.tasks.node import DDGPPIformer

        self.torch, self.device = torch, device
        self.repository = repository
        torch.set_num_threads(2)
        # Released Lightning files serialize author module objects in hparams.
        self.models = [DDGPPIformer.load_from_checkpoint(repository/'weights/ddg_regression'/f'{i}.ckpt',
                       map_location='cpu', weights_only=False).to(device).eval() for i in range(3)]

    def prepare(self, record, work):
        from ppiformer.data.transforms import PDBToPyGPretransform

        path = prepare_ppiformer_interface(record, work)
        # The released val pretransform uses k=None and all three vector features.
        return PDBToPyGPretransform(k=None)(path)

    def predict(self, payload, mutations):
        from torch_geometric.data import Batch

        from topoppi_affinity.ppiformer_finetune import labelled_graph

        batch = Batch.from_data_list([labelled_graph(payload, mutations)], follow_batch=['mut_ddg']).to(self.device)
        # Original API always materializes dense attention for all three models,
        # even when return_attn=False. Only scalar predictions are needed here.
        with self.torch.inference_mode():
            values = [model(batch, return_attn=False) for model in self.models]
        return float(self.torch.stack(values).mean().cpu())


def run_mutation_reference(manifest, output, repository, *, method, device='cuda', limit=None):
    """Manifest-wide author inference; resume only per-complex completed outputs."""
    frame = pd.read_csv(manifest, keep_default_na=False)
    if limit is not None:
        frame = frame.iloc[:limit].copy()
    frame = frame.reset_index(drop=True)
    output, repository = Path(output), Path(repository)
    output.mkdir(parents=True, exist_ok=True)
    rows = declared_rows(frame, method)
    predictor = None
    started = time.monotonic()
    classes = {'saambe': SaambePredictor, 'pythia': PythiaPredictor, 'ppiformer': PpiformerPredictor}
    for index, (complex_id, subset) in enumerate(frame.groupby('complex_id', sort=False), 1):
        work = output/'complexes'/complex_id.replace(':', '_')
        work.mkdir(parents=True, exist_ok=True)
        result_path = work/'predictions.csv'
        if result_path.exists():
            cached = pd.read_csv(result_path, keep_default_na=False)
            if set(cached.sample_id) == set(subset.sample_id):
                cached['predicted'] = pd.to_numeric(cached.predicted, errors='coerce')
                aligned = cached.set_index('sample_id').loc[subset.sample_id].reset_index()
                for column in cached:
                    rows.loc[subset.index, column] = aligned[column].to_numpy()
                continue
        eligible = subset.loc[rows.loc[subset.index, 'status'].ne('not_applicable')]
        if len(eligible):
            if predictor is None:
                predictor = classes[method](repository, device)
            try:
                payload = predictor.prepare(eligible.iloc[0].to_dict(), work)
            except Exception as error:
                rows.loc[eligible.index, 'error'] = f'{type(error).__name__}: {error}'
            else:
                for row in eligible.itertuples():
                    try:
                        value = predictor.predict(payload, row.mutations)
                        if not math.isfinite(value):
                            raise ValueError('Nonfinite original model output')
                        rows.loc[row.Index, ['predicted', 'status']] = [value, 'ok']
                    except Exception as error:
                        rows.loc[row.Index, 'error'] = f'{type(error).__name__}: {error}'
        rows.loc[subset.index].to_csv(result_path, index=False)
        print(f'{method}: {index} complexes, {rows.status.eq("ok").sum()} successful rows', flush=True)
    rows.to_csv(output/'predictions.csv', index=False)
    metadata = {**METHODS[method], 'repository': str(repository), 'python': sys.executable,
                'declared': len(rows), 'counts': rows.status.value_counts().to_dict(),
                'seconds': time.monotonic()-started, 'targets_read': False,
                'native_target': 'binding MT-WT ddG; positive destabilizing',
                'protocol': {'saambe': 'Original file_loop and regression_v01.model; numbered canonical residues; source REMARK inputs preserved',
                             'pythia': 'Original binding checkpoint and inference feature function; native complete-backbone mapping; selected manifest substitutions only',
                             'ppiformer': 'Original three-checkpoint ensemble; author PPIRef joined 10A interfaces on manifest chains; no missing-site mean imputation'}[method],
                'multichain_note': 'Original chain-based features include contacts within a multi-chain partner; not a new partner-aware model',
                'training_exposure': 'Released supervised reference; do not interpret as unseen local-fold training'}
    (output/'provenance.json').write_text(json.dumps(metadata, indent=2)+'\n')
    return metadata
