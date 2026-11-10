"""DDAffinity / GearBind references using released structures, never FoldX jobs.

Author checkpoints and published out-of-fold predictions retain their original
supervision scope. They are not local-split retraining results.
"""
from __future__ import annotations

import copy
import json
import math
import pickle
import random
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from topoppi_affinity.data import canonical_mutations, complex_key


def variant_key(pdb, first, second, mutations):
    # The author preprocessing uses set(mutstr.split(',')); a few source rows
    # repeat the same substitution verbatim.
    mutations = ','.join(dict.fromkeys(mutations.split(',')))
    return complex_key(pdb, list(first), list(second)), canonical_mutations(mutations)


def record_key(record):
    return variant_key(record['pdb_id'], record['chains_a'].split(','),
                       record['chains_b'].split(','), record['mutations'])


class ReleasedCache:
    """Read the author's existing LMDB. No preprocessing constructor is called."""

    def __init__(self, directory):
        import lmdb

        self.directory = Path(directory)
        with (self.directory/'entries.pkl').open('rb') as stream:
            self.entries = pickle.load(stream)
        self.index = defaultdict(list)
        for i, entry in enumerate(self.entries):
            key = variant_key(entry['complex'], entry['group_ligand'],
                              entry['group_receptor'], entry['mutstr'])
            self.index[key].append(i)
        self.db = lmdb.open(str(self.directory/'structures.lmdb'), subdir=False,
                            readonly=True, create=False, lock=False, readahead=False)

    def resolve(self, record):
        indices = self.index.get(record_key(record), [])
        if not indices:
            raise KeyError('No author cache entry for this partner pair and cleaned mutation set')
        # Replicate measurements may have separate FoldX realizations. Select the
        # first original author row, never the closest experimental target.
        return min(indices, key=lambda i: self.entries[i]['id']), len(indices)

    def pair(self, index):
        entry = self.entries[index]
        with self.db.begin() as txn:
            return tuple(pickle.loads(txn.get((state+'_'+entry['pdbcode']).encode()))[0]
                         for state in ('wt', 'mt'))

    def close(self):
        self.db.close()


def cache_coverage(manifest, cache_dir, output):
    frame = pd.read_csv(manifest, keep_default_na=False)
    cache = ReleasedCache(cache_dir)
    rows = []
    try:
        for record in frame.to_dict('records'):
            row = {'sample_id': record['sample_id'], 'status': 'error', 'error': ''}
            try:
                index, matches = cache.resolve(record)
                entry = cache.entries[index]
                row.update(status='ok', cache_entry=entry['pdbcode'], cache_index=index,
                           author_replicates=matches, n_mutations=entry['num_muts'])
            except (KeyError, ValueError) as error:
                row['error'] = str(error)
            rows.append(row)
    finally:
        cache.close()
    result = pd.DataFrame(rows)
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)
    return result


class DDAffinityPredictor:
    def __init__(self, repository, cache, checkpoint, device='cuda'):
        import torch

        sys.path.insert(0, str(repository))
        from rde.datasets.skempi_parallel import SkempiDataset_lmdb
        from rde.models.protein_mpnn_network_2 import ProteinMPNN_NET
        from rde.utils.data_skempi_mpnn import PaddingCollate
        from rde.utils.transforms import get_transform

        self.torch, self.device, self.cache = torch, device, cache
        ckpt = torch.load(checkpoint, map_location='cpu', weights_only=False)
        self.models = []
        for state in ckpt['model']['models']:
            model = ProteinMPNN_NET(ckpt['config'].model)
            model.load_state_dict(state, strict=True)
            self.models.append(model.to(device).eval())
        complexes = sorted({entry['complex'] for entry in cache.entries})
        # SkempiDatasetManager does not pass config.train.seed; the released
        # dataset constructor therefore uses split_seed=2022.
        random.Random(2022).shuffle(complexes)
        size = math.ceil(len(complexes)/len(self.models))
        self.author_fold = {name: i//size for i, name in enumerate(complexes)}
        # Preserve the original item transform while excluding the constructor's
        # cache-generation branch and all true affinity labels.
        self.dataset = SkempiDataset_lmdb.__new__(SkempiDataset_lmdb)
        self.dataset.entries = [dict(entry, ddG=np.float32(0)) for entry in cache.entries]
        self.dataset.cache_dir = str(cache.directory)
        self.dataset.db_conn, self.dataset.db_keys = cache.db, None
        self.dataset.transform = get_transform(ckpt['config'].data.val.transform)
        self.collate = PaddingCollate(ckpt['config'].model.patch_size)

    def predict(self, index, record):
        random.seed(2024 + int(self.cache.entries[index]['id']))
        batch = self.collate([self.dataset[index]])

        def move(value):
            if isinstance(value, self.torch.Tensor):
                return value.to(self.device)
            if isinstance(value, dict):
                return {k: move(v) for k, v in value.items()}
            return value

        batch = move(batch)
        fold = self.author_fold[self.cache.entries[index]['complex']]
        with self.torch.inference_mode():
            predicted = float(self.models[fold](batch)[1]['ddG_pred'].reshape(-1)[0].cpu())
        return predicted, {'ensemble_size': 1, 'author_fold': fold}


def run_cached_reference(manifest, output, repository, cache_dir, checkpoint, *,
                         method, device='cuda', limit=None, dda_repository=None):
    frame = pd.read_csv(manifest, keep_default_na=False)
    if limit is not None:
        frame = frame.iloc[:limit]
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    prediction_path = output/'predictions.csv'
    names = {'ddaffinity': 'DDAffinity-author-fold-cache',
             'gearbind': 'GearBind-author-DDAffinity-cache'}
    cached = {}
    if prediction_path.exists():
        old = pd.read_csv(prediction_path, keep_default_na=False)
        cached = {r['sample_id']: r for r in old.to_dict('records')
                  if r['status'] == 'ok' and r['method'] == names[method]}
    cache = ReleasedCache(cache_dir)
    rows, predictor = [], None
    started = time.monotonic()
    try:
        for record in frame.to_dict('records'):
            sid = record['sample_id']
            if sid in cached:
                rows.append(cached[sid])
                continue
            row = {'sample_id': sid, 'complex_id': record['complex_id'],
                   'method': names[method], 'target': 'ddg', 'predicted': np.nan,
                   'status': 'error', 'error': '', 'comparison_scope': 'author_supervised_reference',
                   'training_overlap': ('released author PDB-fold; not local family-isolated retraining'
                                        if method == 'ddaffinity' else
                                        'SKEMPI supervision; author ensemble is not local-fold independent')}
            try:
                index, matches = cache.resolve(record)
                row.update(cache_entry=cache.entries[index]['pdbcode'], author_replicates=matches)
                if predictor is None:
                    if method == 'ddaffinity':
                        predictor = DDAffinityPredictor(repository, cache, checkpoint, device)
                    else:
                        predictor = GearBindPredictor(repository, dda_repository, cache, checkpoint, device)
                predicted, details = predictor.predict(index, record)
                if not math.isfinite(predicted):
                    raise ValueError('Nonfinite original-model prediction')
                row.update(predicted=predicted, status='ok', **details)
            except Exception as error:
                row['error'] = f'{type(error).__name__}: {error}'
            rows.append(row)
            # Predictions are expensive; retain completed work if interrupted.
            if len(rows) % 25 == 0:
                pd.DataFrame(rows).to_csv(prediction_path, index=False)
    finally:
        cache.close()
    result = pd.DataFrame(rows)
    result.to_csv(prediction_path, index=False)
    (output/'provenance.json').write_text(json.dumps({
        'method': names[method], 'declared': len(result),
        'success': int(result.status.eq('ok').sum()), 'seconds': time.monotonic()-started,
        'repository': str(repository), 'checkpoint': str(checkpoint), 'cache': str(cache_dir),
        'source_url': 'https://github.com/'+('ak422/DDAffinity' if method == 'ddaffinity' else 'DeepGraphLearning/GearBind'),
        'structure_source': 'DDAffinity released SKEMPI2_cache, FoldX5-generated WT/MT heavy-atom tensors',
        'foldx_executed': False, 'labels_provided_to_model': 'zero placeholders only',
        'replicate_policy': 'first author row by id; no label matching or calibration',
        'weight_selection': ('original author PDB-fold, split_seed=2022, 10 folds'
                             if method == 'ddaffinity' else 'all five released supervised folds; exposed on SKEMPI'),
        'input_protocol': ('original LMDB and validation transform; fixed per-entry chain shuffle'
                           if method == 'ddaffinity' else
                           'author GearBind graphs built from retained DDAffinity cached atoms; chain identity recovered from original WT coordinates'),
    }, indent=2)+'\n')
    return result


def import_gearbind_table(manifest, source, output):
    """Published author-fold values, explicitly separate from local inference."""
    frame = pd.read_csv(manifest, keep_default_na=False)
    table = pd.read_csv(source, keep_default_na=False)
    lookup = defaultdict(list)
    for row in table.to_dict('records'):
        pdb, first, second = row['#Pdb'].split('_')
        lookup[variant_key(pdb, first, second, row['Mutation(s)_cleaned'])].append(row)
    rows = []
    for record in frame.to_dict('records'):
        matches = lookup.get(record_key(record), [])
        valid = [r for r in matches if pd.notna(pd.to_numeric(r['GearBind'], errors='coerce'))]
        rows.append({'sample_id': record['sample_id'], 'method': 'GearBind-published-author-fold',
                     'predicted': np.mean([float(r['GearBind']) for r in valid]) if valid else np.nan,
                     'target': 'ddg', 'status': 'ok' if valid else 'error',
                     'error': '' if valid else 'No finite published value for exact partner/mutation key',
                     'comparison_scope': 'published_author_split_reference',
                     'author_fold': ','.join(sorted({str(r['fold']) for r in valid})),
                     'author_replicates': len(valid), 'local_model_executed': False})
    result = pd.DataFrame(rows)
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)
    Path(output).with_suffix('.json').write_text(json.dumps({
        'method': 'GearBind-published-author-fold', 'source': str(source),
        'source_url': 'https://github.com/DeepGraphLearning/GearBind',
        'column': 'GearBind', 'target': 'MT-WT ddG in kcal/mol',
        'matching': 'PDB, unordered partner chain sets, cleaned mutation set',
        'replicates': 'arithmetic mean of published values for identical keys',
        'local_model_executed': False, 'local_fold_independence_claimed': False,
    }, indent=2)+'\n')
    return result


class GearBindPredictor:
    def __init__(self, repository, dda_repository, cache, checkpoint_dir, device='cuda'):
        import torch
        import yaml

        sys.path[:0] = [str(repository), str(dda_repository)]
        from gearbind import dataset, layer, model, task  # noqa: F401
        from torchdrug import core

        self.torch, self.device, self.cache = torch, device, cache
        config = yaml.safe_load((Path(repository)/'config/predict/HER2_GearBind.yaml').read_text())
        self.models = []
        for i in range(5):
            net = core.Configurable.load_config_dict(copy.deepcopy(config['task']))
            state = torch.load(Path(checkpoint_dir)/f'gearbind{i}.pth', map_location='cpu', weights_only=False)['model']
            for name in ('mean', 'std', 'weight'):
                net.register_buffer(name, state[name].clone())
            net.load_state_dict(state, strict=True)
            self.models.append(net.to(device).eval())

    def predict(self, index, record):
        from torchdrug import data

        wt, mt = self.cache.pair(index)
        mapping = recover_cached_chains(wt, record['structure_path'])
        validate_cached_mutations(wt, mt, mapping, record['mutations'])
        proteins = [cached_gearbind_protein(part, mapping, record) for part in (wt, mt)]
        # Original graph construction modifies the input dictionary; pack anew
        # for each ensemble member rather than transforming a transformed graph.
        predictions = []
        with self.torch.inference_mode():
            for net in self.models:
                batch = {name: data.Protein.pack([protein]).to(self.device)
                         for name, protein in zip(('wild_type', 'mutant'), proteins, strict=True)}
                predictions.append(float(net.predict(batch).reshape(-1)[0].cpu()))
        return float(np.mean(predictions)), {'ensemble_size': 5, 'cache_chain_mapping': json.dumps(mapping)}


def recover_cached_chains(chains, original_pdb):
    """DDAffinity discarded chain letters. Recover them from unchanged WT frame.

    Assignment uses C-alpha coordinates and residue identity, not affinity labels.
    Ambiguous copies or a transformed coordinate frame are reported as failures.
    """
    from Bio.PDB import PDBParser
    from rde.utils.protein.constants import AA
    from scipy.optimize import linear_sum_assignment
    from scipy.spatial import cKDTree

    original = next(PDBParser(QUIET=True).get_structure('reference', original_pdb).get_models())
    original = [chain for chain in original if any('CA' in r for r in chain)]
    chain_ids = list(chains)
    costs = np.full((len(chain_ids), len(original)), 1e6)
    for i, key in enumerate(chain_ids):
        cached = chains[key]
        xyz = cached['pos_heavyatom'][:, 1].numpy()
        types = cached['aa'].numpy()
        for j, chain in enumerate(original):
            residues = [r for r in chain if 'CA' in r and AA.is_aa(r.resname)]
            if not residues:
                continue
            coords = np.stack([r['CA'].coord for r in residues])
            distance, nearest = cKDTree(coords).query(xyz)
            identities = np.array([int(AA(r.resname)) for r in residues])
            if len(set(nearest.tolist())) == len(nearest) and np.array_equal(types, identities[nearest]):
                costs[i, j] = np.sqrt(np.mean(distance**2))
    source, destination = linear_sum_assignment(costs)
    if len(source) != len(chain_ids) or (costs[source, destination] > 2.0).any():
        raise ValueError('Cannot recover cached chain identities from WT frame within 2 Angstrom RMSD')
    for i, j in zip(source, destination, strict=True):
        if np.count_nonzero(np.isclose(costs[i], costs[i, j], atol=0.01)) != 1:
            raise ValueError('Ambiguous chain identity in author cache')
    return {int(chain_ids[i]): original[j].id for i, j in zip(source, destination, strict=True)}


def validate_cached_mutations(wt, mt, mapping, mutations):
    from Bio.Data.PDBData import protein_letters_3to1
    from rde.utils.protein.constants import AA

    from topoppi_affinity.data import parse_mutations

    observed = set()
    for key in wt:
        before, after = wt[key], mt[key]
        if not np.array_equal(before['resseq'], after['resseq']):
            raise ValueError('WT/MT cache residue correspondence differs')
        for i in (before['aa'] != after['aa']).nonzero().flatten().tolist():
            observed.add((protein_letters_3to1[AA(int(before['aa'][i])).name], mapping[key],
                          int(before['resseq'][i]), '', protein_letters_3to1[AA(int(after['aa'][i])).name]))
    if observed != set(parse_mutations(mutations)):
        raise ValueError('Cached WT/MT substitutions do not match requested chain and positions')


def cached_gearbind_protein(chains, mapping, record):
    """Repack retained atoms without relaxation, atom completion or PDB rounding."""
    import torch
    from rde.utils.protein.constants import AA, restype_to_heavyatom_names
    from torch.nn import functional as F
    from torchdrug import data

    from topoppi_affinity.data import parse_mutations

    a, b = set(record['chains_a'].split(',')), set(record['chains_b'].split(','))
    positions, atom_types, atom_names, atom_to_residue = [], [], [], []
    residue_types, residue_numbers, chain_names, entity_a, entity_b, mutation = [], [], [], [], [], []
    mutation_sites = {(chain, number) for _, chain, number, _, _ in parse_mutations(record['mutations'])}
    for key in sorted(chains):
        chain = mapping[key]
        if chain not in a | b:
            continue
        part = chains[key]
        for i, aa in enumerate(part['aa'].tolist()):
            resname = AA(aa).name
            number = int(part['resseq'][i])
            residue_index = len(residue_types)
            residue_types.append(data.Protein.residue2id[resname])
            residue_numbers.append(number)
            chain_names.append(data.Protein.alphabet2id[chain])
            entity_a.append(chain in a)
            entity_b.append(chain in b)
            mutation.append((chain, number) in mutation_sites)
            for j, name in enumerate(restype_to_heavyatom_names[AA(aa)]):
                if name and part['mask_heavyatom'][i, j] and name in data.Protein.atom_name2id:
                    positions.append(part['pos_heavyatom'][i, j])
                    atom_types.append(data.feature.atom_vocab.get(name[0], 0))
                    atom_names.append(data.Protein.atom_name2id[name])
                    atom_to_residue.append(residue_index)
    if sum(mutation) != len(mutation_sites):
        raise ValueError('Cached GearBind mutation positions missing or duplicated')
    protein = data.Protein([[0, 0, 0]], atom_type=atom_types, bond_type=[0],
                           residue_type=residue_types, atom_name=atom_names, atom2residue=atom_to_residue,
                           residue_number=residue_numbers, chain_id=chain_names,
                           insertion_code=[-1]*len(residue_types), node_position=torch.stack(positions),
                           num_node=len(positions), num_residue=len(residue_types))
    residue_feature = F.one_hot(protein.residue_type, len(data.Protein.residue2id)+1)
    with protein.node():
        protein.entity_a = torch.tensor(entity_a)[protein.atom2residue]
        protein.entity_b = torch.tensor(entity_b)[protein.atom2residue]
        protein.is_mutation = torch.tensor(mutation)[protein.atom2residue]
        protein.node_feature = torch.cat([F.one_hot(protein.atom_name, 37), residue_feature[protein.atom2residue]], dim=-1)
    with protein.residue():
        protein.residue_feature = residue_feature
    return protein
