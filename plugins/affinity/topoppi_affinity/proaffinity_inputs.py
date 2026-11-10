"""FASTA-aware ProAffinity-GNN inputs using the released preprocessing rules."""
from __future__ import annotations

import ast
import io
import json
import pickle
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd


def upstream_functions(path):
    tree = ast.parse(Path(path).read_text())
    nodes = [n for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef))
             or isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'aminoacid_abbr' for t in n.targets)]
    namespace = {}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), namespace)
    return namespace


def read_fasta_chains(text):
    """RCSB entity FASTA names auth-chain identifiers explicitly when different."""
    from Bio import SeqIO

    result = {}
    for record in SeqIO.parse(io.StringIO(text), 'fasta'):
        group = record.description.split('|')[1]
        group = re.sub(r'^Chains?\s+', '', group)
        for token in group.split(','):
            match = re.search(r'\[auth\s+([^\]]+)\]', token)
            chain = match.group(1).strip() if match else token.strip()
            result[chain] = str(record.seq)
    return result


def prefetch_fastas(pdb_ids, directory):
    """Overlap independent RCSB downloads; GPU graph construction stays serial."""
    import requests

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    pending = [pdb for pdb in dict.fromkeys(str(p).lower() for p in pdb_ids)
               if not (directory/(pdb+'.fasta')).exists()]

    def fetch(pdb):
        destination = directory/(pdb+'.fasta')
        try:
            response = requests.get(f'https://www.rcsb.org/fasta/entry/{pdb.upper()}/display', timeout=60)
            response.raise_for_status()
            temporary = destination.with_suffix('.fasta.download')
            temporary.write_text(response.text)
            temporary.replace(destination)
            return None
        except requests.RequestException as error:
            return {'pdb_id': pdb, 'error': str(error)}

    with ThreadPoolExecutor(max_workers=4) as pool:
        failures = [result for result in pool.map(fetch, pending) if result is not None]
    return {'requested': len(pending), 'downloaded': len(pending)-len(failures), 'failures': failures}


def align_observed_residues(residues, full_sequence):
    # Same greedy subsequence convention as match_seq_with_fasta, with the
    # missing end-of-FASTA bound fixed for a complete observed chain.
    cursor = 0
    for residue in residues:
        position = full_sequence.find(residue['type'], cursor)
        if position < 0:
            raise ValueError('Observed PDBQT chain does not match supplied complete FASTA')
        residue['number'] = position
        cursor = position + 1


def atom_histograms(pairs, atom_pairs, distance_scale):
    """Vectorized original 28 atom-pair classes × 10 distance bins."""
    import torch

    lookup = {}
    for i, pair in enumerate(atom_pairs):
        a, b = pair.split('_')
        lookup[(a, b)] = lookup[(b, a)] = i
    output = []
    for first, second, _ in pairs:
        atoms_a, atoms_b = first['atoms'], second['atoms']
        xyz_a = np.array([[float(a[k]) for k in ('x', 'y', 'z')] for a in atoms_a])
        xyz_b = np.array([[float(a[k]) for k in ('x', 'y', 'z')] for a in atoms_b])
        distances = np.sqrt(np.square(xyz_a[:, None]-xyz_b[None, :]).sum(-1))
        bins = np.minimum(np.ceil(distances/(distance_scale/10)).astype(int), 10)-1
        classes = np.array([[lookup[(a['pdbqt_type'], b['pdbqt_type'])] for b in atoms_b] for a in atoms_a])
        # Negative bin at exactly zero follows the author's numpy indexing.
        indices = (bins*len(atom_pairs)+classes).reshape(-1) % (10*len(atom_pairs))
        output.append(np.bincount(indices, minlength=10*len(atom_pairs)))
    return torch.tensor(np.stack(output), dtype=torch.float64) if output else torch.empty((0, 10*len(atom_pairs)), dtype=torch.float64)


def interaction_pairs(first, second, radius, *, intra=False):
    """Exact original any-atom <= radius criterion, preserving nested-loop order."""
    from scipy.spatial import cKDTree

    def coordinates(residues):
        xyz, owners = [], []
        for i, residue in enumerate(residues):
            for atom in residue['atoms']:
                xyz.append([float(atom[k]) for k in ('x', 'y', 'z')])
                owners.append(i)
        return np.asarray(xyz), np.asarray(owners)

    xyz_a, owner_a = coordinates(first)
    xyz_b, owner_b = coordinates(second)
    neighbors = cKDTree(xyz_a).sparse_distance_matrix(cKDTree(xyz_b), radius, output_type='coo_matrix')
    residue_pairs = sorted(set(zip(owner_a[neighbors.row].tolist(), owner_b[neighbors.col].tolist(), strict=True)))
    pairs = []
    for i, j in residue_pairs:
        a, b = first[i], second[j]
        if intra and a['number'] == b['number']:
            continue
        ca_a = next(atom for atom in a['atoms'] if atom['type'] == 'CA')
        ca_b = next(atom for atom in b['atoms'] if atom['type'] == 'CA')
        distance = sum((float(ca_a[k])-float(ca_b[k]))**2 for k in ('x', 'y', 'z'))**.5
        pairs.append((a, b, distance))
    return pairs


def graph_inputs(repository, pdbqt, chain_groups, sequences, embed):
    import torch
    from torch_geometric.data import Data

    inter = upstream_functions(Path(repository)/'pdb2graph.py')
    intra = upstream_functions(Path(repository)/'pdb2graph_individual.py')
    source = ast.parse((Path(repository)/'graph_construct.py').read_text())
    atom_pairs = next(ast.literal_eval(n.value) for n in source.body
                      if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'atom_pair' for t in n.targets))
    residues, features = [], []
    for group in chain_groups:
        parts = inter['get_residue_list_from_file'](str(pdbqt), group)
        if len(parts) != len(group) or any(not part for part in parts):
            raise ValueError('PDBQT lacks a requested chain with C-alpha atoms')
        sequence_list = [sequences[chain] for chain in group]
        for part, sequence in zip(parts, sequence_list, strict=True):
            align_observed_residues(part, sequence)
        residues.append(inter['concat_reslist'](parts, list(range(len(group))), sequence_list))
        features.append(torch.cat([embed(sequence) for sequence in sequence_list], dim=0))
    pairs = interaction_pairs(residues[0], residues[1], 15)
    edges = inter['get_edge_index'](pairs, len(features[0])).reshape(2, -1).long()
    hist = atom_histograms(pairs, atom_pairs, 15)
    graphs = [Data(x=torch.cat(features), edge_index=edges, edge_attr=torch.cat([hist, hist]))]
    for part, x in zip(residues, features, strict=True):
        pairs = interaction_pairs(part, part, 3.5, intra=True)
        edges = intra['get_edge_index'](pairs).reshape(2, -1).long()
        graphs.append(Data(x=x, edge_index=edges, edge_attr=atom_histograms(pairs, atom_pairs, 3.5)))
    return graphs


def normalize_pdbqt_terminators(source, output):
    # ADFR emits ATOM-only PDBQT. The original parser flushes its last residue
    # only at TER, so supply one per chain (coordinates and atom types unchanged).
    by_chain = {}
    for line in Path(source).read_text().splitlines():
        if line.startswith('ATOM'):
            by_chain.setdefault(line[21], []).append(line)
    lines = []
    for chain, atoms in by_chain.items():
        lines.extend(atoms)
        lines.append('TER'.ljust(21)+chain)
    Path(output).write_text('\n'.join(lines)+'\n')


def write_proaffinity_pdb(record, output):
    """Keep original fixed-column atoms for the selected first-model partners.

    The installed BioPython writer shifts the element columns left by one.
    ADFR then reads calcium as A and selenium as E. Selecting the original
    records preserves coordinates, atom names and elements without reformatting.
    """
    chains = set(record['chains_a'].split(',') + record['chains_b'].split(','))
    selected = []
    for line in Path(record['structure_path']).read_text().splitlines():
        if line.startswith('ENDMDL'):
            break
        if line.startswith(('ATOM  ', 'HETATM', 'TER   ')) and line[21:22] in chains:
            selected.append(line)
    Path(output).write_text('\n'.join(selected) + '\nEND\n')


def prepare_pdbqt(record, output, prepare_receptor, pdbqt_dir=None):
    output = Path(output)
    key = record['pdb_id'].lower()+'_'+record['chains_a'].replace(',', '')+'_'+record['chains_b'].replace(',', '')
    pdbqt = output/'pdbqt'/(key+'.pdbqt')
    if pdbqt.exists():
        return pdbqt
    supplied = Path(pdbqt_dir)/(record['pdb_id'].lower()+'_atom_processed.pdbqt') if pdbqt_dir else None
    if supplied and supplied.exists():
        raw = supplied
    else:
        pdb = output/'pdb'/(key+'.input.pdb')
        write_proaffinity_pdb(record, pdb)
        raw = output/'pdbqt'/(key+'.preparation.raw.pdbqt')
        with (output/'logs'/(key+'.preparation.log')).open('w') as log:
            subprocess.run([str(prepare_receptor), '-r', str(pdb), '-A', 'hydrogens', '-o', str(raw)],
                           stdout=log, stderr=subprocess.STDOUT, check=True)
    temporary = pdbqt.with_suffix('.pdbqt.preparing')
    normalize_pdbqt_terminators(raw, temporary)
    temporary.replace(pdbqt)
    return pdbqt


def prefetch_pdbqts(records, output, prepare_receptor, pdbqt_dir=None):
    """Two independent ADFR processes; deduplicate repeated PPB structure inputs."""
    records = list({(r['pdb_id'].lower(), r['chains_a'], r['chains_b']): r for r in records}.values())

    def prepare(record):
        try:
            prepare_pdbqt(record, output, prepare_receptor, pdbqt_dir)
            return None
        except Exception as error:
            return {'pdb_id': record['pdb_id'], 'chains_a': record['chains_a'], 'chains_b': record['chains_b'],
                    'error': f'{type(error).__name__}: {error}'}

    with ThreadPoolExecutor(max_workers=2) as pool:
        failures = [result for result in pool.map(prepare, records) if result is not None]
    return {'requested': len(records), 'prepared': len(records)-len(failures), 'failures': failures}


def prepare_proaffinity(manifest, output, repository, esm_checkpoint, *,
                        prepare_receptor, fasta_dir=None, pdbqt_dir=None, device='cuda', limit=None):
    import requests
    import torch
    from transformers import AutoTokenizer, EsmModel

    frame = pd.read_csv(manifest, keep_default_na=False)
    if limit is not None:
        frame = frame.iloc[:limit].copy()
    output = Path(output)
    for child in ('graph/inter_graph', 'graph/individual_graph', 'pdbqt', 'fasta', 'pdb', 'logs'):
        (output/child).mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    fasta_downloads = prefetch_fastas(frame.pdb_id, Path(fasta_dir) if fasta_dir else output/'fasta')
    pdbqt_preparation = prefetch_pdbqts(frame.to_dict('records'), output, prepare_receptor, pdbqt_dir)
    rows, partner_lines = [], {}
    tokenizer, encoder = None, None

    def embed(sequence):
        nonlocal tokenizer, encoder
        if encoder is None:
            tokenizer = AutoTokenizer.from_pretrained(str(esm_checkpoint), local_files_only=True)
            encoder = EsmModel.from_pretrained(str(esm_checkpoint), add_pooling_layer=False, local_files_only=True).to(device).eval()
        with torch.inference_mode():
            tokens = tokenizer(sequence, return_tensors='pt').to(device)
            return encoder(**tokens).last_hidden_state[0, 1:-1].cpu()

    for record in frame.to_dict('records'):
        row_started = time.monotonic()
        groups = [record['chains_a'].split(','), record['chains_b'].split(',')]
        key = record['pdb_id'].lower()+'_'+'_'.join(''.join(g) for g in groups)
        partner_lines[key] = key+'\t'+'; '.join(','.join(group) for group in groups)+';\n'
        paths = [output/'graph/inter_graph'/key,
                 output/'graph/individual_graph'/(key+'_1'), output/'graph/individual_graph'/(key+'_2')]
        reused = all(path.exists() for path in paths)
        row = {'sample_id': record['sample_id'], 'author_graph_key': key, 'status': 'error', 'error': '',
               'reused_graph': reused}
        try:
            if not reused:
                fasta = (Path(fasta_dir) if fasta_dir else output/'fasta')/(record['pdb_id'].lower()+'.fasta')
                if not fasta.exists():
                    url = f"https://www.rcsb.org/fasta/entry/{record['pdb_id'].upper()}/display"
                    response = requests.get(url, timeout=60)
                    response.raise_for_status()
                    fasta.parent.mkdir(parents=True, exist_ok=True)
                    fasta.write_text(response.text)
                sequences = read_fasta_chains(fasta.read_text())
                pdbqt = prepare_pdbqt(record, output, prepare_receptor, pdbqt_dir)
                graphs = graph_inputs(repository, pdbqt, groups, sequences, embed)
                for path, graph in zip(paths, graphs, strict=True):
                    with path.open('wb') as stream:
                        pickle.dump(graph, stream)
            row['status'] = 'ok'
        except Exception as error:
            row['error'] = f'{type(error).__name__}: {error}'
        row['seconds'] = time.monotonic() - row_started
        rows.append(row)
        pd.DataFrame(rows).to_csv(output/'coverage.csv', index=False)
    coverage = pd.DataFrame(rows)
    frame = frame.drop(columns=['author_graph_key'], errors='ignore').merge(
        coverage[['sample_id', 'author_graph_key']], on='sample_id', validate='one_to_one')
    frame.to_csv(output/'manifest.csv', index=False)
    (output/'partners.txt').write_text(''.join(partner_lines.values()))
    (output/'provenance.json').write_text(json.dumps({
        'source_url': 'https://github.com/legendzzy/ProAffinity-GNN',
        'protocol': 'full FASTA, original inter/intra pair selection and 280-bin edges, original ESM2 residue states',
        'runtime_changes': ['KD-tree equivalent atom-contact selection', 'vectorized atom-pair histogram', 'complete-FASTA end bound',
                            'TER lines for ADFR ATOM-only output', 'ESM2 without unused pooler; no gradients',
                            'Original fixed-column PDB records preserve elements for ADFR'],
        'fasta_source': str(fasta_dir) if fasta_dir else 'RCSB entity FASTA; auth chain identifiers',
        'esm_checkpoint': str(esm_checkpoint), 'prepare_receptor': str(prepare_receptor),
        'fasta_downloads': fasta_downloads,
        'pdbqt_preparation': pdbqt_preparation,
        'declared': len(frame), 'success': int(coverage.status.eq('ok').sum()),
        'seconds': time.monotonic()-started, 'foldx_executed': False,
    }, indent=2)+'\n')
    return coverage
