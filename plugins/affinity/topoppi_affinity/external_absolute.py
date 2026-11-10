"""Run released absolute-affinity methods, preserving upstream models and units."""
from __future__ import annotations

import json
import math
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd

MIRAGE_REVISION = '01e872fad77ded61f05a514068fcfddbf1219e66'


def download_external_sources(method: str, source_root: Path) -> None:
    """Acquire original small source/metadata and, for MIRAGE, author ProtT5."""
    import io
    import urllib.request
    import zipfile

    source_root.mkdir(parents=True, exist_ok=True)

    def fetch(url: str, destination: Path):
        if destination.exists():
            return
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_suffix(destination.suffix + '.part')
        with urllib.request.urlopen(url, timeout=120) as response, temporary.open('wb') as stream:
            while chunk := response.read(4 * 1024 * 1024):
                stream.write(chunk)
        temporary.replace(destination)

    if method == 'prodigy':
        target = source_root / 'prodigy' / 'elife-07454-supp2.csv'
        if target.exists():
            return
        target.parent.mkdir(parents=True, exist_ok=True)
        url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4523921/supplementaryFiles'
        with urllib.request.urlopen(url, timeout=120) as response:
            archive = zipfile.ZipFile(io.BytesIO(response.read()))
        target.write_bytes(archive.read('elife-07454-supp2.csv'))
    elif method == 'mirage':
        repository = source_root / 'MIRAGE_author'
        archive = source_root / 'MIRAGE_author.zip'
        fetch(f'https://codeload.github.com/ShiweiWu-545/MIRAGE/zip/{MIRAGE_REVISION}', archive)
        if not (repository / 'data' / 'model.pt').is_file():
            with zipfile.ZipFile(archive) as source:
                for member in source.infolist():
                    parts = Path(member.filename).parts[1:]
                    if member.is_dir() or not parts:
                        continue
                    destination = repository.joinpath(*parts)
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    destination.write_bytes(source.read(member))
        prottrans = source_root.parent / 'downloads' / 'prottrans.zip'
        fetch('https://zenodo.org/records/17563574/files/prottrans.zip?download=1', prottrans)
        if not (repository / 'prottrans' / 'prott5_embedder.py').is_file():
            with zipfile.ZipFile(prottrans) as source:
                source.extractall(repository)
    else:
        raise ValueError(f'Unknown released method: {method}')

def native_dg_to_pkd(dg_kcal_mol: float, temperature_K: float = 298.15,
                     *, positive_affinity_score: bool = False,
                     gas_constant: float = 0.001987) -> float:
    """Convert documented energy conventions without fitting a calibration."""
    if not math.isfinite(dg_kcal_mol) or temperature_K <= 0:
        raise ValueError('Prediction must be finite and conversion temperature positive.')
    sign = 1 if positive_affinity_score else -1
    return sign * dg_kcal_mol / (gas_constant * temperature_K * math.log(10))


def select_records(frame: pd.DataFrame, subset: str = 'all', limit: int | None = None) -> pd.DataFrame:
    if subset == 'heldout':
        frame = frame[(frame.split_complex == 'test') | (frame.split_sequence == 'test')]
    elif subset != 'all':
        raise ValueError(f'Unsupported subset: {subset}')
    frame = frame.sort_values('sample_id').copy()
    return frame.iloc[:limit] if limit is not None else frame


def write_partner_structure(record: dict, destination: Path) -> None:
    """Select first-model partner groups from the supplied manifest."""
    from Bio.PDB import PDBIO, PDBParser, Select
    selected = set(record['chains_a'].split(',') + record['chains_b'].split(','))
    model = next(PDBParser(QUIET=True).get_structure(record['pdb_id'], record['structure_path']).get_models())
    if not selected.issubset({chain.id for chain in model}):
        raise ValueError('Manifest partner chains are absent from structure.')

    class Partners(Select):
        def accept_chain(self, chain):
            return chain.id in selected

    destination.parent.mkdir(parents=True, exist_ok=True)
    writer = PDBIO()
    writer.set_structure(model)
    writer.save(str(destination), Partners())


def read_prodigy_training_ids(path: Path) -> set[str]:
    """Read original publication Supplementary file 2, preserving its PDB IDs."""
    if not path.is_file():
        return set()
    import re
    # The publication CSV contains descriptive header rows; PDB IDs begin rows.
    ids = set()
    for line in path.read_text(errors='replace').splitlines():
        match = re.match(r'^\s*"?([0-9][A-Za-z0-9]{3})(?:[._,;\s"]|$)', line)
        if match:
            ids.add(match[1].lower())
    return ids


def _prodigy_record(arguments) -> dict:
    record, work_dir, training_ids = arguments
    started = time.monotonic()
    output = {'sample_id': record['sample_id'], 'predicted': np.nan, 'target': 'pkd',
              'method': 'PRODIGY', 'status': 'error', 'error': '',
              'training_overlap': 'known' if record['pdb_id'].lower() in training_ids else 'possible',
              'exact_training_pdb_overlap': record['pdb_id'].lower() in training_ids,
              'reference_temperature_K': 298.15}
    try:
        from prodigy_prot.modules.parsers import parse_structure
        from prodigy_prot.modules.prodigy import Prodigy
        path = Path(work_dir) / record['pdb_id'] / 'partners.pdb'
        write_partner_structure(record, path)
        models, _, _ = parse_structure(str(path))
        predictor = Prodigy(models[0], name=record['pdb_id'],
                            selection=[record['chains_a'], record['chains_b']], temp=25.0)
        predictor.predict(distance_cutoff=5.5, acc_threshold=0.05)
        output.update(predicted=-math.log10(predictor.kd_val),
                      native_dg_kcal_mol=float(predictor.ba_val), native_kd_M=float(predictor.kd_val), status='ok')
        (path.parent / 'native_prediction.json').write_text(json.dumps(predictor.as_dict(), indent=2) + '\n')
    except Exception as error:
        output['error'] = f'{type(error).__name__}: {error}'
    output['seconds'] = round(time.monotonic() - started, 3)
    return output


def run_prodigy(manifest: Path, output: Path, source_root: Path, subset: str = 'all',
                limit: int | None = None, jobs: int = 4) -> pd.DataFrame:
    frame = select_records(pd.read_csv(manifest, keep_default_na=False), subset, limit)
    training_path = source_root / 'prodigy' / 'elife-07454-supp2.csv'
    training_ids = read_prodigy_training_ids(training_path)
    work_dir = output.parent / 'prodigy_work'
    output.parent.mkdir(parents=True, exist_ok=True)
    arguments = [(record, str(work_dir), training_ids) for record in frame.to_dict('records')]
    rows = []
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        for index, row in enumerate(pool.map(_prodigy_record, arguments)):
            rows.append(row)
            if (index + 1) % 10 == 0 or index == len(arguments) - 1:
                print(f'PRODIGY {index + 1}/{len(arguments)}: {sum(r["status"] == "ok" for r in rows)} successful', flush=True)
    results = pd.DataFrame(rows)
    results.to_csv(output, index=False)
    provenance = {'method': 'PRODIGY', 'target': 'pkd', 'protocol': 'frozen_empirical',
        'training_overlap': 'known' if results.exact_training_pdb_overlap.any() else 'possible',
        'source_url': 'https://github.com/haddocking/prodigy', 'upstream_version': version('prodigy-prot'),
        'license': 'Apache-2.0', 'package': 'prodigy-prot', 'source_manifest': str(manifest.resolve()),
        'subset': subset, 'requested': len(frame), 'successful': int((results.status == 'ok').sum()),
        'temperature_C': 25.0, 'distance_cutoff_angstrom': 5.5, 'accessibility_threshold': 0.05,
        'upstream_native_output': 'thermodynamic binding DeltaG kcal/mol; negative means favorable',
        'output_conversion': 'pKd=-log10(upstream kd_val); upstream R=0.0019858775 kcal/mol/K at298.15K',
        'original_training_dataset': 'Vangone and Bonvin 2015 cleaned affinity benchmark, Supplementary file2',
        'training_id_source': str(training_path), 'training_pdb_ids_available': len(training_ids),
        'exact_training_pdb_overlap_count': int(results.exact_training_pdb_overlap.sum()),
        'notes': 'Native released fixed coefficients and preprocessing. First model and manifest partner chains selected before native parser; upstream insertion-code and nonstandard-residue handling retained. No calibration/retraining. Exact PDB matches flagged; absent exact match does not exclude homologous training systems.',
        'created_utc': datetime.now(timezone.utc).isoformat()}
    output.with_suffix('.json').write_text(json.dumps(provenance, indent=2) + '\n')
    return results


def _load_mirage(source_root: Path):
    import importlib
    import sys
    repository = source_root / 'MIRAGE_author'
    if str(repository) not in sys.path:
        sys.path.insert(0, str(repository))
    return importlib.import_module('run_predict')


def _mirage_region(arguments) -> dict:
    record, work_dir, source_root = arguments
    path = Path(work_dir) / record['pdb_id'] / 'partners.pdb'
    try:
        write_partner_structure(record, path)
        upstream = _load_mirage(Path(source_root))
        region = path.with_suffix('.region')
        if not region.exists():
            upstream.generate_region_file(path, record['chains_a'].split(','), record['chains_b'].split(','), region)
        return {'sample_id': record['sample_id'], 'path': str(path), 'error': ''}
    except Exception as error:
        return {'sample_id': record['sample_id'], 'path': str(path), 'error': f'{type(error).__name__}: {error}'}


def prepare_mirage_regions(manifest: Path, output: Path, source_root: Path,
                           subset: str = 'all', limit: int | None = None,
                           jobs: int = 4) -> list[dict]:
    frame = select_records(pd.read_csv(manifest, keep_default_na=False), subset, limit)
    work_dir = output.parent / 'mirage_work'
    arguments = [(record, str(work_dir), str(source_root)) for record in frame.to_dict('records')]
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        return list(pool.map(_mirage_region, arguments))




def remap_region_positions(region: pd.DataFrame, accepted_old_positions: list[int]) -> pd.DataFrame:
    """Retain native SASA/labels while mapping CA-based positions to atom14 rows."""
    positions = {old: index + 1 for index, old in enumerate(accepted_old_positions)}
    selected = region[region['pos'].isin(positions)].copy()
    selected['pos'] = selected['pos'].map(positions)
    if sorted(selected['pos']) != list(range(1, len(accepted_old_positions) + 1)):
        raise ValueError('Native region rows do not cover the atom14 residue set exactly.')
    return selected


def align_mirage_sequence_region(record: dict, path: Path, original_region: Path,
                                 source_root: Path) -> tuple[dict[str, str], dict]:
    """Match sequence and region row identities to the unchanged author parser."""
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import index_to_one
    upstream = _load_mirage(source_root)
    chains = record['chains_a'].split(',') + record['chains_b'].split(',')
    parsed = upstream.parse_pdb_atom14(path, chains=chains)
    sequences = {}
    for chain, aa in zip(parsed['chain_id'], parsed['aa'], strict=True):
        key = 'PDB|' + chain
        sequences[key] = sequences.get(key, '') + index_to_one(int(aa))
    model = next(PDBParser(QUIET=True).get_structure('alignment', path).get_models())
    old_position = 0
    accepted_old_positions, omitted = [], []
    for chain in model:
        if chain.id not in chains:
            continue
        for residue in chain:
            supported = upstream.is_supported_amino_acid(residue.resname)
            if supported and 'CA' in residue:
                old_position += 1
                if 'N' in residue and 'C' in residue:
                    accepted_old_positions.append(old_position)
                    continue
            omitted.append({'chain': chain.id, 'residue_number': residue.id[1],
                            'insertion_code': residue.id[2].strip(), 'resname': residue.resname,
                            'reason': 'unsupported_residue' if not supported else 'incomplete_N_CA_C'})
    if len(accepted_old_positions) != len(parsed['aa']):
        raise ValueError('Alignment traversal disagrees with original atom14 parser.')
    region = remap_region_positions(pd.read_csv(original_region), accepted_old_positions)
    region.to_csv(path.with_suffix('.region'), index=False)
    return sequences, {'parsed_residues': len(parsed['aa']), 'native_region_residues': old_position,
                       'omitted_residues': omitted, 'sequence_source': 'original parse_pdb_atom14 aa and chain_id arrays',
                       'region_mapping': 'original CA-based positions mapped to original atom14 N/CA/C residue order; SASA/region labels unchanged'}

def prepare_mirage_embeddings(frame: pd.DataFrame, output: Path, source_root: Path, device: str, *, align_to_atom14: bool = False) -> None:
    """Call the original embedder once; only identifiers/device/batch size differ.

    Author weights, tokenizer, float32 precision and residue extraction are
    unchanged. Batch size one keeps the released model within laptop VRAM.
    """
    import contextlib
    import importlib.util
    import os

    import h5py
    import torch

    repository = source_root / 'MIRAGE_author'
    prottrans = repository / 'prottrans'

    def load_script(name, path):
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    fasta_module = load_script('_mirage_pdb_fasta', prottrans / 'pdb_tofasta.py')
    embedder = load_script('_mirage_prott5', prottrans / 'prott5_embedder.py')
    embedder.device = torch.device(device)
    work_name = 'mirage_aligned_work' if align_to_atom14 else 'mirage_work'
    combined = output.parent / work_name / 'combined.fasta'
    combined.parent.mkdir(parents=True, exist_ok=True)
    assignments = []
    fasta_parts = []
    for index, record in enumerate(frame.to_dict('records')):
        work_dir = combined.parent / record['pdb_id']
        path = work_dir / 'partners.pdb'
        write_partner_structure(record, path)
        target = path.with_suffix('.h5')
        if target.exists():
            continue
        fasta_path = path.with_suffix('.fasta')
        if align_to_atom14:
            sequences, audit = align_mirage_sequence_region(record, path, output.parent / 'mirage_work' / record['pdb_id'] / 'partners.region', source_root)
            fasta_path.write_text(''.join(f'>{key}\n{sequence}\n' for key, sequence in sequences.items()))
            (work_dir / 'alignment.json').write_text(json.dumps(audit, indent=2) + '\n')
        else:
            with path.open() as structure:
                fasta_path.write_text(''.join(fasta_module.run(structure, multi=True)))
        sequences = embedder.read_fasta(fasta_path)
        keys = []
        for original_key, sequence in sequences.items():
            combined_key = f'sample{index:04d}|{original_key}'
            fasta_parts.append(f'>{combined_key}\n{sequence}\n')
            keys.append((combined_key, original_key))
        assignments.append((target, keys))
    if not assignments:
        return
    combined.write_text(''.join(fasta_parts))
    combined_h5 = combined.with_suffix('.h5')
    cwd = Path.cwd()
    try:
        os.chdir(prottrans)
        with (combined.parent / 'embedding.log').open('w') as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            embedder.get_embeddings(combined, combined_h5, None, per_protein=False, max_batch=1)
    finally:
        os.chdir(cwd)
    with h5py.File(combined_h5, 'r') as source:
        for target, keys in assignments:
            missing = [key for key, _ in keys if key not in source]
            if missing:
                raise RuntimeError(f'Original embedder did not return requested chains: {missing}')
            with h5py.File(target, 'w') as destination:
                for combined_key, original_key in keys:
                    destination.create_dataset(original_key, data=source[combined_key][:])
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

def run_mirage(manifest: Path, output: Path, source_root: Path, subset: str = 'all',
               limit: int | None = None, device: str = 'cuda', native_results: Path | None = None) -> pd.DataFrame:
    """Use the original author feature pipeline and strict released checkpoint."""
    import contextlib
    import os

    import torch

    frame = select_records(pd.read_csv(manifest, keep_default_na=False), subset, limit)
    output.parent.mkdir(parents=True, exist_ok=True)
    repository = source_root / 'MIRAGE_author'
    checkpoint = repository / 'data' / 'model.pt'
    a1741 = set(pd.read_excel(repository / 'datasets' / 'A1741.xlsx')['PDB'].str.lower())
    ppb = set(pd.read_excel(repository / 'datasets' / 'PPB-Affinity.xlsx')['name'].str.lower())
    upstream = _load_mirage(source_root)
    # Fixed seed belongs to the released checkpoint; CPU thread limits avoid
    # multiplying BLAS pools when its embedding subprocesses start.
    torch.set_num_threads(4)
    os.environ.setdefault('OMP_NUM_THREADS', '4')
    os.environ.setdefault('MKL_NUM_THREADS', '4')
    os.environ.setdefault('TRANSFORMERS_OFFLINE', '1')
    native = pd.read_csv(native_results, keep_default_na=False) if native_results is not None else None
    if native is not None and set(native['sample_id']) != set(frame['sample_id']):
        raise ValueError('The native MIRAGE CSV must cover exactly the selected manifest subset; use the same --subset/--limit.')
    native_successes = {} if native is None else {row['sample_id']: row for row in native.to_dict('records') if row['status'] == 'ok'}
    pending = frame[~frame['sample_id'].isin(native_successes)]
    aligned = native is not None
    method = 'MIRAGE-aligned-inputs' if aligned else 'MIRAGE'
    work_name = 'mirage_aligned_work' if aligned else 'mirage_work'
    prepare_mirage_embeddings(pending, output, source_root, device, align_to_atom14=aligned)
    rows = []
    for index, record in enumerate(frame.to_dict('records')):
        if record['sample_id'] in native_successes:
            result = dict(native_successes[record['sample_id']])
            result['method'] = method
            result['input_alignment'] = 'native_success_retained'
            rows.append(result)
            continue
        started = time.monotonic()
        pdb_id = record['pdb_id'].lower()
        result = {'sample_id': record['sample_id'], 'predicted': np.nan, 'target': 'pkd',
                  'method': method, 'status': 'error', 'error': '',
                  'training_overlap': 'possible', 'in_author_A1741': pdb_id in a1741,
                  'in_author_PPB_Affinity': pdb_id in ppb, 'reference_temperature_K': 298.15}
        result['input_alignment'] = 'sequence_and_region_matched_to_atom14' if aligned else 'native'
        work_dir = output.parent / work_name / pdb_id
        work_dir.mkdir(parents=True, exist_ok=True)
        path = work_dir / 'partners.pdb'
        try:
            write_partner_structure(record, path)
            with (work_dir / 'native_inference.log').open('w') as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                region_path, h5_path = upstream.ensure_feature_files(path, record['chains_a'].split(','), record['chains_b'].split(','))
                score, metadata = upstream.predict_dg(checkpoint, path, region_path, h5_path,
                    chains=record['chains_a'].split(',') + record['chains_b'].split(','), device_arg=device)
            result.update(predicted=native_dg_to_pkd(score, positive_affinity_score=True),
                          native_positive_affinity_kcal_mol=score, native_dg_kcal_mol=-score, status='ok')
            (work_dir / 'native_prediction.json').write_text(json.dumps({'score': score, 'metadata': metadata}, indent=2, default=str) + '\n')
        except Exception as error:
            result['error'] = f'{type(error).__name__}: {error}'
        result['seconds'] = round(time.monotonic() - started, 3)
        rows.append(result)
        pd.DataFrame(rows).to_csv(output, index=False)
        print(f'MIRAGE {index + 1}/{len(frame)} {record["sample_id"]}: {result["status"]} ({result["seconds"]}s) {result["error"]}', flush=True)
    results = pd.DataFrame(rows)
    results.to_csv(output, index=False)
    provenance = {'method': method, 'target': 'pkd', 'protocol': 'native_pretrained',
        'training_overlap': 'possible', 'source_url': 'https://github.com/ShiweiWu-545/MIRAGE',
        'upstream_revision': MIRAGE_REVISION,
        'license': 'No explicit license in author repository snapshot; local research use, upstream material not redistributed in TopoPPI',
        'raw_native_results': str(native_results) if native_results is not None else None,
        'native_coverage_before_alignment': f'{len(native_successes)}/{len(frame)}' if aligned else None,
        'coverage_after_alignment': f'{int((results.status == "ok").sum())}/{len(frame)}' if aligned else None,
        'recovered_input_failures': int((results.status == 'ok').sum()) - len(native_successes) if aligned else None,
        'input_alignment': 'Original atom14 aa/chain arrays define FASTA; original region labels/SASA preserved with exact residue-position remapping. Raw successes retained and only failures rerun.' if aligned else 'unmodified author FASTA/region processing',
        'checkpoint': str(checkpoint), 'checkpoint_loading': 'author predict_dg with strict=True',
        'embedding_execution': 'Original author get_embeddings called once on combined named-chain FASTA, max_batch=1, requested device, original float32 weights and tokenizer. Per-complex HDF5 preserves original chain keys.',
        'runtime_versions': {name: version(name) for name in ['torch', 'transformers', 'biopython', 'h5py', 'sentencepiece']},
        'source_manifest': str(manifest.resolve()), 'subset': subset, 'requested': len(frame),
        'successful': int((results.status == 'ok').sum()), 'device': device,
        'upstream_native_output': 'positive affinity score DG=-RT ln(Kd), kcal/mol; thermodynamic DeltaG=-DG',
        'output_conversion': 'pKd=positive_score/(0.001987*298.15*ln(10)); explicit reference temperature, not temperature extrapolation',
        'reference_temperature_K': 298.15,
        'author_A1741_exact_pdb_overlap': int(results.in_author_A1741.sum()),
        'author_PPB_exact_pdb_overlap': int(results.in_author_PPB_Affinity.sum()),
        'training_split_provenance': 'Default released checkpoint config contains PPB-Affinity path but no unique training split membership. A1741 and PPB candidate-training membership reported per row; contamination cannot be excluded.',
        'notes': 'Native author ShrakeRupley interface selection, original ProtT5 embedder and released graph model. First model/complete manifest partners supplied. No evaluation labels supplied, no model or calibration refit. Potentially training-overlapping inference is not a clean unseen-complex benchmark.',
        'created_utc': datetime.now(timezone.utc).isoformat()}
    output.with_suffix('.json').write_text(json.dumps(provenance, indent=2) + '\n')
    return results
