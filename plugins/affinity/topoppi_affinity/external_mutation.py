"""Adapters for original released RDE-Net and StaB-ddG models.

Upstream neural architectures, weights and prediction rules are imported from
their downloaded source trees. These adapters only align inputs and outputs.
"""

from __future__ import annotations

import copy
import json
import os
import pickle
import random
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from topoppi_affinity.data import parse_mutations
from topoppi_affinity.external_inputs import external_partner_groups, prepare_numbered_complex, remap_mutations
from topoppi_affinity.runtime import check_cancel


def read_mutation_inputs(manifest, limit_complexes=None):
    """Read exactly the supplied inputs, optionally limiting complete complexes."""
    frame = pd.read_csv(manifest, keep_default_na=False)
    if limit_complexes is not None:
        frame = frame[frame.complex_id.isin(sorted(frame.complex_id.unique())[:limit_complexes])]
    return frame.copy()


def _save_predictions(records, output_dir, stem, provenance):
    result = pd.DataFrame(records).sort_values("sample_id")
    result.to_csv(output_dir / f"{stem}.csv", index=False)
    overlap_counts = result.training_overlap.value_counts().to_dict()
    # Known overlap in any row takes precedence; unresolved rows remain visible.
    overlap = next(state for state in ("known", "unknown", "possible", "excluded") if state in overlap_counts)
    provenance.update({"n_attempted": len(result), "n_success": int((result.status == "ok").sum()),
                       "training_overlap": overlap, "training_overlap_counts": overlap_counts})
    (output_dir / f"{stem}.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return result


def run_rde(manifest, repository, output_dir, *, device="cpu", batch_size=8, limit_complexes=None, records_frame=None):
    import torch
    from Bio.PDB import PDBParser

    repository, output_dir = Path(repository).resolve(), Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(repository))
    # The author's public checkpoints serialize their EasyDict configuration.
    os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
    from rde.models.rde_ddg import DDG_RDE_Network
    from rde.utils.data import PaddingCollate
    from rde.utils.protein.constants import AA
    from rde.utils.protein.parsers import parse_biopython_structure
    from rde.utils.transforms import get_transform

    torch.set_num_threads(4)
    torch.manual_seed(42)
    checkpoint_paths = list((repository / "trained_models").rglob("DDG_RDE_Network_30k.pt"))
    if len(checkpoint_paths) != 1:
        raise ValueError("Expected one original DDG_RDE_Network_30k.pt checkpoint")
    checkpoint_path = checkpoint_paths[0]
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    config = checkpoint["config"]
    config.model.rde_checkpoint = str(checkpoint_path.with_name("RDE.pt"))
    models = []
    for state in checkpoint["model"]["models"]:
        model = DDG_RDE_Network(config.model)
        model.load_state_dict(state, strict=True)
        models.append(model.to(device).eval())
    transform = get_transform(config.data.transform)
    collate = PaddingCollate()
    frame = read_mutation_inputs(manifest, limit_complexes) if records_frame is None else records_frame.copy()
    records, structures = [], {}
    batch, rows = [], []
    start = time.monotonic()

    def predict_batch():
        if not batch:
            return
        tensors = collate(batch)
        tensors = {key: value.to(device) if isinstance(value, torch.Tensor) else value for key, value in tensors.items()}
        with torch.inference_mode():
            predictions = torch.stack([model(tensors)[1]["ddG_pred"] for model in models]).mean(0).cpu().numpy()
        if not np.isfinite(predictions).all():
            raise ValueError("Original RDE model emitted a nonfinite prediction")
        records.extend({"sample_id": row.sample_id, "predicted": float(value), "target": "ddg", "method": "RDE-Net",
                        "status": "ok", "error": "", "training_overlap": "possible"}
                       for row, value in zip(rows, predictions, strict=True))
        print(f"RDE {len(records)}/{len(frame)} observations; {time.monotonic() - start:.1f}s", flush=True)
        batch.clear()
        rows.clear()

    for _, row in frame.iterrows():
        try:
            if row.complex_id not in structures:
                original = PDBParser(QUIET=True).get_structure("rde", row.structure_path)[0]
                selected = [original[c] for c in sorted(set(row.chains_a.split(",") + row.chains_b.split(",")))]
                structures[row.complex_id] = parse_biopython_structure(selected)
            data, sequence_map = copy.deepcopy(structures[row.complex_id])
            data["aa_mut"] = data["aa"].clone()
            mutations = parse_mutations(row.mutations)
            for wt, chain, number, insertion, mt in mutations:
                key = (chain, number, insertion or " ")
                if key not in sequence_map:
                    raise ValueError(f"RDE parser omitted mutation residue {key}")
                index = sequence_map[key]
                if int(data["aa"][index]) != int(AA(wt)):
                    raise ValueError(f"RDE WT identity mismatch at {key}")
                data["aa_mut"][index] = int(AA(mt))
            data["mut_flag"] = data["aa"] != data["aa_mut"]
            data["ddG"] = 0.0  # Original forward returns a loss; labels never affect its prediction.
            data = transform(data)
            if int(data["mut_flag"].sum()) != len(mutations):
                raise ValueError("RDE 128-residue patch does not retain every mutation")
            batch.append(data)
            rows.append(row)
            if len(batch) >= batch_size:
                predict_batch()
        except (ValueError, KeyError, IndexError) as exc:
            records.append({"sample_id": row.sample_id, "predicted": np.nan, "target": "ddg", "method": "RDE-Net",
                            "status": "error", "error": str(exc), "training_overlap": "possible"})
    predict_batch()
    provenance = {"method": "RDE-Net", "target": "ddg", "protocol": "native_pretrained",
        "source_url": "https://github.com/luost26/RDE-PPI", "license": "Apache-2.0", "checkpoint": str(checkpoint_path),
        "checkpoint_iteration": int(checkpoint["iteration"]), "fold_models": len(models), "ensemble": "mean of released three models, as author pointmut inference",
        "preprocessing": "Original parser and checkpoint transforms; complete partner chains; original residue numbering; 128-residue mutation-centered patch",
        "sign": "Native Gmut-Gwt, positive weaker; no sign reversal", "units": "kcal/mol", "device": device,
        "elapsed_seconds": time.monotonic() - start,
        "notes": "The three released cross-validation models were fitted on SKEMPI. Training overlap with supplied structures is possible and has not been checked. Inference uses structure and mutation inputs."}
    return _save_predictions(records, output_dir, "rde", provenance)


def run_stabddg(manifest, repository, output_dir, *, checkpoint_name="stabddg.pt", device="cpu",
                mc_samples=20, batch_tokens=2000, limit_complexes=None, records_frame=None,
                cancel_event=None, log=None):
    import torch

    log = log or (lambda text: print(text, file=sys.stderr, flush=True))
    check_cancel(cancel_event)
    repository, output_dir = Path(repository).resolve(), Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(repository))
    from stabddg.model import StaBddG
    from stabddg.mpnn_utils import ProteinMPNN, parse_PDB
    from stabddg.ppi_dataset import SKEMPIDataset

    if checkpoint_name not in {"stabddg.pt", "stability_finetuned.pt"}:
        raise ValueError("Use the released binding or stability-finetuned checkpoint")
    stability_only = checkpoint_name == "stability_finetuned.pt"
    method = "StaB-ddG-stability-only" if stability_only else "StaB-ddG"
    stem = "stabddg_stability" if stability_only else "stabddg"
    frame = read_mutation_inputs(manifest, limit_complexes) if records_frame is None else records_frame.copy()
    pdb_dir = output_dir / "pdbs"
    pdb_dir.mkdir(exist_ok=True)
    mapping_rows, failures = [], []
    for complex_id, group in frame.groupby("complex_id", sort=True):
        check_cancel(cancel_event)
        row = group.iloc[0]
        name = complex_id.replace(":", "-").replace(",", "")
        try:
            mapping = prepare_numbered_complex(row, pdb_dir / f"{name}.pdb")
            binder_a, binder_b = external_partner_groups(row.chains_a, row.chains_b)
            for _, observation in group.iterrows():
                mapping_rows.append({"sample_id": observation.sample_id, "#Pdb": f"{name}_{binder_a}_{binder_b}",
                    "Mutation(s)_cleaned": remap_mutations(observation.mutations, mapping), "ddG": 0., "pdb_id": row.pdb_id})
        except (ValueError, KeyError) as exc:
            failures.extend({"sample_id": sample, "predicted": np.nan, "target": "ddg", "method": method,
                             "status": "error", "error": str(exc), "training_overlap": "unknown"} for sample in group.sample_id)
    mapped = pd.DataFrame(mapping_rows)
    mapped.to_csv(output_dir / "input_mapping.csv", index=False)
    # Independent experimental repetitions share a structural prediction.
    unique = mapped.drop_duplicates(["#Pdb", "Mutation(s)_cleaned"])
    unique.to_csv(output_dir / "input.csv", index=False)
    class InputDataset(SKEMPIDataset):
        def preprocess_structures(self, pdb_dir="", pdb_dict_cache_path=""):
            # Original parse_PDB supports an explicit chain list. Restrict its
            # 352 candidate-chain file scans to chains actually present, retaining
            # the author's chain ordering and identical parsed numerical data.
            import string

            chain_order = list(string.ascii_uppercase + string.ascii_lowercase) + [str(i) for i in range(300)]
            structures = []
            for name in self.pdb_file_names:
                check_cancel(cancel_event)
                path = Path(pdb_dir) / f"{name}.pdb"
                present = {line[21] for line in path.read_text().splitlines() if line.startswith("ATOM  ")}
                parsed = parse_PDB(str(path), input_chain_list=[chain for chain in chain_order if chain in present])[0]
                parsed["masked_list"] = [key.split("_")[-1] for key in parsed if key.startswith("seq_chain_")]
                parsed["visible_list"] = []
                structures.append(parsed)
            return structures

    dataset = InputDataset(csv_path=str(output_dir / "input.csv"), pdb_dir=str(pdb_dir), pdb_dict_cache_path="")
    torch.set_num_threads(4)
    torch.manual_seed(42)
    np.random.seed(42)
    random.seed(42)
    pmpnn = ProteinMPNN(node_features=128, edge_features=128, hidden_dim=128, num_encoder_layers=3,
                       num_decoder_layers=3, k_neighbors=48, dropout=0., augment_eps=0.)
    checkpoint_path = repository / "model_ckpts" / checkpoint_name
    state = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    pmpnn.load_state_dict(state.get("model_state_dict", state), strict=True)
    model = StaBddG(pmpnn=pmpnn, use_antithetic_variates=True, noise_level=.1, device=device).to(device).eval()
    with (repository / "data/SKEMPI/train_pdb.pkl").open("rb") as handle:
        author_train = {name.split("_")[0].lower() for name in pickle.load(handle)}
    records = list(failures)
    start = time.monotonic()
    for number, sample in enumerate(dataset, start=1):
        check_cancel(cancel_event)
        complex_sequences = sample["complex_mut_seqs"].to(device)
        a_sequences = sample["binder1_mut_seqs"].to(device)
        b_sequences = sample["binder2_mut_seqs"].to(device)
        prediction_samples = []
        size = max(1, batch_tokens // complex_sequences.shape[1])
        with torch.inference_mode():
            for draw in range(mc_samples):
                check_cancel(cancel_event)
                log(f"{method}: complex {number}/{len(dataset)}, draw {draw + 1}/{mc_samples}")
                values = []
                for offset in range(0, len(complex_sequences), size):
                    check_cancel(cancel_event)
                    part = slice(offset, offset + size)
                    values.append(model(sample["complex"], sample["binder1"], sample["binder2"],
                        complex_sequences[part], a_sequences[part], b_sequences[part]).detach().cpu())
                prediction_samples.append(torch.cat(values))
        # The author CLI reverses the internal folding-score sign here.
        predicted = -torch.stack(prediction_samples).mean(0).numpy()
        if not np.isfinite(predicted).all():
            raise ValueError(f"StaB-ddG produced nonfinite values for {sample['name']}")
        for mutations, value in zip(sample["mutation_list"], predicted, strict=True):
            matched = mapped[(mapped["#Pdb"] == sample["name"]) & (mapped["Mutation(s)_cleaned"] == ",".join(mutations))]
            for _, row in matched.iterrows():
                overlap = "excluded" if stability_only else "known" if row.pdb_id.lower() in author_train else "possible"
                records.append({"sample_id": row.sample_id, "predicted": float(value), "target": "ddg", "method": method,
                                "status": "ok", "error": "", "training_overlap": overlap,
                                "author_training_pdb": row.pdb_id.lower() in author_train})
        log(f"{method} {number}/{len(dataset)} complexes; {len(records)} observations; {time.monotonic()-start:.1f}s")
        pd.DataFrame(records).to_csv(output_dir / f"{stem}_progress.csv", index=False)
    provenance = {"method": method, "target": "ddg", "protocol": "native_pretrained",
        "source_url": "https://github.com/LDeng0205/StaB-ddG", "license": "MIT", "checkpoint": str(checkpoint_path),
        "training_data": "ProteinMPNN structures and Megascale folding stability" + ("" if stability_only else "; SKEMPI binding fine-tuning"),
        "overlap_definition": "Stability-only excludes supervised SKEMPI binding labels, not structural pretraining homology; binding-supervised model flags exact author training PDBs and leaves other homologous overlap possible",
        "input_scope": "supplied manifest rows" if records_frame is None else "explicit caller-provided records",
        "mc_samples": mc_samples, "backbone_noise_A": .1, "antithetic_variates": True, "batch_tokens": batch_tokens,
        "seed": 42, "device": device, "units": "kcal/mol", "sign": "Negative internal binding output, as author CLI; positive weaker",
        "preprocessing": "Complete declared partner groups; observed residues renumbered per chain with WT checks and saved mapping; parser restricted to observed chains; original model and mutation dataset logic",
        "elapsed_seconds": time.monotonic() - start,
        "notes": "Stability-supervised author checkpoint; binding-supervised weights have a separate training scope." if stability_only else
                 "The released binding-finetuned model was trained on SKEMPI. Per-row flags identify exact author training PDBs; overlap counts summarize the supplied inputs."}
    return _save_predictions(records, output_dir, stem, provenance)
