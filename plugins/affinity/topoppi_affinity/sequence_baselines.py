"""Geometry-free readouts of frozen observed-chain ESMC representations.

The registry originates from observed structural residues. Reading it does not
make these full experimental-construct sequences: unresolved residues remain
absent. This module opens no coordinate, contact, surface, or footprint file.
Multi-mutants use the jointly mutated chain embedding from the registry.
"""
from __future__ import annotations

import json
import time
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from topoppi_affinity.features import parse_mutations

ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
INPUT_COLUMNS = ("sample_id", "source", "complex_id", "chains_a", "chains_b", "mutations")
SEQUENCE_SCOPE = "observed_structural_sequence; unresolved residues absent; not complete experimental construct"


def _composition(sequence):
    return np.asarray([sequence.count(aa) for aa in ALPHABET], dtype=float) / len(sequence)


def _one_hot(aa):
    result = np.zeros(len(ALPHABET))
    result[ALPHABET.index(aa if aa in ALPHABET else "X")] = 1.
    return result


def _symmetric(left, right):
    return np.concatenate((left + right, np.abs(left - right), left * right))


def sequence_state_features(partners):
    """Pool every residue of each chain, then length-weight chains per partner.

    Each partner is a list of (sequence, embedding) pairs. No interface mask,
    chain selection based on proximity, or truncation is accepted.
    """
    pooled = []
    for chains in partners:
        lengths = np.asarray([len(sequence) for sequence, _ in chains])
        chain_vectors = [np.concatenate((embedding.mean(axis=0), _composition(sequence)))
                         for sequence, embedding in chains]
        mean = np.average(chain_vectors, axis=0, weights=lengths)
        pooled.append(np.concatenate((mean, [np.log1p(lengths.sum()), np.log1p(len(chains))])))
    dimension = partners[0][0][1].shape[1]
    names = [f"embedding.{i}" for i in range(dimension)]
    names += [f"composition.{aa}" for aa in ALPHABET] + ["log_length", "log_chain_count"]
    names = [f"sequence.global.{operation}.{name}" for operation in ("sum", "absdiff", "product")
             for name in names]
    return _symmetric(*pooled).astype(np.float32), np.asarray(names)


def mutation_context_features(wild, mutant, sites, *, context_radius=8):
    """An odd readout of substitution identity and sequence-local state changes.

    ``sites`` contains (chain, zero-based observed-sequence index). Reversing
    wild/mutant negates every coordinate. A local window is defined by sequence
    order, not author-number distance or any spatial neighborhood.
    """
    dimension = next(iter(wild.values()))[1].shape[1]
    site_delta, window_delta = [], []
    identity = np.zeros(len(ALPHABET))
    oriented_pair = np.zeros((len(ALPHABET), len(ALPHABET)))
    identity_context = np.zeros_like(oriented_pair)
    position_identity = np.zeros((2, len(ALPHABET)))
    for chain, position in sites:
        wt_sequence, wt_embedding = wild[chain]
        mt_sequence, mt_embedding = mutant[chain]
        delta = _one_hot(mt_sequence[position]) - _one_hot(wt_sequence[position])
        identity += delta
        oriented_pair += np.outer(_one_hot(wt_sequence[position]), _one_hot(mt_sequence[position]))
        oriented_pair -= np.outer(_one_hot(mt_sequence[position]), _one_hot(wt_sequence[position]))
        left, right = max(0, position - context_radius), min(len(wt_sequence), position + context_radius + 1)
        context = (_composition(wt_sequence[left:right]) + _composition(mt_sequence[left:right])) / 2.
        identity_context += np.outer(delta, context)
        position_identity[0] += delta * position / max(len(wt_sequence) - 1, 1)
        position_identity[1] += delta * min(position, len(wt_sequence) - 1 - position) / len(wt_sequence)
        site_delta.append(mt_embedding[position] - wt_embedding[position])
        window_delta.append(mt_embedding[left:right].mean(axis=0) - wt_embedding[left:right].mean(axis=0))
    # Identity substitutions have no changed sites and must encode exactly zero.
    site_delta = np.asarray(site_delta) if site_delta else np.zeros((1, dimension))
    window_delta = np.asarray(window_delta) if window_delta else np.zeros((1, dimension))
    pair_i, pair_j = np.triu_indices(len(ALPHABET), k=1)
    blocks = [site_delta.sum(axis=0), site_delta.mean(axis=0), window_delta.mean(axis=0), identity,
              oriented_pair[pair_i, pair_j], identity_context.ravel(), position_identity.ravel()]
    names = []
    for scope in ("site_sum", "site_mean", "window_mean"):
        names.extend(f"sequence.context.{scope}.embedding.{i}" for i in range(dimension))
    names.extend(f"sequence.context.identity.{aa}" for aa in ALPHABET)
    names.extend(f"sequence.context.substitution.{ALPHABET[i]}_{ALPHABET[j]}" for i, j in zip(pair_i, pair_j, strict=True))
    names.extend(f"sequence.context.identity_window.{aa}.{bb}" for aa in ALPHABET for bb in ALPHABET)
    names.extend(f"sequence.context.{scope}.{aa}" for scope in ("relative_position", "terminal_distance") for aa in ALPHABET)
    return np.concatenate(blocks).astype(np.float32), np.asarray(names)


class SequenceFeatureStore:
    """Read only existing chain mappings and frozen sequence embeddings."""

    def __init__(self, sequence_dir):
        self.sequence_dir = Path(sequence_dir)
        registry = pd.read_csv(self.sequence_dir / "sequences.csv", keep_default_na=False)
        self.sequences = dict(zip(registry.sequence_id, registry.sequence, strict=True))
        wild = pd.read_csv(self.sequence_dir / "wild_type_chains.csv", keep_default_na=False)
        self.wild = {(str(row.source), str(row.complex_id), str(row.chain)):
                     (str(row.sequence_id), json.loads(row.residue_keys)) for row in wild.itertuples()}
        mutant_path = self.sequence_dir / "mutant_chains.csv"
        mutant = pd.read_csv(mutant_path, keep_default_na=False) if mutant_path.exists() else pd.DataFrame()
        self.mutants = {(str(row.sample_id), str(row.chain)): str(row.sequence_id) for row in mutant.itertuples()}
        self.models = set()
        self._embedding = lru_cache(maxsize=16)(self._load_embedding)

    def _load_embedding(self, sequence_id):
        sequence = self.sequences[sequence_id]
        with np.load(self.sequence_dir / "embeddings" / f"{sequence_id}.npz", allow_pickle=False) as source:
            if str(source["sequence"]) != sequence:
                raise ValueError(f"{sequence_id}: embedding sequence differs from registry")
            embedding = source["embedding"].astype(np.float32)
            self.models.add(str(source["model_name"]))
        if embedding.ndim != 2 or len(embedding) != len(sequence) or not np.isfinite(embedding).all():
            raise ValueError(f"{sequence_id}: invalid residue embedding")
        return sequence, embedding

    def features(self, row, *, task, context_radius=8):
        """Return global and optional mutation-context views for one row."""
        groups = [tuple(sorted(set(str(row[name]).split(",")))) for name in ("chains_a", "chains_b")]
        wild, mutant, positions = {}, {}, {}
        mutations = parse_mutations(row["mutations"]) if task == "mutation" else []
        mutated_chains = {m["chain"] for m in mutations if m["wt"] != m["mt"]}
        for chain in sorted(set(groups[0] + groups[1])):
            sequence_id, residue_keys = self.wild[(str(row["source"]), str(row["complex_id"]), chain)]
            wild[chain] = self._embedding(sequence_id)
            positions[chain] = {str(key).rsplit(":", 1)[-1]: i for i, key in enumerate(residue_keys)}
            if len(residue_keys) != len(wild[chain][0]):
                raise ValueError(f"{row['sample_id']}: residue numbering differs from sequence length")
            # A missing true-mutant mapping is a failure, never a WT fallback.
            mutant_id = self.mutants[(str(row["sample_id"]), chain)] if chain in mutated_chains else sequence_id
            mutant[chain] = self._embedding(mutant_id)
        wild_vector, names = sequence_state_features([[wild[chain] for chain in group] for group in groups])
        if task == "absolute":
            return {"global": (wild_vector, names)}
        sites = []
        expected = {chain: list(values[0]) for chain, values in wild.items()}
        for mutation in mutations:
            chain = mutation["chain"]
            position = positions[chain][f"{mutation['position']}{mutation['insertion']}"]
            if wild[chain][0][position] != mutation["wt"]:
                raise ValueError(f"{row['sample_id']}: WT substitution identity differs from sequence")
            expected[chain][position] = mutation["mt"]
            if mutation["wt"] != mutation["mt"]:
                sites.append((chain, position))
        for chain in wild:
            if mutant[chain][0] != "".join(expected[chain]):
                raise ValueError(f"{row['sample_id']}: joint mutant sequence does not match substitutions")
        mutant_vector, _ = sequence_state_features([[mutant[chain] for chain in group] for group in groups])
        global_delta = mutant_vector - wild_vector
        context, context_names = mutation_context_features(wild, mutant, sorted(set(sites)), context_radius=context_radius)
        return {"global": (global_delta, names),
                "context": (np.concatenate((global_delta, context)), np.concatenate((names, context_names)))}


def prepare_sequence_features(manifest, sequence_dir, output, *, task, limit=None, context_radius=8):
    """Write benchmark-compatible NPZ matrices plus full declared coverage.

    Only INPUT_COLUMNS are read even when the manifest contains affinity labels
    or structure paths. Fitting and choosing hyperparameters are separate.
    """
    if task not in ("absolute", "mutation"):
        raise ValueError("task must be absolute or mutation")
    frame = (manifest.loc[:, INPUT_COLUMNS].copy() if isinstance(manifest, pd.DataFrame) else
             pd.read_csv(manifest, usecols=list(INPUT_COLUMNS), keep_default_na=False))
    if limit is not None:
        frame = frame.iloc[:limit]
    if frame.sample_id.duplicated().any():
        raise ValueError("Sequence feature sample IDs must be unique")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    store = SequenceFeatureStore(sequence_dir)
    matrices, names, ids, coverage = {}, {}, [], []
    started = time.monotonic()
    for index, row in enumerate(frame.to_dict("records")):
        record = {"sample_id": row["sample_id"]}
        try:
            views = store.features(row, task=task, context_radius=context_radius)
            for view, (vector, columns) in views.items():
                if not np.isfinite(vector).all():
                    raise ValueError("Nonfinite sequence readout")
                if view in names and not np.array_equal(names[view], columns):
                    raise ValueError("Sequence embedding dimension changed")
            for view, (vector, columns) in views.items():
                names[view] = columns
                matrices.setdefault(view, []).append(vector)
            ids.append(row["sample_id"])
            record.update(status="ok", error="")
        except (KeyError, ValueError, IndexError, OSError) as error:
            record.update(status="error", error=f"{type(error).__name__}: {error}")
        coverage.append(record)
        if (index + 1) % 100 == 0:
            print(f"Sequence-only {task}: {index + 1}/{len(frame)}; {time.monotonic()-started:.1f}s", flush=True)
    pd.DataFrame(coverage).to_csv(output / "sequence_coverage.csv", index=False)
    target = "pkd" if task == "absolute" else "ddg"
    paths = {}
    for view, vectors in matrices.items():
        path = output / f"{target}_sequence_{view}.npz"
        np.savez_compressed(path, X=np.asarray(vectors), sample_ids=np.asarray(ids, dtype=str), feature_names=names[view])
        paths[view] = str(path)
    report = {
        "task": task, "n_declared": len(frame), "n_success": len(ids), "n_failed": len(frame) - len(ids),
        "input_scope": SEQUENCE_SCOPE, "observed_affinity_labels_loaded": False,
        "coordinate_files_read": False, "contact_masks_used": False, "footprints_used": False,
        "sequence_directory": str(sequence_dir), "embedding_models": sorted(store.models),
        "absolute_readout": "Per-chain all-residue mean, length-weighted partner mean, sum/absolute-difference/product",
        "mutation_readout": "Joint mutant complete-chain state minus WT; optional odd site/window/identity context",
        "context_window_radius_observed_residues": context_radius,
        "state_difference_does_not_imply_experimental_epistasis_validation": True,
        "features": {view: {"path": path, "dimension": len(names[view])} for view, path in paths.items()},
        "elapsed_seconds": time.monotonic() - started,
    }
    (output / "sequence_protocol.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def fit_sequence_regressor(matrix, labels, groups, *, target, family="ridge", alpha=10., gamma_factor=.1):
    """Fit one supplied training fold; the caller owns splits and inner tuning.

    Uses the existing complex weighting, train-only scaling and regression
    implementations. Mutation Ridge/RBF predictions are exactly odd.
    """
    if target not in ("pkd", "ddg") or family not in ("ridge", "rbf"):
        raise ValueError("Expected target pkd/ddg and regression family ridge/rbf")
    gamma = gamma_factor / matrix.shape[1]
    if target == "ddg":
        from topoppi_affinity.benchmark_models import fit_odd_regressor

        return fit_odd_regressor(matrix, labels, groups, family=family, alpha=alpha, gamma=gamma)
    from topoppi_affinity.embedding_baselines import fit_frozen_regressor

    return fit_frozen_regressor(matrix, labels, groups, family=family, params={"alpha": alpha, "gamma": gamma})
