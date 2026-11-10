"""Observed-sequence correspondence and mapped native contact differences."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.SeqUtils import seq1
from topoppi.atlas.footprints import source_atom_residue_labels
from topoppi.io.io_loader import PDBLoader
from topoppi.sequence_alignment import align_protein_sequences

from topoppi_affinity.features import _chains, _group_data

from .context import atom_context


def observed(row):
    loader = PDBLoader(row["structure_path"])
    groups = sorted((_chains(row["chains_a"]), _chains(row["chains_b"])))
    chains = {}
    for chain in sum((list(g) for g in groups), []):
        xyz, atoms = loader.get_chain_data(chain)
        keys = list(dict.fromkeys(map(str, source_atom_residue_labels(atoms))))
        chains[chain] = ("".join(seq1(k.split(":")[1]) for k in keys), keys)
    residues, edges = atom_context([_group_data(loader, g) for g in groups])
    return groups, chains, residues, edges


def map_structures(left, right):
    ga, ca, _, _ = left
    gb, cb, _, _ = right
    matches = {}
    failures = []
    for a, (sa, _ka) in ca.items():
        for b, (sb, _kb) in cb.items():
            if min(len(sa), len(sb)) / max(len(sa), len(sb)) < 0.8:
                continue
            try:
                pairs, metrics = align_protein_sequences(sa, sb)
            except ValueError as exc:
                failures.append(str(exc))
                continue
            if (
                metrics["alignment_identity"] >= 0.3
                and min(metrics["reference_coverage"], metrics["mobile_coverage"]) >= 0.8
            ):
                matches[a, b] = (pairs, metrics)
    assignments = []
    for orientation in ((0, 1), (1, 0)):
        if any(len(ga[s]) != len(gb[orientation[s]]) for s in (0, 1)):
            continue
        choices = []
        for side in (0, 1):
            choices.append(
                [
                    list(zip(ga[side], perm, strict=True))
                    for perm in itertools.permutations(gb[orientation[side]])
                    if all((a, b) in matches for a, b in zip(ga[side], perm, strict=True))
                ]
            )
        for first, second in itertools.product(*choices):
            assignments.append(first + second)
    if len(assignments) != 1:
        return None, dict(
            status="ambiguous_chain_assignment" if assignments else "no_complete_mapping",
            assignments=len(assignments),
            alignment_failures=len(failures),
        )
    records = []
    for a, b in assignments[0]:
        pairs, metrics = matches[a, b]
        if metrics["selected_pair_consensus_fraction"] != 1:
            return None, dict(status="ambiguous_residue_alignment", chain_a=a, chain_b=b)
        for i, j in pairs:
            records.append(
                dict(
                    reference_key=ca[a][1][i],
                    mobile_key=cb[b][1][j],
                    reference_chain=a,
                    mobile_chain=b,
                    identity=metrics["alignment_identity"],
                    reference_coverage=metrics["reference_coverage"],
                    mobile_coverage=metrics["mobile_coverage"],
                )
            )
    return records, dict(status="automatic_mapping", independent_review="not_performed", assignments=1, mapped_residues=len(records))


def contact_delta(mapping, left_edges, right_edges, cutoff=6.0):
    lookup = {r["reference_key"]: r["mobile_key"] for r in mapping}
    reverse = {v: k for k, v in lookup.items()}
    left = {
        tuple(sorted((e["source"], e["target"]))): e["min_distance"]
        for e in left_edges
        if e["source"] in lookup and e["target"] in lookup
    }
    right = {
        tuple(sorted((reverse[e["source"]], reverse[e["target"]]))): e["min_distance"]
        for e in right_edges
        if e["source"] in reverse and e["target"] in reverse
    }
    rows = []
    for edge in sorted(set(left) | set(right)):
        a, b = left.get(edge, float("inf")), right.get(edge, float("inf"))
        if min(a, b) > cutoff:
            continue
        status = "shared" if max(a, b) <= cutoff else "reference_only" if a <= cutoff else "mobile_only"
        rows.append(
            dict(
                reference_source=edge[0],
                reference_target=edge[1],
                mobile_source=lookup[edge[0]],
                mobile_target=lookup[edge[1]],
                reference_min_distance=a if np.isfinite(a) else None,
                mobile_min_distance=b if np.isfinite(b) else None,
                distance_missing_meaning="above exported 6.5 A distance limit",
                cutoff_A=cutoff,
                status=status,
            )
        )
    return rows


def run(manifest, output, target=24):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(manifest, keep_default_na=False).drop_duplicates("complex_id")
    rng = np.random.default_rng(20260926)
    families = sorted(frame.family_group.unique())
    rng.shuffle(families)
    cache = {}
    registry = []
    pairs = []
    for family in families:
        group = frame[frame.family_group == family].sort_values("sample_id").to_dict("records")
        if len(group) < 2:
            continue
        candidates = list(itertools.combinations(range(len(group)), 2))
        rng.shuffle(candidates)
        for i, j in candidates:
            a, b = group[i], group[j]
            if a["pdb_id"] == b["pdb_id"]:
                continue
            try:
                for row in (a, b):
                    if row["sample_id"] not in cache:
                        cache[row["sample_id"]] = observed(row)
                left, right = cache[a["sample_id"]], cache[b["sample_id"]]
                mapping, report = map_structures(left, right)
            except (ValueError, KeyError, OSError) as exc:
                mapping = None
                report = dict(status="source_error", error=str(exc))
            record = dict(reference=a["sample_id"], mobile=b["sample_id"], family_group=family, **report)
            registry.append(record)
            if mapping is None:
                continue
            pair_id = f"pair_{len(pairs) + 1:03d}"
            folder = output / pair_id
            folder.mkdir(exist_ok=True)
            pd.DataFrame(mapping).to_csv(folder / "correspondence.csv", index=False)
            deltas = [r for cutoff in (5.5, 6.0, 6.5) for r in contact_delta(mapping, left[3], right[3], cutoff)]
            pd.DataFrame(deltas).to_csv(folder / "contact_differences.csv", index=False)
            for label, row, structure in [("reference", a, left), ("mobile", b, right)]:
                (folder / f"{label}.json").write_text(
                    json.dumps(dict(source=row, groups=structure[0], residues=structure[2], contacts=structure[3]))
                    + "\n"
                )
            pairs.append(
                dict(
                    pair_id=pair_id,
                    **record,
                    reference_structure=a["structure_path"],
                    mobile_structure=b["structure_path"],
                )
            )
            print(f"{pair_id}: {a['pdb_id']} / {b['pdb_id']}, {family}", flush=True)
            break  # One pair per candidate family; source structures never reused.
        if len(pairs) >= target:
            break
    pd.DataFrame(registry).to_csv(output / "candidate_registry.csv", index=False)
    pd.DataFrame(pairs).to_csv(output / "pairs.csv", index=False)
    (output / "status.json").write_text(
        json.dumps(
            dict(
                automatically_mapped_pairs=len(pairs),
                independent_candidate_families=len(pairs),
                reviewed_pairs=0,
                human_performance_observations=0,
                selection="seed 20260926; random family/pair order; first eligible pair per family; stop at target",
                target_pairs=target,
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("manifest")
    p.add_argument("output")
    p.add_argument("--target", type=int, default=24)
    a = p.parse_args()
    run(a.manifest, a.output, a.target)
