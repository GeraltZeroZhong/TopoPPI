"""Native residue context shared by all review conditions."""

from __future__ import annotations

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from topoppi.atlas.footprints import source_atom_residue_labels
from topoppi.io.io_loader import PDBLoader

from topoppi_affinity.features import _chains, _group_data


def atom_context(groups, cutoff=6.5):
    coords = [g[0] for g in groups]
    labels = [source_atom_residue_labels(g[1]).astype(str) for g in groups]
    distances, counts = {}, {}
    for i, indices in enumerate(cKDTree(coords[1]).query_ball_point(coords[0], cutoff)):
        for j in indices:
            key = (str(labels[0][i]), str(labels[1][j]))
            d = float(np.linalg.norm(coords[0][i] - coords[1][j]))
            distances[key] = min(d, distances.get(key, float("inf")))
            if d <= 6.0:
                counts[key] = counts.get(key, 0) + 1
    residues = []
    for side in (0, 1):
        for key in sorted(set(labels[side])):
            xyz = coords[side][labels[side] == key].mean(axis=0)
            residues.append(dict(residue_key=key, side=side, x=float(xyz[0]), y=float(xyz[1]), z=float(xyz[2])))
    edges = [
        dict(source=a, target=b, min_distance=d, atom_pair_count_6A=counts.get((a, b), 0))
        for (a, b), d in sorted(distances.items())
    ]
    return residues, edges


def export_job(args):
    row, output = args
    loader = PDBLoader(row["structure_path"])
    groups = sorted((_chains(row["chains_a"]), _chains(row["chains_b"])))
    residue, edges = atom_context([_group_data(loader, g) for g in groups])
    obj = dict(
        sample_id=row["sample_id"],
        structure_path=row["structure_path"],
        groups=groups,
        family_group=row.get("family_group"),
        residues=residue,
        contacts=edges,
        contact_cutoff_A=6.0,
        exported_distance_limit_A=6.5,
    )
    path = Path(output) / (row["sample_id"].replace(":", "_") + ".json")
    path.write_text(json.dumps(obj) + "\n")
    contact = [e for e in edges if e["min_distance"] <= 6.0]
    return dict(
        sample_id=row["sample_id"],
        contact_count=len(contact),
        chain_count=sum(map(len, groups)),
        contacting_residues=len({e[k] for e in contact for k in ("source", "target")}),
        context=str(path),
    )


def run(manifest, output, workers=4):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(manifest, keep_default_na=False).drop_duplicates("sample_id")
    with ProcessPoolExecutor(workers) as pool:
        rows = list(pool.map(export_job, [(r, str(output)) for r in frame.to_dict("records")]))
    pd.DataFrame(rows).to_csv(output / "native_covariates.csv", index=False)
    print(f"Exported native context for {len(rows)} structures")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("manifest")
    p.add_argument("output")
    p.add_argument("--workers", type=int, default=4)
    a = p.parse_args()
    run(a.manifest, a.output, a.workers)


def independent_contact_sets(structure_path, groups, cutoff=6.0):
    """Answer-key reference: original atoms, explicit identities, dense distances.

    This path does not read the saved contact map or the display's graph.
    """
    from scipy.spatial.distance import cdist

    loader = PDBLoader(structure_path)
    coordinates = []
    identities = []
    for group in groups:
        xyz, atoms = _group_data(loader, group)
        coordinates.append(xyz)
        identities.append(
            np.array(
                [
                    f"{a.get_parent().get_parent().id}:{a.get_parent().resname.strip()}:{a.get_parent().id[1]}{a.get_parent().id[2].strip()}"
                    for a in atoms
                ]
            )
        )
    result = {key: set() for labels in identities for key in set(labels)}
    # Blocked dense evaluation bounds memory while retaining a separate distance path.
    for start in range(0, len(coordinates[0]), 256):
        a, b = np.nonzero(cdist(coordinates[0][start : start + 256], coordinates[1]) <= cutoff)
        for i, j in zip(a + start, b, strict=True):
            left, right = str(identities[0][i]), str(identities[1][j])
            result[left].add(right)
            result[right].add(left)
    return result
