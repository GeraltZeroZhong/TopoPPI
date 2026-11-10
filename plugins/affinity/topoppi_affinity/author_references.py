"""Released author references with explicit input and supervision provenance.

ProAffinity-GNN consumes the author's complete, FASTA-aware three-graph inputs.
The convenient coordinate-sequence inference variant is a different protocol and
is deliberately not relabelled as the original method here. Graph preparation is
the upstream ``pdb2graph*.py`` / ``graph_construct*.py`` pipeline.
"""

from __future__ import annotations

import ast
import json
import math
import pickle
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def run_bindpred_features(manifest, feature_npz, checkpoint, output):
    """Released general BindPred head on explicitly supplied author-format features.

    The public inference notebook takes a 1280-dimensional ligand vector followed
    by a 1280-dimensional receptor vector. It does not define extraction from a
    multi-chain complex. This adapter does not quietly replace ESM2 with ESMC or
    assume an undocumented chain-pooling convention.
    """
    import numpy as np
    from catboost import CatBoostRegressor

    from topoppi_affinity.benchmark_models import LoadedFeatures

    frame = pd.read_csv(manifest, usecols=['sample_id'], keep_default_na=False)
    features = LoadedFeatures.read(feature_npz)
    if features.X.shape[1] != 2560:
        raise ValueError('BindPred expects ligand 1280 + receptor 1280 ESM2 features')
    matrix = features.aligned(frame.sample_id)
    good = np.isfinite(matrix).all(axis=1)
    model = CatBoostRegressor(thread_count=1)
    model.load_model(str(checkpoint))
    result = frame.copy()
    result['predicted'] = np.nan
    started = time.monotonic()
    if good.any():
        # The author's train.py fits log10(Kd), whereas the benchmark uses pKd.
        result.loc[good, 'predicted'] = -model.predict(matrix[good])
    result['target'], result['method'] = 'pkd', 'BindPred-author-features'
    result['status'] = np.where(good & np.isfinite(result.predicted), 'ok', 'error')
    result['error'] = np.where(result.status.eq('ok'), '', 'Missing/nonfinite supplied feature or prediction')
    result['comparison_scope'] = 'author_supervised_reference'
    result['training_overlap'] = 'possible_PPB_supervised_exposure'
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)
    (output.with_suffix('.json')).write_text(json.dumps({
        'method': 'BindPred', 'checkpoint': str(checkpoint), 'feature_source': str(feature_npz),
        'source_url': 'https://huggingface.co/hbp5181/BindPred',
        'paper_url': 'https://doi.org/10.1093/bioinformatics/btag309',
        'native_target': 'log10(Kd_M)', 'conversion_to_pkd': 'negative sign',
        'input_protocol': 'Author-supplied ligand 1280 followed by receptor 1280, unchanged',
        'sequence_to_feature_preprocessing_reproduced': False,
        'declared': len(result), 'success': int(result.status.eq('ok').sum()),
        'seconds': time.monotonic()-started,
    }, indent=2)+'\n')
    return result


def read_author_partners(path: Path) -> dict[str, tuple[frozenset, frozenset]]:
    """Read the ProAffinity data index, whose chain groups are semicolon separated."""
    result = {}
    for line in path.read_text().splitlines():
        pdb, groups = line.split("\t", 1)
        first, second = [part for part in groups.split(";") if part.strip()]
        result[pdb[:4].lower()+pdb[4:]] = tuple(
            frozenset(chain.strip() for chain in group.split(","))
            for group in (first, second)
        )
    return result


def match_author_partners(record: dict, partners: tuple[frozenset, frozenset]) -> None:
    """A PDB identifier alone does not identify the interaction measured in a row."""
    declared = tuple(frozenset(record[key].split(",")) for key in ("chains_a", "chains_b"))
    if declared != partners and declared != partners[::-1]:
        raise ValueError("Manifest partners differ from the author graph partners")


def _load_proaffinity_model(repository: Path, graph, device: str):
    """Import the two original model classes without executing upstream test data I/O."""
    import torch
    import torch.nn.functional as functional
    from torch_geometric.nn.models import AttentiveFP

    source = repository / "test.py"
    parsed = ast.parse(source.read_text())
    definitions = [node for node in parsed.body
                   if isinstance(node, ast.ClassDef)
                   and node.name in {"AttentiveFPModel", "GraphNetwork"}]
    namespace = {"torch": torch, "F": functional, "AttentiveFP": AttentiveFP}
    exec(compile(ast.Module(body=definitions, type_ignores=[]), str(source), "exec"), namespace)
    model = namespace["GraphNetwork"](
        graph.num_node_features, 256, 64, graph.num_edge_features, 3, 2, 0.5, 32, 1,
    ).to(device)
    state = torch.load(repository / "model/model.pkl", map_location=device, weights_only=True)
    # The author's *_compatible.py makes this same PyG naming update. Both old
    # names reference one shared projection; retain strict checkpoint loading.
    renamed = {}
    for key, value in state.items():
        new_key = key.replace(".lin_src.", ".lin.").replace(".lin_dst.", ".lin.")
        if new_key in renamed and not torch.equal(renamed[new_key], value):
            raise ValueError("ProAffinity checkpoint has distinct source/destination projections")
        renamed[new_key] = value
    model.load_state_dict(renamed)
    return model.eval()


def run_proaffinity_graphs(
    manifest: Path,
    output: Path,
    repository: Path,
    *,
    graph_root: Path | None = None,
    partner_index: Path | None = None,
    device: str = "cuda",
    limit: int | None = None,
) -> pd.DataFrame:
    """Run the released full model on original-format three-graph inputs.

    ``author_graph_key`` may specify the graph basename; otherwise ``pdb_id`` is
    used. The corresponding author-format partner index is mandatory. The four
    bundled author examples are useful for smoke tests but are not an unseen
    test set. User-generated graphs retain their own partner index.
    """
    import torch
    from torch_geometric.data import Batch

    frame = pd.read_csv(manifest, keep_default_na=False)
    if limit is not None:
        frame = frame.iloc[:limit]
    if frame.sample_id.duplicated().any():
        raise ValueError("Reference manifest sample_id values must be unique")
    graph_root = graph_root or repository / "data/graph"
    partner_index = partner_index or repository / "data/index_example.txt"
    partners = read_author_partners(partner_index)
    output.parent.mkdir(parents=True, exist_ok=True)
    rows, model = [], None
    started = time.monotonic()
    for record in frame.to_dict("records"):
        row = {
            "sample_id": record["sample_id"], "predicted": float("nan"), "target": "pkd",
            "method": "ProAffinity-GNN", "status": "error", "error": "",
            "comparison_scope": "author_supervised_reference",
            "training_overlap": "possible", "exposure": "author_supervised_checkpoint_possible_overlap",
            "input_protocol": "author_full_fasta_three_graphs",
            "input_origin": "author_bundled_fixture" if graph_root == repository / "data/graph"
                            else "upstream_preprocessed_graphs",
        }
        key = record.get("author_graph_key") or record["pdb_id"].lower()
        row["author_graph_key"] = key
        try:
            match_author_partners(record, partners[key])
            paths = [graph_root / "inter_graph" / key,
                     graph_root / "individual_graph" / f"{key}_1",
                     graph_root / "individual_graph" / f"{key}_2"]
            graphs = []
            for path in paths:
                with path.open("rb") as stream:
                    graph = pickle.load(stream)
                graph.edge_attr = graph.edge_attr.float()
                graphs.append(graph)
            if model is None:
                model = _load_proaffinity_model(repository, graphs[0], device)
            batch = [Batch.from_data_list([graph]).to(device) for graph in graphs]
            with torch.inference_mode():
                prediction = float(model(*batch).reshape(-1)[0].cpu())
            if not math.isfinite(prediction):
                raise ValueError("Author model produced a nonfinite prediction")
            row.update(predicted=prediction, status="ok",
                       author_partners=";".join(",".join(sorted(group)) for group in partners[key]))
        except Exception as error:
            row["error"] = f"{type(error).__name__}: {error}"
        rows.append(row)
    results = pd.DataFrame(rows) if rows else pd.DataFrame(columns=[
        "sample_id", "predicted", "target", "method", "status", "error", "training_overlap", "exposure",
    ])
    results.to_csv(output, index=False)
    metadata = {
        "method": "ProAffinity-GNN", "target": "pkd", "requested": len(frame),
        "successful": int((results.status == "ok").sum()), "failed": int((results.status != "ok").sum()),
        "source_manifest": str(manifest.resolve()), "repository": str(repository.resolve()),
        "source_url": "https://github.com/legendzzy/ProAffinity-GNN",
        "paper_url": "https://doi.org/10.1021/acs.jcim.4c01850",
        "checkpoint": str(repository / "model/model.pkl"),
        "graph_root": str(graph_root), "partner_index": str(partner_index),
        "protocol": "Original complete FASTA-aware three-graph checkpoint; original test.py classes",
        "runtime_adaptation": "Author-compatible PyG lin_src/lin_dst to lin rename; shared tensors verified equal; strict load",
        "native_output": "pKa affinity score in author code; Kd-labelled source records interpret it as pKd",
        "output_conversion": "identity; no fitting or calibration",
        "training_overlap": "possible",
        "evaluation_role": "author-reference only until per-row supervision and homology exposure are resolved",
        "input_limit": "Requires author-preprocessed graphs; this runner does not substitute ESMC or repacking tools",
        "seconds": time.monotonic() - started, "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    output.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    return results
