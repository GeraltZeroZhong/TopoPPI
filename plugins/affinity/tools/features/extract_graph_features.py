#!/usr/bin/env python3
"""Extract label-free native residue and surface graphs from declared manifests."""
from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from topoppi_affinity.graph_features import GRAPH_PROTOCOL, extract_complex_graph, save_graph


def _extract(task):
    row, output, feature_cache, max_points = task
    start = time.monotonic()
    result = {"complex_id": row["complex_id"], "source": row["source"], "structure_path": row["structure_path"],
              "chains_a": row["chains_a"], "chains_b": row["chains_b"], "graph_path": str(output)}
    try:
        if output.exists():
            with np.load(output, allow_pickle=False) as graph:
                metadata = json.loads(str(graph["metadata_json"]))
                if metadata["protocol"] == GRAPH_PROTOCOL and metadata["max_surface_points_per_side"] == max_points:
                    return {**result, "status": "ok", "cached": True, "n_nodes": len(graph["x"]),
                            "n_edges": graph["edge_index"].shape[1], "n_points": len(graph["point_pos"]),
                            "seconds": time.monotonic() - start}
        graph = extract_complex_graph(row["structure_path"], row["chains_a"], row["chains_b"],
                                      feature_cache, max_points)
        save_graph(graph, output)
        return {**result, "status": "ok", "cached": False, "n_nodes": len(graph["x"]),
                "n_edges": graph["edge_index"].shape[1], "n_points": len(graph["point_pos"]),
                "seconds": time.monotonic() - start}
    except Exception as error:
        return {**result, "status": "error", "error": f"{type(error).__name__}: {error}",
                "seconds": time.monotonic() - start}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifests", nargs="+", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--feature-cache", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--max-surface-points", type=int, default=2048)
    parser.add_argument("--manifest-name", default="manifest.csv", help="Output index filename for the supplied manifests")
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    frame = pd.concat([pd.read_csv(path, keep_default_na=False) for path in args.manifests], ignore_index=True)
    frame = frame.drop_duplicates(["source", "complex_id", "structure_path"])
    if args.limit:
        frame = frame.head(args.limit)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tasks = [(row, args.output_dir / row["source"] / (row["complex_id"].replace(":", "_").replace(",", "+") + ".npz"),
              args.feature_cache, args.max_surface_points) for row in frame.to_dict("records")]
    results = []
    started = time.monotonic()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(_extract, task) for task in tasks]
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            pd.DataFrame(results).to_csv(args.output_dir / args.manifest_name, index=False)
            print(f"{len(results)}/{len(tasks)} {result['complex_id']} {result['status']} "
                  f"{result['seconds']:.1f}s {result.get('error', '')}", flush=True)
    summary = {"protocol": GRAPH_PROTOCOL, "n_complexes": len(results),
               "n_ok": sum(row["status"] == "ok" for row in results), "seconds": time.monotonic() - started,
               "inputs": [str(path) for path in args.manifests], "max_surface_points_per_side": args.max_surface_points,
               "supervised_labels_in_graphs": False}
    summary_name = "summary.json" if args.manifest_name == "manifest.csv" else Path(args.manifest_name).stem + "_summary.json"
    (args.output_dir / summary_name).write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
