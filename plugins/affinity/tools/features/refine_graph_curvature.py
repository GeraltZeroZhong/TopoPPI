#!/usr/bin/env python3
"""Create a separate graph cache replacing only two native curvature columns."""
import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from topoppi_affinity.analytic_curvature import GRAPH_PROTOCOL, replace_graph_curvature
from topoppi_affinity.graph_features import load_graph, save_graph


def refine(task):
    row, path = task
    started = time.monotonic()
    result = {**row, "old_graph_path": row["graph_path"], "graph_path": str(path)}
    try:
        if path.resolve() == Path(row["graph_path"]).resolve():
            raise ValueError("Analytic curvature output must use a separate cache")
        if path.exists():
            graph = load_graph(path)
        else:
            graph = replace_graph_curvature(load_graph(row["graph_path"]), row["structure_path"])
            save_graph(graph, path)
        result.update(status="ok", curvature_protocol=json.loads(str(graph["metadata_json"]))["protocol"],
                      max_abs_curvature_mean=float(abs(graph["x"][:, 34]).max()),
                      max_curvature_std=float(graph["x"][:, 35].max()), error="")
    except (ValueError, KeyError, OSError) as error:
        result.update(status="error", error=f"{type(error).__name__}: {error}")
    result["seconds"] = time.monotonic() - started
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    frame = pd.read_csv(args.manifest, keep_default_na=False)
    frame = frame[frame.status == "ok"].drop_duplicates(["source", "complex_id"])
    if args.limit:
        frame = frame.iloc[:args.limit]
    args.output.mkdir(parents=True, exist_ok=True)
    tasks = [(row, args.output / row["source"] / (row["complex_id"].replace(":", "_").replace(",", "+") + ".npz"))
             for row in frame.to_dict("records")]
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for future in as_completed([executor.submit(refine, task) for task in tasks]):
            row = future.result()
            rows.append(row)
            pd.DataFrame(rows).to_csv(args.output / "manifest.csv", index=False)
            print(f"Analytic curvature {len(rows)}/{len(tasks)} {row['complex_id']} {row['status']} {row['seconds']:.1f}s", flush=True)
    summary = {"protocol": GRAPH_PROTOCOL, "n_graphs": len(rows), "n_ok": sum(row["status"] == "ok" for row in rows),
               "original_manifest": str(args.manifest), "changed_columns_only": [34, 35], "target_labels_used": False}
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
