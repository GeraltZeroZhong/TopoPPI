"""Prepare and inspect two explicitly mapped interfaces with linked selection."""

from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from topoppi.interface import prepare_interface
from topoppi.interface_io import save_interface

from topoppi_affinity.attribution_features import buried_sasa_weights
from topoppi_affinity.features import extract_complex_features
from topoppi_affinity.workflow import _generate_layout, _residue_table, load_analysis

from .session import ReviewSession


def prepare_one(args):
    path, binary, predictions = args
    path = Path(path)
    obj = json.loads(path.read_text())
    row = obj["source"]
    out = path.parent / (path.stem + "_analysis")
    out.mkdir(exist_ok=True)
    if (out / "preparation.json").exists():
        return json.loads((out / "preparation.json").read_text())
    start = time.monotonic()
    try:
        prepared = prepare_interface(
            row["structure_path"], row["chains_a"], row["chains_b"], source=row["source"], complex_id=row["complex_id"]
        )
        save_interface(prepared, out / "interface.npz")
        bundle = extract_complex_features(row["structure_path"], row["chains_a"], row["chains_b"], prepared=prepared)
        residues = _residue_table(bundle, prepared)
        residues["buried_sasa_A2"] = buried_sasa_weights(
            row, dict(residue_key=residues.residue_key.to_numpy(dtype=str), x=np.zeros((len(residues), 1)))
        )
        residues.to_csv(out / "residues.csv", index=False)
        value = predictions.get(row["sample_id"])
        pd.DataFrame(
            [
                dict(
                    variant_id="absolute",
                    sample_id=row["sample_id"],
                    predicted=value,
                    units="pKd",
                    target="pkd",
                    sites="[]",
                    mutations="",
                    model="Absolute-affinity family-held-out model",
                    status="ok" if value is not None else "unavailable",
                )
            ]
        ).to_csv(out / "predictions.csv", index=False)
        report = dict(
            format="topoppi-affinity-analysis",
            version=1,
            analysis_id=prepared.analysis_id,
            interface=prepared.metadata(),
            files=dict(interface="interface.npz", predictions="predictions.csv", residues="residues.csv"),
            units="pKd",
            model=dict(
                name="TopoPPI absolute affinity", scope="Retrospective family-held-out absolute predictions", target="pkd"
            ),
            warnings=["Cross-structure residue correspondence is sequence-based; independent manual validation is unavailable."],
            layout={},
            layout_status="pending",
            prediction_status="ok" if value is not None else "unavailable",
        )
        (out / "analysis.json").write_text(json.dumps(report, indent=2) + "\n")
        _generate_layout(prepared, report, out, binary, 600.0, lambda s: None, None, None)
        result = dict(
            pair=path.parent.name,
            side=path.stem,
            status="prepared",
            layout_status=report["layout_status"],
            seconds=time.monotonic() - start,
        )
    except Exception as exc:
        result = dict(
            pair=path.parent.name,
            side=path.stem,
            status="error",
            error=f"{type(exc).__name__}: {exc}",
            seconds=time.monotonic() - start,
        )
    (out / "preparation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(result, flush=True)
    return result


def prepare(root, binary, predictions, workers=3):
    root = Path(root)
    frame = pd.read_csv(predictions)
    frame = frame[(frame.method == "core::topoppi_surface") & (frame.task == "absolute_family")]
    values = dict(zip(frame.sample_id, frame.predicted, strict=True))
    jobs = [
        (str(path), binary, values)
        for path in sorted(root.glob("pair_*/*.json"))
        if path.stem in ("reference", "mobile")
    ]
    with ProcessPoolExecutor(workers) as pool:
        results = list(pool.map(prepare_one, jobs))
    pd.DataFrame(results).to_csv(root / "view_preparation.csv", index=False)


def launch(pair, condition="A", layout="primary"):
    import tkinter as tk
    from tkinter import ttk

    from .gui import ReviewWindow

    pair = Path(pair)
    mapping = pd.read_csv(pair / "correspondence.csv", keep_default_na=False)
    forward = dict(zip(mapping.reference_key, mapping.mobile_key, strict=True))
    reverse = {v: k for k, v in forward.items()}
    root = tk.Tk()
    root.title("Mapped interface comparison")
    ttk.Label(
        root, text="Mapping: automatic sequence correspondence; not independently reviewed. Each atlas retains its own UV coordinates."
    ).pack()
    windows = []
    syncing = [False]

    def selected(event, index):
        if syncing[0] or event["event"] != "selection" or len(windows) != 2:
            return
        lookup = forward if index == 0 else reverse
        target = {lookup[k] for k in event["residue_keys"] if k in lookup}
        other = windows[1 - index]
        if target != other.session.selected:
            syncing[0] = True
            try:
                other.choose(target)
            finally:
                syncing[0] = False

    for index, side in enumerate(("reference", "mobile")):
        session = ReviewSession(load_analysis(pair / (side + "_analysis")), condition=condition, layout=layout)
        windows.append(
            ReviewWindow(root, session, condition_locked=True, event_callback=lambda e, i=index: selected(e, i))
        )
    text = tk.Text(root, width=110, height=20)
    text.pack(fill="both", expand=True)
    delta = pd.read_csv(pair / "contact_differences.csv")
    delta = delta[(delta.cutoff_A == 6.0) & (delta.status != "shared")]
    text.insert("1.0", delta.to_string(index=False))
    text.configure(state="disabled")
    root.mainloop()


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="command", required=True)
    a = sub.add_parser("prepare")
    a.add_argument("root")
    a.add_argument("predictions")
    a.add_argument("--binary", required=True)
    a.add_argument("--workers", type=int, default=3)
    b = sub.add_parser("view")
    b.add_argument("pair")
    b.add_argument("--condition", choices=["D", "G", "A"], default="A")
    b.add_argument("--layout", choices=["primary", "residue_aware"], default="primary")
    args = p.parse_args()
    if args.command == "prepare":
        prepare(args.root, args.binary, args.predictions, args.workers)
    else:
        launch(args.pair, args.condition, args.layout)
