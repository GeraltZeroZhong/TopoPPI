"""Fit a finite graph readout using caller-supplied training, validation and test groups."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="CSV: path, split; each path locates a graph-probe NPZ and JSON")
    parser.add_argument("--view", default="native", choices=["self", "native", "ordinary", "aware", "aware_uniform"])
    parser.add_argument("--width", type=int, default=32)
    parser.add_argument("--depth", type=int, default=1)
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--patience", type=int, default=20)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    import pandas as pd
    import torch

    from topoppi_representation.analysis.graphs import FEATURE_NAMES, TARGET_NAMES, ReaderConfig, load_probe_sample
    from topoppi_representation.analysis.readers import fit_reader, predict_reader

    manifest = pd.read_csv(args.manifest, keep_default_na=False)
    partitions = {name: [] for name in ("train", "validation", "test")}
    groups = {}
    ids = set()
    for row in manifest.itertuples():
        if row.split not in partitions:
            raise ValueError("split must be train, validation or test")
        path = Path(row.path)
        sample = load_probe_sample(path if path.is_absolute() else args.manifest.parent / path)
        if sample.sample_id in ids:
            raise ValueError("Each sample must occur once")
        ids.add(sample.sample_id)
        if sample.family_id in groups and groups[sample.family_id] != row.split:
            raise ValueError("A relationship group crosses data partitions")
        groups[sample.family_id] = row.split
        partitions[row.split].append(sample)
    if any(not values for values in partitions.values()):
        raise ValueError("Training, validation and test partitions must be nonempty")
    config = ReaderConfig(
        width=args.width,
        depth=args.depth,
        max_epochs=args.epochs,
        patience=args.patience,
        seed=args.seed,
        cpu_threads=args.threads,
        device=args.device,
    )
    model, detail = fit_reader(partitions["train"], partitions["validation"], args.view, config)
    predicted = predict_reader(model, partitions["test"], args.view)
    output = pd.DataFrame(predicted, columns=TARGET_NAMES)
    output.insert(0, "sample_id", [sample.sample_id for sample in partitions["test"]])
    args.output.mkdir(parents=True, exist_ok=False)
    output.to_csv(args.output / "predictions.csv", index=False)
    torch.save(
        {
            "state_dict": model.cpu().state_dict(),
            "config": detail["config"],
            "feature_names": FEATURE_NAMES,
            "target_names": TARGET_NAMES,
            "view": args.view,
        },
        args.output / "reader.pt",
    )
    (args.output / "fit.json").write_text(json.dumps(detail, indent=2) + "\n")


if __name__ == "__main__":
    main()
