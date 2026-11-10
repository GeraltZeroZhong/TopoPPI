"""Fit a physical-target reader from named features and explicit data partitions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("features", type=Path, help="NPZ with X, feature_names and sample_ids")
    parser.add_argument("targets", type=Path, help="CSV with sample_id, group, split and 13 physical target columns")
    parser.add_argument("--kind", choices=["ridge", "rbf"], default="ridge")
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--gamma", type=float, default=0.01)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    import joblib
    import numpy as np
    import pandas as pd

    from topoppi_representation.analysis.graphs import TARGET_NAMES
    from topoppi_representation.analysis.physical_recovery import fit_physical_reader, validation_loss

    table = pd.read_csv(args.targets, dtype={"sample_id": str, "group": str})
    if table.sample_id.duplicated().any():
        raise ValueError("Each sample must have one target row")
    if not table.split.isin(["train", "validation", "test"]).all():
        raise ValueError("split must be train, validation or test")
    if table.groupby("group").split.nunique().gt(1).any():
        raise ValueError("A relationship group crosses data partitions")
    with np.load(args.features, allow_pickle=False) as archive:
        ids, features, names = archive["sample_ids"].astype(str), archive["X"], archive["feature_names"].astype(str)
    if len(set(ids)) != len(ids) or set(ids) != set(table.sample_id):
        raise ValueError("Feature and target sample identities must match exactly")
    table = table.set_index("sample_id").loc[ids].reset_index()
    targets = table[list(TARGET_NAMES)].to_numpy(float)
    mask = np.isfinite(targets)
    samples = [
        SimpleNamespace(sample_id=row.sample_id, family_id=row.group, target=y, target_mask=m)
        for row, y, m in zip(table.itertuples(), targets, mask, strict=False)
    ]
    training = np.flatnonzero(table.split.eq("train"))
    testing = np.flatnonzero(table.split.eq("test"))
    if not len(training) or not len(testing):
        raise ValueError("Provide nonempty training and test partitions")
    specification = {"kind": args.kind, "alpha": args.alpha, "gamma": args.gamma}
    model = fit_physical_reader(features[training], [samples[i] for i in training], specification)
    predicted = model.predict(features[testing])
    result = pd.DataFrame(predicted, columns=TARGET_NAMES)
    result.insert(0, "sample_id", ids[testing])
    args.output.mkdir(parents=True, exist_ok=False)
    result.to_csv(args.output / "predictions.csv", index=False)
    joblib.dump(
        {"model": model, "feature_names": names.tolist(), "target_names": list(TARGET_NAMES)},
        args.output / "reader.joblib",
    )
    (args.output / "report.json").write_text(
        json.dumps(
            {
                "specification": specification,
                "train_samples": len(training),
                "test_samples": len(testing),
                "group_equal_physical_loss": validation_loss([samples[i] for i in testing], predicted),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
