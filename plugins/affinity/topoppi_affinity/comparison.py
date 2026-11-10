"""Compare external predictions on the exact saved test observations."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from topoppi_affinity.evaluation import regression_metrics


def paired_mae_difference(truth, topoppi, external, groups, *, repetitions=1000, seed=42):
    """Positive differences favor TopoPPI; resample independent split groups."""
    truth, topoppi, external = (np.asarray(x, dtype=float) for x in (truth, topoppi, external))
    differences = np.abs(external - truth) - np.abs(topoppi - truth)
    labels = np.asarray(groups).astype(str)
    if not len(truth) or topoppi.shape != truth.shape or external.shape != truth.shape or labels.shape != truth.shape:
        raise ValueError("Paired comparison requires aligned nonempty observations")
    if not np.isfinite(differences).all():
        raise ValueError("Paired comparison requires finite predictions")
    indices = [np.flatnonzero(labels == name) for name in np.unique(labels)]
    result = {"external_minus_topoppi_mae": float(differences.mean()), "n_groups": len(indices)}
    if len(indices) >= 2 and repetitions:
        rng = np.random.default_rng(seed)
        draws = [float(differences[np.concatenate([indices[i] for i in rng.integers(len(indices), size=len(indices))])].mean())
                 for _ in range(repetitions)]
        result["paired_group_bootstrap_95_ci"] = np.quantile(draws, [.025, .975]).tolist()
    return result


def read_external(path):
    path = Path(path)
    frame = pd.read_csv(path, keep_default_na=False)
    required = {"sample_id", "predicted", "target", "method", "status", "error"}
    if not required.issubset(frame):
        raise ValueError(f"{path}: missing columns {sorted(required - set(frame))}")
    if frame.sample_id.duplicated().any() or frame.method.nunique() != 1 or frame.target.nunique() != 1:
        raise ValueError(f"{path}: use one method/target and unique sample IDs per file")
    metadata = json.loads(path.with_suffix(".json").read_text())
    for key in ("protocol", "training_overlap", "source_url"):
        if key not in metadata:
            raise ValueError(f"{path}: provenance must specify {key}")
    if "training_overlap" not in frame:
        frame["training_overlap"] = metadata["training_overlap"]
    frame["predicted"] = pd.to_numeric(frame.predicted, errors="coerce")
    frame["valid_prediction"] = (frame.status == "ok") & np.isfinite(frame.predicted)
    return frame, metadata


def common_cohort_metrics(matched_records):
    """Use the intersection of successful IDs for a multi-method comparison."""
    rows = []
    for experiment, frame in matched_records.groupby("experiment", sort=False):
        methods = list(frame.method.unique())
        shared = set.intersection(*(set(frame.loc[frame.method == method, "sample_id"]) for method in methods))
        common = frame[frame.sample_id.isin(shared)]
        if common.empty:
            continue
        first = common[common.method == methods[0]]
        top_metrics = regression_metrics(first.observed, first.topoppi_prediction, first.complex_id)
        rows.append({"experiment": experiment, "method": "TopoPPI", "n_common": len(first),
                     "n_external_methods": len(methods), **top_metrics})
        for method, subset in common.groupby("method", sort=False):
            metrics = regression_metrics(subset.observed, subset.predicted, subset.complex_id)
            rows.append({"experiment": experiment, "method": method, "n_common": len(subset),
                         "n_external_methods": len(methods), **metrics})
    return pd.DataFrame(rows)


def compare_methods(run_dir, prediction_paths, output_dir, *, repetitions=1000, seed=42):
    run_dir, output_dir = Path(run_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = json.loads((run_dir / "summary.json").read_text())
    external = [read_external(path) for path in prediction_paths]
    rows, reports, matched_records = [], [], []
    for experiment, specification in summary["experiments"].items():
        subset = pd.read_csv(run_dir / experiment / "split.csv", keep_default_na=False)
        subset = subset[subset.used_test].copy()
        saved = pd.read_csv(run_dir / experiment / "predictions.csv")
        selected = saved[saved.model == specification["selected_topoppi_model"]]
        base = subset.merge(selected[["sample_id", "observed", "predicted"]], on="sample_id", validate="one_to_one")
        if len(base) != len(subset):
            raise ValueError(f"{experiment}: selected TopoPPI model lacks test predictions")
        base = base.rename(columns={"predicted": "topoppi_prediction"})
        split_unit = {"split_complex": "interaction_group", "split_sequence": "sequence_group"}.get(
            specification["split_column"], "complex_id")
        for candidates, provenance in external:
            method, target = candidates.method.iloc[0], candidates.target.iloc[0]
            if target != specification["target"]:
                continue
            attempted = candidates[candidates.sample_id.isin(base.sample_id)]
            valid = attempted[attempted.valid_prediction]
            matched = base.merge(valid[["sample_id", "predicted", "training_overlap"]], on="sample_id", validate="one_to_one")
            result = {"experiment": experiment, "method": method, "target": target,
                      "protocol": provenance["protocol"], "training_overlap": provenance["training_overlap"],
                      "n_requested": len(base), "n_attempted": len(attempted), "n_matched": len(matched),
                      "n_failed": int((~attempted.valid_prediction).sum()),
                      "n_not_attempted": len(base) - len(attempted), "bootstrap_unit": split_unit,
                      "topoppi_selected_model": specification["selected_topoppi_model"]}
            report = {**result, "provenance": provenance, "failure_records": attempted[~attempted.valid_prediction][
                ["sample_id", "status", "error"]].to_dict("records")}
            if len(matched):
                truth, top, other = matched.observed.to_numpy(), matched.topoppi_prediction.to_numpy(), matched.predicted.to_numpy()
                top_metrics = regression_metrics(truth, top, matched.complex_id)
                other_metrics = regression_metrics(truth, other, matched.complex_id)
                paired = paired_mae_difference(truth, top, other, matched[split_unit], repetitions=repetitions, seed=seed)
                report.update({"topoppi": top_metrics, "external": other_metrics, "paired": paired,
                               "overlap_counts": matched.training_overlap.value_counts().to_dict()})
                for prefix, metrics in [("topoppi", top_metrics), ("external", other_metrics)]:
                    result.update({f"{prefix}_{key}": metrics[key] for key in ("mae", "rmse", "macro_mae", "spearman", "pearson")})
                result["mae_difference"] = paired["external_minus_topoppi_mae"]
                interval = paired.get("paired_group_bootstrap_95_ci", [None, None])
                result["mae_difference_ci_low"], result["mae_difference_ci_high"] = interval
                matched["experiment"], matched["method"], matched["protocol"] = experiment, method, provenance["protocol"]
                matched_records.append(matched[["experiment", "method", "protocol", "sample_id", "complex_id", "interaction_group",
                    "sequence_group", "n_mutations", "observed", "topoppi_prediction", "predicted", "training_overlap"]])
            rows.append(result)
            reports.append(report)
    pd.DataFrame(rows).to_csv(output_dir / "comparison.csv", index=False)
    if matched_records:
        combined = pd.concat(matched_records, ignore_index=True)
        combined.to_csv(output_dir / "matched_predictions.csv", index=False)
        common_cohort_metrics(combined).to_csv(output_dir / "common_cohort_metrics.csv", index=False)
    result = {"scope": "external-method comparison on the supplied evaluation cohort", "seed": seed,
              "interpretation": "Paired differences compare identical test IDs. Released models with known, possible or unknown training overlap do not establish independent generalization or a fair retrained-method ranking.",
              "comparisons": reports}
    (output_dir / "comparison.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return result
