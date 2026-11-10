#!/usr/bin/env python3
"""Extract native TopoPPI descriptors and fixed-WT mutation controls."""
from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np

from topoppi_affinity.features import extract_complex_features, feature_group, mutation_features


def write_npz(path, rows, names, **extra):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, sample_ids=np.asarray([r[0] for r in rows], dtype=str),
                        X=np.asarray([r[1] for r in rows], dtype=np.float32).reshape(len(rows), len(names)),
                        feature_names=np.asarray(names, dtype=str),
                        feature_groups=np.asarray([feature_group(k) for k in names], dtype=str), **extra)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cache-dir", required=True, type=Path)
    parser.add_argument("--limit", type=int, help="Optional diagnostic row limit; not a training split.")
    args = parser.parse_args()
    with args.manifest.open(newline="") as handle:
        records = list(csv.DictReader(handle))
    if args.limit:
        records = records[:args.limit]
    mutation_task = any(str(r.get("mutations", "")).strip() for r in records)
    rows, singleton_rows, singleton_parents, singleton_mutations = [], [], [], []
    failures, coverage, timings, geometry_warnings, names = [], [], [], [], None
    started = time.monotonic()
    previous_key, previous_bundle = None, None
    for index, record in enumerate(records):
        start = time.monotonic()
        sample_id = record["sample_id"]
        try:
            structure_key = (record["structure_path"], record["chains_a"], record["chains_b"])
            if structure_key != previous_key:
                previous_bundle = extract_complex_features(*structure_key, args.cache_dir)
                previous_key = structure_key
            bundle = previous_bundle
            features = dict(bundle["absolute_features"])
            empty_sides = [i for i, side in enumerate(bundle["metadata"]["sides"])
                           if side["accepted_patch_count"] == 0]
            if empty_sides:
                geometry_warnings.append({"sample_id": sample_id, "empty_interface_sides": empty_sides,
                                          "interface_area": features["surface.interface_area_sum"]})
            mutations = record.get("mutations", "").strip()
            if mutation_task:
                features.update(mutation_features(bundle, mutations))
            if names is None:
                names = list(features)
            if list(features) != names:
                raise ValueError("Feature schema differs across records.")
            values = [features[k] for k in names]
            if not np.isfinite(values).all():
                raise ValueError("Non-finite feature value.")
            rows.append((sample_id, values))
            if mutation_task:
                for mutation in mutations.split(","):
                    if not mutation.strip():
                        continue
                    single = dict(bundle["absolute_features"])
                    single.update(mutation_features(bundle, mutation))
                    singleton_rows.append((f"{sample_id}:{mutation}", [single[k] for k in names]))
                    singleton_parents.append(sample_id)
                    singleton_mutations.append(mutation)
                coverage.append({"sample_id": sample_id, "complex_id": record["complex_id"],
                    "n_mutations": features["mutation.n_mutations"],
                    "surface_covered_count": features["mutation.surface_covered_count"],
                    "interface_count": features["mutation.contact_interface_count"],
                    "near_count": features["mutation.contact_near_count"],
                    "distal_count": features["mutation.contact_distal_count"],
                    "geometry_protocol": "fixed_WT_structure_identity_substitution"})
            timings.append({"sample_id": sample_id, "seconds": round(time.monotonic() - start, 3)})
        except (ValueError, KeyError, OSError, IndexError) as exc:
            failures.append({"sample_id": sample_id, "complex_id": record.get("complex_id", ""),
                             "error": str(exc)})
        if index == 0 or (index + 1) % 10 == 0 or index + 1 == len(records):
            print(f"{index + 1}/{len(records)} samples: {len(rows)} extracted, {len(failures)} rejected; "
                  f"elapsed {time.monotonic() - started:.1f}s", flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary = {"manifest": str(args.manifest.resolve()), "samples_requested": len(records),
               "samples_extracted": len(rows), "samples_rejected": len(failures),
               "features": len(names or []), "elapsed_seconds": round(time.monotonic() - started, 3),
               "representation": "native_Gaussian_surface_3D_footprints_no_UV_optimization",
               "geometry_protocol": ("fixed_WT_structure_identity_substitution" if mutation_task else
                                     "supplied_complex_native_Gaussian_surface"), "timings": timings,
               "failures": failures, "geometry_warnings": geometry_warnings}
    args.output.with_suffix(".json").write_text(json.dumps(summary, indent=2))
    if not rows:
        raise SystemExit(f"No valid samples; see {args.output.with_suffix('.json')}")
    write_npz(args.output, rows, names)
    if mutation_task:
        write_npz(args.output.with_name(args.output.stem + "_singletons.npz"), singleton_rows, names,
                  parent_sample_ids=np.asarray(singleton_parents, dtype=str),
                  mutations=np.asarray(singleton_mutations, dtype=str))
        with args.output.with_name(args.output.stem + "_coverage.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(coverage[0]))
            writer.writeheader()
            writer.writerows(coverage)
    print(json.dumps({k: v for k, v in summary.items() if k not in {"timings", "failures", "geometry_warnings"}}, indent=2))


if __name__ == "__main__":
    main()
