"""Run original EvoEF 1.1 repair/build/binding on an explicitly selected cohort."""

from __future__ import annotations

import json
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from topoppi_affinity.external_inputs import prepare_numbered_complex, remap_mutations


def parse_binding_energy(output):
    values = re.findall(r"^Total\s*=\s*([-+0-9.eE]+)\s*$", output, flags=re.MULTILINE)
    if len(values) != 1 or not np.isfinite(float(values[0])):
        raise ValueError("Expected exactly one finite group-split EvoEF binding energy")
    return float(values[0])


def _run_command(binary, directory, logfile, *arguments):
    result = subprocess.run([str(binary), *arguments], cwd=directory, capture_output=True, text=True, timeout=1800)
    (directory / logfile).write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"EvoEF exited {result.returncode}; see {directory / logfile}")
    return result.stdout


def _run_complex(binary, output_dir, frame):
    row = frame.iloc[0]
    case = output_dir / "cases" / row.complex_id.replace(":", "_").replace(",", "")
    cache = case / "predictions.csv"
    if cache.exists():
        cached = pd.read_csv(cache, keep_default_na=False)
        if cached.sample_id.tolist() == frame.sample_id.tolist() and (cached.status == "ok").all():
            return cached.to_dict("records")
    case.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    records = []
    try:
        mapping = prepare_numbered_complex(row, case / "wt.pdb")
        unique = dict.fromkeys(frame.mutations)
        translated = {value: remap_mutations(value, mapping) for value in unique}
        (case / "mutations.txt").write_text("".join(value + ";\n" for value in translated.values()))
        _run_command(binary, case, "repair.log", "--command=RepairStructure", "--pdb=wt.pdb", "--num_of_runs=3")
        _run_command(binary, case, "build.log", "--command=BuildMutant", "--pdb=wt_Repair.pdb",
                     "--mutant_file=mutations.txt", "--num_of_runs=10")
        split = row.chains_a.replace(",", "") + "," + row.chains_b.replace(",", "")
        values = {}
        for i, mutation in enumerate(unique, start=1):
            prefix = f"wt_Repair_Model_{i:04d}"
            wt = parse_binding_energy(_run_command(binary, case, prefix + "_WT.energy.log", "--command=ComputeBinding",
                                      f"--pdb={prefix}_WT.pdb", f"--split={split}"))
            mt = parse_binding_energy(_run_command(binary, case, prefix + ".energy.log", "--command=ComputeBinding",
                                      f"--pdb={prefix}.pdb", f"--split={split}"))
            values[mutation] = (mt - wt, wt, mt)
        for _, observation in frame.iterrows():
            ddg, wt, mt = values[observation.mutations]
            records.append({"sample_id": observation.sample_id, "predicted": ddg, "target": "ddg", "method": "EvoEF1",
                            "status": "ok", "error": "", "wt_binding_energy": wt, "mt_binding_energy": mt,
                            "mapped_mutations": translated[observation.mutations], "training_overlap": "possible"})
    except (ValueError, KeyError, OSError, RuntimeError, subprocess.TimeoutExpired) as exc:
        records = [{"sample_id": sample, "predicted": np.nan, "target": "ddg", "method": "EvoEF1",
                    "status": "error", "error": str(exc), "training_overlap": "possible"} for sample in frame.sample_id]
    pd.DataFrame(records).to_csv(cache, index=False)
    (case / "run.json").write_text(json.dumps({"complex_id": row.complex_id, "elapsed_seconds": time.monotonic() - start,
        "repair_rounds": 3, "build_rounds": 10, "reference": "paired WT emitted by original BuildMutant2"}, indent=2))
    return records


def run_evoef(manifest, binary, output_dir, *, workers=4, limit_complexes=None):
    """Run the supplied manifest rows, optionally limiting the number of complexes."""
    frame = pd.read_csv(manifest, keep_default_na=False)
    output_dir, binary = Path(output_dir), Path(binary).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    groups = [g.copy() for _, g in frame.groupby("complex_id", sort=True)]
    if limit_complexes is not None:
        groups = groups[:limit_complexes]
    records = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        jobs = [executor.submit(_run_complex, binary, output_dir, group) for group in groups]
        for completed, job in enumerate(as_completed(jobs), start=1):
            records.extend(job.result())
            print(f"EvoEF {completed}/{len(groups)} complexes; {len(records)} observations", flush=True)
    predictions = pd.DataFrame(records).sort_values("sample_id")
    predictions.to_csv(output_dir / "evoef.csv", index=False)
    provenance = {"method": "EvoEF1", "target": "ddg", "protocol": "frozen_empirical", "training_overlap": "possible",
        "source_url": "https://github.com/tommyhuangthu/EvoEF", "version": "1.1", "license": "MIT",
        "binary": str(binary), "manifest": str(Path(manifest).resolve()), "n_attempted": len(predictions),
        "cohort": "all_manifest_rows", "limit_complexes": limit_complexes,
        "n_success": int((predictions.status == "ok").sum()), "repair_rounds": 3, "build_rounds": 10,
        "workers": workers,
        "sign": "mutant binding energy minus paired WT binding energy; positive weaker",
        "units": "EvoEF empirical ddG scale, compared to experimental kcal/mol without test calibration",
        "numbering": "Observed canonical residues renumbered per chain; coordinates retained; complete declared partner groups",
        "reference_protocol": "Original current BuildMutant2 produces matched optimized mutant and WT; ComputeBinding --split on each",
        "notes": "EvoEF1 parameters were optimized using mutation thermodynamic datasets. No training-overlap exclusion is claimed. Native interaction energies are not treated as absolute experimental binding free energies."}
    (output_dir / "evoef.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return predictions
