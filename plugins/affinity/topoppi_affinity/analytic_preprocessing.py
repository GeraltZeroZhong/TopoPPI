"""Analytic-curvature feature preprocessing for supplied structures."""
from __future__ import annotations

import json
import subprocess
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from topoppi_affinity_tools import tool_file

from topoppi_affinity.analytic_curvature import update_pooled_curvature
from topoppi_affinity.embedding_baselines import zero_surface_slots
from topoppi_affinity.frozen_probes import load_features
from topoppi_affinity.graph_features import load_graph

CURVATURE = "analytic_Gaussian_Goldman"
PREPROCESSING = {
    "curvature": CURVATURE,
    "field": "native unit-peak Gaussian sum; sigma=1 Angstrom; active axis-aligned 4-sigma atom support",
    "geometry": "same native unsmoothed mesh, vertex area and residue footprint identity",
    "derivatives": "analytic gradient and Hessian; Goldman implicit Gaussian curvature formula",
    "replaced_graph_columns": ["surface.gaussian_curvature_mean", "surface.gaussian_curvature_std"],
    "pooled_readout": "replace curvature entries in surface.chemical_delta_context; keep all ESM-C and other slots",
    "command": "python -m topoppi_affinity_tools.features.prepare_analytic_mutation_inputs --input INPUT.csv --graph-manifest GRAPHS.csv --features ORIGINAL_SURFACE.npz --contact-features ORIGINAL_CONTACT.npz --output OUTPUT_DIRECTORY",
    "model_input": "OUTPUT_DIRECTORY/ddg_plm_analytic.npz; original discrete-curvature surface vectors are a different protocol",
}


def prepare_analytic_mutation_inputs(input_manifest, graph_manifest, features, contact_features, output):
    """Only structure/mutation identity and unlabeled pooled vectors are read."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(input_manifest, usecols=["sample_id", "source", "complex_id", "mutations"], keep_default_na=False)
    graphs = pd.read_csv(graph_manifest, keep_default_na=False)
    graphs = frame[["source", "complex_id"]].drop_duplicates().merge(
        graphs, on=["source", "complex_id"], validate="one_to_one")
    graphs.to_csv(output / "input_graphs.csv", index=False)
    subprocess.run([sys.executable, tool_file('features/refine_graph_curvature.py'), "--manifest",
                    str(output / "input_graphs.csv"), "--output", str(output / "graphs"), "--workers", "2"], check=True)
    refined = pd.read_csv(output / "graphs/manifest.csv", keep_default_na=False)
    paths = {(row.source, row.complex_id): row.graph_path for row in refined.itertuples() if row.status == "ok"}
    original, names = load_features(features, frame)
    contact, contact_names = load_features(contact_features, frame)
    if not np.array_equal(names, contact_names):
        raise ValueError("Supplied original surface/contact schemas differ")

    @lru_cache(maxsize=16)
    def graph(key):
        return load_graph(paths[key])

    values, selected, coverage = [], [], []
    for index, row in enumerate(frame.itertuples()):
        record = {"sample_id": row.sample_id, "source": row.source, "complex_id": row.complex_id}
        try:
            values.append(update_pooled_curvature(original[index], names, graph((row.source, row.complex_id)),
                                                  target="ddg", mutations=row.mutations))
            selected.append(index)
            record.update(status="ok", error="")
        except (ValueError, KeyError, OSError) as error:
            record.update(status="error", error=f"{type(error).__name__}: {error}")
        coverage.append(record)
    matrix = np.asarray(values)
    difference = float(np.max(np.abs(zero_surface_slots(matrix, names) - contact[selected])))
    if difference != 0.:
        raise ValueError("Analytic curvature changed the original contact control")
    np.savez_compressed(output / "ddg_plm_analytic.npz", X=matrix,
                        sample_ids=frame.iloc[selected].sample_id.to_numpy(dtype=str),
                        feature_names=names, curvature_protocol=np.asarray(CURVATURE))
    pd.DataFrame(coverage).to_csv(output / "coverage.csv", index=False)
    report = {**PREPROCESSING, "input_manifest": str(input_manifest), "original_surface_features": str(features),
              "n_graphs": len(graphs), "successful_graphs": len(paths), "n_rows": len(frame),
              "successful_rows": len(selected), "dimension": matrix.shape[1],
              "changed_pooled_columns": int(np.any(original[selected] != matrix, axis=0).sum()),
              "contact_max_absolute_difference": difference, "labels_read": False, "model_prediction_run": False}
    (output / "protocol.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    return report


