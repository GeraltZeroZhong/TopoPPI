"""Alanine-scan ranking from shared mutation physics and decoded UV surfaces."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Data.PDBData import protein_letters_3to1


def rank_mutation_surface(prepared, selected, model_dir, output_dir, *, optcuts_bin=None, context=None):
    """Predict WT-to-Ala binding changes using the bundled surface regressor."""
    from topoppi_affinity.data import parse_mutations
    from topoppi_affinity.features import extract_complex_features
    from topoppi_affinity.mutation_physics import AtomicQueryContext, compute_residue_context
    from topoppi_affinity.production import mutation_physical_features

    from topoppi import save_interface
    from topoppi.residues import source_atom_residue_labels
    from topoppi_representation import prepare_spatial_descriptors

    from .mutation_model import MutationSurfaceRegressor

    def progress(stage, message):
        if context is not None:
            context.check_cancelled()
            context.emit(stage, message=message)

    model = MutationSurfaceRegressor(Path(model_dir) / "model.json")
    groups = tuple(sorted(partner.chains for partner in prepared.partners))
    source = prepared.source["path"]
    progress("features", "Calculating shared mutation descriptors and solvent accessibility.")
    bundle = extract_complex_features(source, *groups, prepared=prepared)
    records = bundle["residues"]
    interface_keys = [key for key, row in records.items() if row["contact_degree"] > 0]
    selected_keys = [key for key in interface_keys if records[key]["chain"] in selected]
    candidates = [key for key in selected_keys if records[key]["aa"] in "CDEFGHIKLMNPQRSTVWY"
                  and key.split(":")[-2] in protein_letters_3to1]
    if not candidates:
        raise ValueError("The selected partner has no standard non-alanine residues within the 6 A interface.")
    atomic_table = compute_residue_context(prepared)
    atomic_context = AtomicQueryContext(atomic_table)
    progress("encoding", "Preparing the residue-aware UV representation.")
    spatial, encoded, encoding_report = prepare_spatial_descriptors(
        prepared, sasa_residues=atomic_table,
        budget_ratio=model.metadata["encoding"]["budget_ratio"], optcuts_bin=optcuts_bin,
        cancel_event=context.cancel_event if context is not None else None,
        log=lambda message: progress("encoding", message),
    )
    has_ca = set()
    for partner in prepared.partners:
        for atom, key in zip(partner.atoms, source_atom_residue_labels(partner.atoms), strict=True):
            if atom.get_name() == "CA":
                has_ca.add(str(key))
    feature_rows, rows = [], []
    temperature = float(model.metadata.get("temperature_K", 298.15))
    for key in candidates:
        progress("prediction", f"Predicting alanine substitution at {key}.")
        row = records[key]
        mutation = f"{row['chain']}:{row['aa']}{row['residue_number']}{row['insertion']}A"
        values = mutation_physical_features(bundle, mutation, temperature, atomic_context=atomic_context)
        z, names, coverage = spatial.describe(parse_mutations(mutation))
        values.update(zip(names, z, strict=True))
        feature_rows.append([values[name] for name in model.feature_names])
        rows.append({"complex_id": prepared.source.get("complex_id", prepared.analysis_id),
                     "residue_key": key, "chain": row["chain"], "residue_number": row["residue_number"],
                     "insertion_code": row["insertion"], "side": groups.index(selected), "wt": row["aa"],
                     "mutations": mutation, "footprint_area_A2": row["footprint_area"], "has_CA": key in has_ca,
                     "spatial_coverage": "missing" if coverage["missing_sites"] else "covered"})
    X = np.asarray(feature_rows, dtype=float)
    predictions = model.predict(X, model.feature_names)
    if not np.isfinite(predictions).all():
        raise ValueError("The mutation-derived model returned a nonfinite prediction.")
    result = pd.DataFrame(rows)
    result["predicted_ddg_kcal_mol"] = predictions
    result["rank_score"] = predictions
    result["rank_min"] = result.rank_score.rank(method="min", ascending=False).astype(int)
    result["rank_max"] = result.rank_score.rank(method="max", ascending=False).astype(int)
    result["full_candidate_percentile"] = (result.rank_score.rank(method="average") - .5) / len(result)
    result["flags"] = [";".join(flag for flag, present in (
        ("missing_native_footprint", row.footprint_area_A2 <= 0), ("missing_CA", not row.has_CA),
        ("missing_spatial_context", row.spatial_coverage == "missing"),
        ("tied_score", row.rank_min != row.rank_max)) if present) for row in result.itertuples()]
    result = result.sort_values(["rank_min", "residue_key"]).reset_index(drop=True)
    coverage = {"partner_chains": ",".join(selected),
                "n_structure_residues": sum(row["chain"] in selected for row in records.values()),
                "n_WT_interface_context_nodes": len(selected_keys), "n_ranked_candidates": len(result),
                "n_native_alanine_context_nodes": sum(records[key]["aa"] == "A" for key in selected_keys),
                "n_other_excluded_interface_nodes": sum(key not in candidates and records[key]["aa"] != "A"
                                                        for key in selected_keys),
                "n_candidates_without_native_footprint": int((result.footprint_area_A2 <= 0).sum()),
                "n_candidates_without_CA": int((~result.has_CA).sum()),
                "n_candidates_without_spatial_context": int(result.spatial_coverage.eq("missing").sum())}
    output = Path(output_dir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    files = {"rankings": "rankings.csv", "candidate_coverage": "candidate_coverage.csv",
             "features": "features.csv", "representation": "representation.bin",
             "interface": "interface.npz", "metadata": "metadata.json"}
    metadata = {
        "model_dir": str(Path(model_dir).resolve()), "model_kind": "mutation_surface",
        "model_id": model.metadata["model_id"], "model_version": model.metadata["version"],
        "model_scope": model.metadata.get("scope"), "representation": model.metadata["representation"],
        "structure": str(source), "complex_id": rows[0]["complex_id"],
        "selected_partner": list(selected), "partner_groups": [list(group) for group in groups],
        "feature_protocol": model.metadata["feature_protocol"], "feature_protocol_source": "model.json",
        "prepared_interface": prepared.metadata(), "native_interface": "interface.npz",
        "feature_columns": list(model.feature_names), "encoding": encoding_report,
        "target": "ddg", "units": "kcal/mol", "temperature_K": temperature,
        "score_type": "predicted alanine-substitution binding delta-delta G (kcal/mol); positive means weaker binding",
        "candidate_coverage": coverage, "n_full_WT_interface_context_nodes": len(interface_keys),
        "candidate_population": "standard non-alanine WT residues of the selected partner within 6 A",
        "feature_context_population": "complete selected WT partners and their decoded interface surface",
        "rank_intervals": "descending minimum/maximum ranks within exact score ties",
        "full_candidate_percentile": "ascending midrank percentile (rank-.5)/candidate_count within partner",
        "labels_used": False, "affinity_plugin_used": True, "embeddings": None, "files": files,
    }
    progress("saving", "Saving predictions, descriptors and interface representation.")
    result.to_csv(output / files["rankings"], index=False)
    pd.DataFrame([coverage]).to_csv(output / files["candidate_coverage"], index=False)
    features = pd.DataFrame(X, columns=model.feature_names, index=pd.Index(candidates, name="residue_key"))
    features.loc[result.residue_key].to_csv(output / files["features"])
    (output / files["representation"]).write_bytes(encoded.blob)
    save_interface(prepared, output / files["interface"])
    (output / files["metadata"]).write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return result, metadata
