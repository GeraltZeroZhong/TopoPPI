"""Label-free inference using fitted task models."""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from topoppi_affinity.features import feature_group, mutation_features
from topoppi_affinity.model_metadata import training_metadata
from topoppi_affinity.models import AffinityRegressor
from topoppi_affinity.runtime import check_cancel
from topoppi_affinity.thermodynamics import R_KCAL, kd_to_pkd, pkd_to_dg, pkd_to_kd


def model_info(bundle):
    from topoppi_affinity.production import PhysicalRegressor
    if isinstance(bundle, PhysicalRegressor):
        return {"kind": "production", "task": bundle.metadata["task"], "target": bundle.target,
                "scope": bundle.metadata["scope"], "curvature": "none",
                "training": training_metadata(bundle.metadata["training"]),
                "model_version": bundle.metadata["model_version"],
                "feature_scheme": bundle.metadata["feature_scheme"]}
    if isinstance(bundle, AffinityRegressor):
        return {"kind": "descriptor", "target": bundle.target, "training": training_metadata(bundle.training_metadata),
                "curvature": "none", "scope": "Conventional interface-descriptor regression"}
    if not isinstance(bundle, dict) or not isinstance(bundle.get("report"), dict):
        raise ValueError("This file is not a supported TopoPPI affinity model; choose an installed task or a compatible model artifact")
    report = bundle["report"]
    if "residual" in bundle and "base" in bundle:
        kind, scope = "known", "Mutation-effect ensemble for systems represented in training"
    elif "base" in bundle and "prior_method" in bundle:
        kind, scope = "unseen", "Mutation-effect ensemble for transfer to unrepresented systems"
    elif "plm_feature_names" in bundle and "legacy_feature_names" in bundle:
        kind, scope = "absolute", "Absolute-affinity regression with sequence and interface descriptor kernels"
    elif "model" in bundle and "feature_names" in bundle:
        kind, scope = "fitted", "Custom fitted affinity regressor"
    else:
        raise ValueError("Unsupported deployment artifact; use an embedded absolute-affinity or mutation-effect model")
    training = training_metadata(report)
    analytic = (report.get("curvature") == "analytic_Gaussian_Goldman" or
                (kind == "fitted" and training.get("representation") == "sequence_and_analytic_surface"))
    curvature = "analytic_gaussian" if analytic else "angle_defect"
    return {"kind": kind, "target": report["target"], "training": training, "scope": scope,
            "curvature": curvature, "training_scope_note": "Training membership and counts describe this fitted estimator; cross-validation metrics describe the evaluation protocol"}


def add_units(frame, target, temperature=298.15, wt_kd=None):
    if not np.isfinite(temperature) or temperature <= 0:
        raise ValueError("Temperature must be finite positive kelvin")
    frame = frame.copy()
    frame["target"] = target
    if target == "pkd":
        if wt_kd is not None:
            raise ValueError("WT Kd anchoring applies to a mutation model")
        frame["predicted_pkd"] = frame.predicted
        frame["predicted_kd_M"] = pkd_to_kd(frame.predicted.to_numpy())
        frame["predicted_dg_at_reference_kcal_mol"] = pkd_to_dg(frame.predicted.to_numpy(), temperature)
        frame["reference_temperature_K"] = temperature
    else:
        frame["predicted_ddg_kcal_mol"] = frame.predicted
        if wt_kd is not None:
            anchor = float(kd_to_pkd(wt_kd))
            frame["anchor_wt_kd_M"] = wt_kd
            frame["anchored_mutant_pkd"] = anchor - frame.predicted / (R_KCAL * temperature * np.log(10))
            frame["anchored_mutant_kd_M"] = pkd_to_kd(frame.anchored_mutant_pkd.to_numpy())
            frame["conversion_temperature_K"] = temperature
            frame["affinity_mode"] = "measured_WT_anchor_plus_predicted_ddg"
    return frame


def predict_features(model, *, features=None, plm_features=None, legacy_features=None, prior=None,
                     temperature=298.15, wt_kd=None):
    """Infer from named NPZ arrays and optional explicit prior CSV, without labels."""
    from topoppi_affinity.artifacts import predict_frozen
    from topoppi_affinity.legacy_kernel import predict_legacy_kernel
    from topoppi_affinity.mutation_ensemble import predict_mutation_ensemble, predict_scalar_prior_ensemble
    from topoppi_affinity.production import load_model
    bundle = load_model(model)
    info = model_info(bundle)
    kind = info["kind"]
    required = {"absolute": {"plm-features": plm_features, "descriptor-features": legacy_features},
                "known": {"plm-features": plm_features, "descriptor-features": legacy_features},
                "unseen": {"plm-features (or features)": plm_features or features}}
    for label, path in required.get(kind, {"features": features}).items():
        if path is None:
            raise ValueError(f"The selected {kind} model requires --{label}")
        if not Path(path).expanduser().is_file():
            raise ValueError(f"Feature file does not exist: {path}")
    if kind == "absolute":
        result = predict_legacy_kernel(model, plm_features, legacy_features, temperature_K=temperature)
    elif kind == "known":
        result = predict_mutation_ensemble(bundle, legacy_features, plm_features, prior)
    elif kind == "unseen":
        result = predict_scalar_prior_ensemble(bundle, plm_features or features, prior)
    elif kind in {"descriptor", "production"}:
        with np.load(features, allow_pickle=False) as values:
            predicted = bundle.predict(values["X"], values["feature_names"])
            result = pd.DataFrame({"sample_id": values["sample_ids"].astype(str), "predicted": predicted})
    else:
        result = predict_frozen(model, features, temperature_K=temperature)
    if "status" not in result:
        result["status"] = "ok"
    if "error" not in result:
        result["error"] = ""
    return add_units(result, info["target"], temperature, wt_kd)


def _save_features(path, rows, names, ids, **extra):
    np.savez_compressed(path, X=np.asarray(rows), feature_names=np.asarray(names, dtype=str),
        sample_ids=np.asarray(ids, dtype=str), **extra)
    return path


def _pooled_vectors(graph, frame, sequence_dir, *, target, use_surface=True):
    from topoppi_affinity.embedding_baselines import pooled_features
    from topoppi_affinity.graph_features import mutant_node_features

    wild = pd.read_csv(sequence_dir / "wild_type_chains.csv", keep_default_na=False, dtype={"chain": str, "source": str, "complex_id": str})
    mutant = pd.read_csv(sequence_dir / "mutant_chains.csv", keep_default_na=False, dtype={"chain": str, "sample_id": str})
    mutations = {(row.sample_id, row.chain): row.sequence_id for row in mutant.itertuples()}
    row = frame.iloc[0]
    wild = wild[(wild.source == row.source) & (wild.complex_id == row.complex_id)]
    wt = np.zeros((len(graph["x"]), 1152), dtype=np.float32)
    mappings = {}
    for chain in np.unique(graph["chain"]):
        record = wild[wild.chain == chain].iloc[0]
        keys = json.loads(record.residue_keys)
        positions = np.flatnonzero(graph["chain"] == chain)
        lookup = {key: index for index, key in enumerate(keys)}
        indices = np.array([lookup[key] for key in graph["residue_key"][positions]])
        with np.load(sequence_dir / "embeddings" / f"{record.sequence_id}.npz", allow_pickle=False) as source:
            if str(source["sequence"]) != record.sequence or str(source["model_name"]) != "esmc_600m":
                raise ValueError("ESMC embedding sequence/model differs from the structural sequence")
            wt[positions] = source["embedding"][indices]
        mappings[chain] = (positions, indices)
    rows, names = [], None
    for row in frame.itertuples():
        mt = wt.copy()
        for chain, (positions, indices) in mappings.items():
            identifier = mutations.get((row.sample_id, chain))
            if identifier:
                with np.load(sequence_dir / "embeddings" / f"{identifier}.npz", allow_pickle=False) as source:
                    if str(source["model_name"]) != "esmc_600m":
                        raise ValueError("Mutant embedding is not ESMC-600M")
                    mt[positions] = source["embedding"][indices]
        x_mut, _, mask = mutant_node_features(graph, row.mutations)
        data = SimpleNamespace(x=graph["x"], x_mut=x_mut, esm=wt, esm_mut=mt,
            side=graph["side"], contact_mask=graph["x"][:, 26] <= 10., mutation_mask=mask)
        vector, names = pooled_features(data, target=target, use_surface=use_surface)
        rows.append(vector)
    return np.asarray(rows), np.asarray(names)


def build_features(prepared, descriptor_bundle, frame, model_bundle, output, *, sequence_dir=None,
                   device="cpu", log=print, cancel_event=None, shared_embeddings=None, event_callback=None, cache_dir=None):
    """Use the shared mesh and exact frozen protocols, not GUI display settings."""
    info = model_info(model_bundle)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if info["kind"] == "production":
        from topoppi_affinity.production import physical_feature_rows
        check_cancel(cancel_event)
        rows, names = physical_feature_rows(prepared, descriptor_bundle, frame, model_bundle)
        path = _save_features(output / "descriptors.npz", rows, names, frame.sample_id)
        return {"features": path}, {"feature_scheme": info["feature_scheme"],
                                   "geometry": "fixed_input_complex", "uv_used": False}
    rows = []
    for row in frame.itertuples():
        values = dict(descriptor_bundle["absolute_features"])
        if info["target"] == "ddg":
            values.update(mutation_features(descriptor_bundle, row.mutations))
        names = list(values)
        rows.append([values[name] for name in names])
    paths = {"legacy_features": _save_features(output / "descriptors.npz", rows, names, frame.sample_id,
        feature_groups=np.asarray([feature_group(name) for name in names]))}
    representation = info["training"].get("representation", "")
    if info["kind"] == "descriptor" or (info["kind"] == "fitted" and representation == "conventional_descriptors"):
        paths["features"] = paths["legacy_features"]
        return paths, {"geometry": "native_gaussian_interface", "sequence": "residue composition", "uv_used": False}

    from topoppi_affinity.graph_features import extract_complex_graph, save_graph
    from topoppi_affinity.sequence_embeddings import extract_embeddings, prepare_sequences

    if cancel_event is not None and cancel_event.is_set():
        raise InterruptedError("Analysis cancelled")
    log("Computing native interface graph descriptors from the shared native interface")
    graph = extract_complex_graph(prepared.source["path"], *[p.chains for p in prepared.partners],
        curvature_method=info["curvature"], prepared=prepared, descriptor_bundle=descriptor_bundle)
    graph_metadata = json.loads(str(graph["metadata_json"]))
    graph_metadata["protocol"] = "native_residue_surface_graph"
    graph_metadata["source"]["protocol"] = "native_gaussian_interface"
    graph["metadata_json"] = np.asarray(json.dumps(graph_metadata))
    save_graph(graph, output / "graph.npz")
    graph_index = pd.DataFrame([{"source": frame.iloc[0].source, "complex_id": frame.iloc[0].complex_id,
        "graph_path": str((output / "graph.npz").resolve()), "status": "ok"}])
    graph_index.to_csv(output / "graphs.csv", index=False)
    sequence_dir = Path(sequence_dir) if sequence_dir is not None else output / "sequences"
    prepare_sequences([output.parent / "inputs.csv"], sequence_dir, prepared=prepared)
    wt = pd.read_csv(sequence_dir / "wild_type_chains.csv", keep_default_na=False, dtype={"chain": str, "source": str, "complex_id": str})
    mt = pd.read_csv(sequence_dir / "mutant_chains.csv", keep_default_na=False, dtype={"chain": str, "sample_id": str})
    first = frame.iloc[0]
    ids = set(wt[(wt.source == first.source) & (wt.complex_id == first.complex_id)].sequence_id)
    ids.update(mt[mt.sample_id.isin(frame.sample_id)].sequence_id)
    log("Reading or computing ESMC-600M embeddings")
    provenance = extract_embeddings(sequence_dir, device=device, sequence_ids=ids,
        shared_cache=shared_embeddings, cancel_event=cancel_event, progress=event_callback, log=log, cache_dir=cache_dir)
    check_cancel(cancel_event)
    custom_representation = representation if info["kind"] == "fitted" else ""
    matrix, names = _pooled_vectors(graph, frame, sequence_dir, target=info["target"],
                                    use_surface=custom_representation != "sequence")
    if custom_representation == "sequence_and_contact":
        from topoppi_affinity.embedding_baselines import zero_surface_slots

        matrix = zero_surface_slots(matrix, names)
    paths["plm_features"] = _save_features(output / "plm.npz", matrix, names, frame.sample_id,
        curvature_protocol=np.asarray(info["curvature"]))
    paths["features"] = paths["plm_features"]
    return paths, {"geometry": "native_gaussian_interface", "curvature": info["curvature"],
        "sequence": provenance, "contact_statistics_cutoff_A": 6., "ESMC_contact_pooling_cutoff_A": 10.,
        "uv_used": False, "structure_protocol": "fixed_WT_identity_substitution"}


def compute_prior(frame, info, output, *, prior_csv=None, repository=None, checkpoint=None, device="cpu",
                  cancel_event=None, log=print):
    """Use the specified frozen prior, preserving the predeclared base-only fallback."""
    kind = info["kind"]
    check_cancel(cancel_event)
    if kind not in {"known", "unseen"}:
        return None, {"status": "not_required"}
    required = (info["training"]["residual"]["prior_method"] if kind == "known"
                else info["training"]["prior_method"])
    if prior_csv is not None:
        prior = pd.read_csv(prior_csv, keep_default_na=False, dtype={"sample_id": str})
        if set(prior.sample_id) != set(frame.sample_id) or prior.sample_id.duplicated().any():
            raise ValueError("Supplied prior CSV must match the analysis sample IDs exactly")
        provenance = {"status": "supplied", "path": str(Path(prior_csv).resolve()), "method": required}
    elif repository is not None:
        try:
            if kind == "known":
                from topoppi_affinity.external_mutation import run_stabddg

                prior = run_stabddg(output / "inputs.csv", repository, output / "prior",
                    checkpoint_name="stability_finetuned.pt", device=device, mc_samples=20, records_frame=frame,
                    cancel_event=cancel_event, log=log)
            else:
                from topoppi_affinity.ba_cycle_prior import run_ba_cycle

                if checkpoint is None:
                    raise ValueError("BA-Cycle requires --prior-checkpoint (soluble ProteinMPNN v_48_020)")
                prior = run_ba_cycle(output, output / "prior", repository, checkpoint, device=device,
                    mc_samples=8, seed=42, manifest=output / "inputs.csv", graph_index=output / "features/graphs.csv",
                    cancel_event=cancel_event, log=log)
            provenance = {"status": "computed", "method": required, "repository": str(repository)}
        except InterruptedError:
            raise
        except (OSError, ValueError, RuntimeError, ImportError) as exc:
            prior = None
            provenance = {"status": "failed", "method": required, "error": f"{type(exc).__name__}: {exc}"}
    else:
        prior = None
        provenance = {"status": "not_provided", "method": required,
                      "error": "Prior not requested; predeclared base_only prediction"}
    if prior is None:
        prior = pd.DataFrame({"sample_id": frame.sample_id, "predicted": np.nan,
            "target": "ddg", "method": required, "status": "error", "error": provenance["error"]})
    if not prior.method.eq(required).all() or not prior.target.eq("ddg").all():
        raise ValueError(f"This model requires the prior {required} in kcal/mol")
    path = output / "prior.csv"
    prior.to_csv(path, index=False)
    return path, provenance
