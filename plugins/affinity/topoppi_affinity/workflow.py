"""Shared-interface analysis, independent prediction/layout branches and review."""
from __future__ import annotations

import importlib.util
import json
import os
import uuid
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
import pandas as pd
from topoppi import __version__ as core_version
from topoppi.config import DEFAULT_RUN_CONFIG
from topoppi.interface import chain_group, load_partner_data, prepare_interface
from topoppi.interface_io import load_interface, save_interface
from topoppi.io.io_loader import PDBLoader

from topoppi_affinity.features import extract_complex_features, mutation_features, parse_mutations
from topoppi_affinity.runtime import check_cancel, resolve_device

RESIDUE_COLUMNS = ["analysis_id", "residue_key", "side", "chain", "residue_number", "insertion", "wt",
                   "contact_degree", "contact_atom_pairs", "partner_distance_A", "surface_area_A2",
                   "footprint_area_A2", "footprint_patch_count", "coverage", "source"]


def _plain(value):
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_plain(v) for v in value]
    if isinstance(value, np.generic):
        return _plain(value.item())
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _write_report(report, output):
    pending = output / ".analysis.json.tmp"
    pending.write_text(json.dumps(_plain(report), indent=2, allow_nan=False) + "\n")
    pending.replace(output / "analysis.json")


def _software_versions():
    versions = {}
    for name in ("numpy", "scipy", "biopython", "trimesh", "scikit-learn", "joblib", "torch", "esm"):
        try:
            versions[name] = version(name)
        except PackageNotFoundError:
            pass
    return versions


def variants_table(mutations=(), variants=None):
    """Each row is a substitution set; multiple substitutions in one row stay joint."""
    if variants is not None:
        frame = pd.read_csv(variants, keep_default_na=False, dtype=str) if isinstance(variants, (str, Path)) else variants.copy()
        if "variant_id" not in frame:
            raise ValueError("Variant CSV needs a variant_id column")
        if "mutations" not in frame:
            required = {"chain", "position", "wt", "mt"}
            if not required.issubset(frame.columns):
                raise ValueError("Variant CSV needs variant_id,mutations or variant_id,chain,position,wt,mt (optional insertion)")
            if "insertion" not in frame:
                frame["insertion"] = ""
            frame["mutations"] = [f"{r.chain}:{r.wt}{r.position}{r.insertion}{r.mt}" for r in frame.itertuples()]
            frame = frame.groupby("variant_id", sort=False, as_index=False).mutations.agg(",".join)
        frame = frame[["variant_id", "mutations"]].fillna("").astype(str).copy()
    else:
        mutations = [mutations] if isinstance(mutations, str) else list(mutations)
        frame = pd.DataFrame([{"variant_id": f"variant_{i+1:03d}", "mutations": m} for i, m in enumerate(mutations)])
        if frame.empty:
            frame = pd.DataFrame([{"variant_id": "complex", "mutations": ""}])
    if frame.empty or frame.variant_id.duplicated().any() or frame.variant_id.str.strip().eq("").any():
        raise ValueError("Variant IDs must be nonempty and unique")
    frame["sample_id"] = frame.variant_id
    return frame


def _structure_residues(atoms):
    from Bio.Data.PDBData import protein_letters_3to1_extended

    records = {}
    for atom in atoms:
        residue = atom.get_parent()
        key = (str(residue.get_parent().id), int(residue.id[1]), str(residue.id[2]).strip())
        records[key] = {"chain": key[0], "position": key[1], "insertion": key[2],
                        "wt": protein_letters_3to1_extended.get(residue.resname, "X")}
    return list(records.values())


def inspect_structure(path):
    """Read chain and author-residue identities for input controls, without surfaces."""
    loader = PDBLoader(str(Path(path).expanduser().resolve()))
    chains, residues = [], []
    for chain in loader.get_protein_chain_ids():
        _, atoms = loader.get_chain_data(chain)
        current = _structure_residues(atoms)
        chains.append({"id": chain, "residue_count": len(current)})
        residues.extend(current)
    return {"chains": chains, "residues": residues}


def _unused_directory(path):
    candidate, index = path, 2
    while candidate.exists():
        candidate = path.with_name(f"{path.name}-{index}")
        index += 1
    return candidate


def validate_inputs(structure, chains_a, chains_b, model=None, output=None, *, mutations=(), variants=None,
                    temperature=298.15, wt_kd=None, device="auto", task=None, prior_csv=None,
                    prior_repository=None, prior_checkpoint=None, mutation_context=None,
                    layout=False, optcuts_bin="OptCuts_bin", layout_timeout=600., retry=False, cache_dir=None, prepared=None, **_):
    """Validate cheap inputs before surfaces/embedding; retain individual invalid rows."""

    from topoppi_affinity.inference import model_info
    from topoppi_affinity.resources import resolve_model

    structure = Path(structure).expanduser().resolve()
    if prepared is None and not structure.is_file():
        raise ValueError(f"Structure file does not exist: {structure}")
    if not np.isfinite(temperature) or temperature <= 0:
        raise ValueError("Temperature must be finite positive kelvin")
    if wt_kd is not None and (not np.isfinite(wt_kd) or wt_kd <= 0):
        raise ValueError("Measured WT Kd must be finite and positive, in M (e.g. 1e-8 for 10 nM)")
    if layout and (not np.isfinite(layout_timeout) or layout_timeout <= 0):
        raise ValueError("Layout timeout must be finite and positive seconds")
    group_a, group_b = chain_group(chains_a), chain_group(chains_b)
    if set(group_a) & set(group_b):
        raise ValueError("Partner chain groups overlap; each chain must belong to one partner")
    inferred_task = "mutation" if mutations or variants is not None else "absolute"
    model = resolve_model(task or inferred_task, model, cache_dir=cache_dir)
    if not model.is_file():
        raise ValueError(f"Model file does not exist: {model}. Prepare models with topoppi-affinity models download")
    from topoppi_affinity.production import load_model
    bundle = load_model(model)
    info = model_info(bundle)
    model_task = info.get("task", info["kind"])
    if task and task != "custom" and model_task != task:
        raise ValueError(f"Selected task {task} conflicts with the model ({info['kind']}); choose its task or Custom model")
    if mutation_context and info["kind"] in {"known", "unseen"} and mutation_context != info["kind"]:
        raise ValueError("Mutation context conflicts with the selected artifact")
    if wt_kd is not None and info["target"] == "pkd":
        raise ValueError("Measured WT Kd is only used with a mutation-effect model")
    if prior_csv and prior_repository:
        raise ValueError("Choose a prior CSV or an author repository, not both")
    for label, path in (("Prior CSV", prior_csv), ("Prior checkpoint", prior_checkpoint)):
        if path and not Path(path).expanduser().is_file():
            raise ValueError(f"{label} does not exist: {path}")
    if prior_repository and not Path(prior_repository).expanduser().is_dir():
        raise ValueError(f"Prior repository does not exist: {prior_repository}")
    if prepared is None:
        loader = PDBLoader(str(structure))
        coords_a, atoms_a = load_partner_data(loader, group_a)
        coords_b, atoms_b = load_partner_data(loader, group_b)
    else:
        left, right = prepared.partner(group_a), prepared.partner(group_b)
        coords_a, atoms_a, coords_b, atoms_b = left.coords, left.atoms, right.coords, right.atoms
    residues = {(r["chain"], r["position"], r["insertion"]): r["wt"] for r in _structure_residues(atoms_a + atoms_b)}
    requests = variants_table(mutations, variants)
    errors = {}
    for index, row in requests.iterrows():
        try:
            if info["target"] == "pkd" and row.mutations:
                raise ValueError("Absolute affinity predicts this structure; select a mutation-effect task for variants")
            if info["target"] == "ddg" and not row.mutations:
                raise ValueError("A mutation-effect task requires at least one substitution")
            seen = set()
            for mutation in parse_mutations(row.mutations):
                key = (mutation["chain"], mutation["position"], mutation["insertion"])
                if key in seen:
                    raise ValueError(f"Duplicate mutation site: {key}; put alternative substitutions in separate variants")
                seen.add(key)
                if key not in residues:
                    raise ValueError(f"Mutation residue absent from supplied structure/selected chains: {key}")
                if residues[key] != mutation["wt"]:
                    raise ValueError(f"WT identity mismatch for {key}: structure={residues[key]}, mutation={mutation['wt']}")
        except (ValueError, KeyError) as exc:
            errors[index] = str(exc)
    if len(errors) == len(requests):
        details = "; ".join(f"{requests.loc[i, 'variant_id']}: {message}" for i, message in errors.items())
        raise ValueError(f"No valid input variants. {details}")
    if prior_csv:
        supplied = pd.read_csv(Path(prior_csv).expanduser(), keep_default_na=False, dtype={"sample_id": str})
        if not {"sample_id", "predicted", "target", "method", "status", "error"}.issubset(supplied.columns):
            raise ValueError("Prior CSV requires sample_id,predicted,target,method,status,error columns")
        valid_ids = set(requests.loc[~requests.index.isin(errors), "sample_id"])
        if supplied.sample_id.duplicated().any() or set(supplied.sample_id) != valid_ids:
            raise ValueError("Prior CSV sample IDs must match the valid input variants exactly")
        if info["kind"] not in {"known", "unseen"}:
            raise ValueError("A prior CSV applies only to a known/unseen mutation ensemble")
        required = (info["training"]["residual"]["prior_method"] if info["kind"] == "known"
                    else info["training"]["prior_method"])
        if not supplied.method.eq(required).all() or not supplied.target.eq("ddg").all():
            raise ValueError(f"The selected model requires {required} ddG priors in kcal/mol")
    neural = not (info["kind"] in {"descriptor", "production"} or (info["kind"] == "fitted" and info["training"].get("representation") == "conventional_descriptors"))
    if neural:
        missing = [name for name in ("torch", "esm") if importlib.util.find_spec(name) is None]
        if missing:
            raise ImportError(f"Missing model dependencies: {', '.join(missing)}. Run topoppi-install-affinity")
    resolved_device = resolve_device(device, neural=neural)
    previous = Path(output).expanduser().resolve() if output else None
    if previous and (previous / "analysis.json").exists() and not retry:
        raise ValueError("This output already contains an analysis; choose a new directory or use --retry")
    if previous and previous.exists() and not previous.is_dir():
        raise ValueError(f"Output must be a directory: {previous}")
    if previous and not retry:
        result_names = ("requests.csv", "predictions.csv", "residues.csv", "inputs.csv", "interface.npz",
                        "features", "prior", "prior.csv", "atlas_0.npz", "atlas_1.npz")
        occupied = [name for name in result_names if (previous / name).exists()]
        if occupied:
            raise ValueError(f"Output already contains result files: {', '.join(occupied)}; "
                             "choose a new directory or use --retry")
    destination = (_unused_directory(previous.with_name(previous.name + "-retry")) if retry and previous and previous.exists()
                   else previous if previous else _unused_directory(Path.cwd() / "topoppi-results" /
                       f"{structure.stem}-{task or info['kind']}-{datetime.now().strftime('%Y%m%d-%H%M%S')}"))
    parent = destination
    while not parent.exists():
        parent = parent.parent
    if not os.access(parent, os.W_OK):
        raise ValueError(f"Output directory is not writable: {parent}")
    return {"structure": structure, "model": model, "model_bundle": bundle, "model_info": info,
            "requests": requests, "row_errors": errors, "output": destination, "previous": previous if retry else None,
            "device": resolved_device, "chain_data": (coords_a, atoms_a, coords_b, atoms_b), "neural": neural}


def _residue_table(bundle, prepared):
    records = []
    side_for_chain = {chain: side for side, p in enumerate(prepared.partners) for chain in p.chains}
    for key, r in bundle["residues"].items():
        records.append({"analysis_id": prepared.analysis_id, "residue_key": key, "side": side_for_chain[r["chain"]],
            "chain": r["chain"], "residue_number": r["residue_number"], "insertion": r["insertion"],
            "wt": r["aa"], "contact_degree": r["contact_degree"], "contact_atom_pairs": r["contact_atom_pairs"],
            "partner_distance_A": r["partner_distance"], "surface_area_A2": r["surface_area"],
            "footprint_area_A2": r["footprint_area"], "footprint_patch_count": r["footprint_patch_count"],
            "coverage": "footprint" if r["footprint_area"] > 0 else "outside_retained_interface",
            "source": "input_descriptor_3D_not_model_attribution"})
    return pd.DataFrame(records, columns=RESIDUE_COLUMNS)


def _mapped_sites(mutations, residues):
    lookup = {(r.chain, int(r.residue_number), r.insertion): r for r in residues.itertuples()}
    sites = []
    for mutation in parse_mutations(mutations):
        r = lookup.get((mutation["chain"], int(mutation["position"]), mutation["insertion"]))
        sites.append({**mutation, "residue_key": None if r is None else r.residue_key,
            "side": None if r is None else int(r.side),
            "coverage": "absent_from_structure" if r is None else r.coverage})
    return sites


def analyze(structure, chains_a, chains_b, model=None, output=None, *, mutations=(), variants=None,
            source="user", complex_id=None, temperature=298.15, wt_kd=None, sequence_dir=None,
            device="auto", prior_csv=None, prior_repository=None, prior_checkpoint=None,
            mutation_context=None, scope_note="", layout=False, optcuts_bin="OptCuts_bin",
            layout_timeout=600., log=print, cancel_event=None, event_callback=None, retry=False, task=None, cache_dir=None, prepared=None):
    """Run one record; save available results even if either branch fails."""
    from topoppi_affinity import __version__ as plugin_version
    from topoppi_affinity.inference import build_features, compute_prior, predict_features
    from topoppi_affinity.resources import cache_root, prior_options

    def notify(stage, status="running", message="", **extra):
        if message:
            log(message)
        if event_callback:
            event_callback({"stage": stage, "status": status, "message": message, **extra})

    check_cancel(cancel_event)
    notify("validation", message="Checking structure, variants, model and runtime before preparing surfaces")
    checked = validate_inputs(structure, chains_a, chains_b, model, output, mutations=mutations, variants=variants,
        temperature=temperature, wt_kd=wt_kd, device=device, task=task, prior_csv=prior_csv,
        prior_repository=prior_repository, prior_checkpoint=prior_checkpoint, mutation_context=mutation_context,
        layout=layout, layout_timeout=layout_timeout, retry=retry, cache_dir=cache_dir, prepared=prepared)
    structure, model, output, device = (checked[key] for key in ("structure", "model", "output", "device"))
    model_bundle, info = checked["model_bundle"], checked["model_info"]
    if prior_csv is None and prior_repository is None:
        options = prior_options(info["kind"], cache_dir=cache_dir)
        prior_repository = options.get("repository")
        prior_checkpoint = prior_checkpoint or options.get("checkpoint")
    # Registries belong to an analysis; only sequence-addressed embeddings are shared.
    previous = checked["previous"]
    if sequence_dir is None:
        sequence_dir = output / "features" / "sequences"
    shared_embeddings = cache_root(cache_dir) / "embeddings" / "esmc_600m"
    output.mkdir(parents=True, exist_ok=True)
    requests = checked["requests"]
    requests.to_csv(output / "requests.csv", index=False)
    pd.DataFrame(columns=RESIDUE_COLUMNS).to_csv(output / "residues.csv", index=False)
    report = {"format": "topoppi-affinity-analysis", "version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(), "core_version": core_version, "plugin_version": plugin_version,
        "output_directory": str(output), "device": device, "model": info,
        "software": _software_versions(),
        "positioning": "Affinity modelling using interface descriptors produced by the TopoPPI workflow",
        "model_path": str(Path(model).expanduser().resolve()), "input": {"path": str(Path(structure).resolve()),
            "chains_a": chains_a, "chains_b": chains_b, "source": source, "complex_id": complex_id,
            "task": task, "temperature": temperature, "wt_kd": wt_kd},
        "prediction_status": "pending", "layout_status": "pending" if layout else "not_requested",
        "interface_status": "pending", "descriptor_status": "pending",
        "files": {"requests": "requests.csv", "predictions": "predictions.csv", "residues": "residues.csv"},
        "layout": {}, "semantics": {"pkd": "-log10(Kd/M); whole-complex prediction",
            "ddg": "Gmut-Gwt in kcal/mol; positive weakens binding; joint variants retain their total effect",
            "display": "Contacts and areas are input descriptors, not model attribution or residue energies",
            "uncertainty": "No per-prediction confidence interval is supplied"}}
    predictions = requests.copy()
    predictions["predicted"] = np.nan
    predictions["status"], predictions["error"], predictions["sites"] = "error", "Prediction not completed", "[]"
    for index, error in checked["row_errors"].items():
        predictions.at[index, "error"] = error
    if previous:
        report["retry_of"] = str(previous)
    report["validation_errors"] = {requests.loc[i, "variant_id"]: error for i, error in checked["row_errors"].items()}
    predictions.to_csv(output / "predictions.csv", index=False)
    _write_report(report, output)
    notify("interface", message=f"Using {device}; results will be saved in {output}")
    try:
        if prepared is None and previous and (previous / "interface.npz").is_file():
            old = load_interface(previous / "interface.npz")
            stat = structure.stat()
            if (old.source["path"] == str(structure) and old.source["mtime_ns"] == stat.st_mtime_ns
                    and old.source["size"] == stat.st_size
                    and [p.chains for p in old.partners] == [chain_group(chains_a), chain_group(chains_b)]):
                prepared = old
                prepared.analysis_id = uuid.uuid4().hex
                prepared.source.update(source=source, complex_id=complex_id or structure.stem)
                report["reused_interface"] = str(previous / "interface.npz")
                notify("interface", message="Reusing the unchanged native interface from the previous analysis")
        check_cancel(cancel_event)
        if prepared is None:
            prepared = prepare_interface(structure, chains_a, chains_b, source=source, complex_id=complex_id,
                                         chain_data=checked["chain_data"], log=log, cancel_event=cancel_event)
        report["analysis_id"] = prepared.analysis_id
        report["interface"] = prepared.metadata()
        save_interface(prepared, output / "interface.npz")
        report["interface_status"] = "ok"
        report["files"]["interface"] = "interface.npz"
    except Exception as exc:
        status = "cancelled" if isinstance(exc, InterruptedError) else "error"
        report.update(interface_status=status, descriptor_status="unavailable", prediction_status=status, layout_status="unavailable",
                      interface_error=f"{type(exc).__name__}: {exc}")
        unfinished = predictions.error.eq("Prediction not completed")
        predictions.loc[unfinished, "error"] = report["interface_error"]
        predictions.loc[unfinished, "status"] = status
        predictions.to_csv(output / "predictions.csv", index=False)
        _write_report(report, output)
        notify("finished", status, report["interface_error"], analysis_json=str(output / "analysis.json"))
        return report

    try:
        bundle = extract_complex_features(structure, chains_a, chains_b, prepared=prepared)
        residues = _residue_table(bundle, prepared)
        residues.to_csv(output / "residues.csv", index=False)
        report["descriptor_status"] = "ok"
        target = info["target"]
        if info["kind"] in {"known", "unseen"}:
            if mutation_context is not None and mutation_context != info["kind"]:
                raise ValueError("Mutation context conflicts with the selected artifact")
            training = info["training"].get("residual", info["training"].get("base", {}))
            members = set(training.get("train_complex_ids", []))
            report["scope_assignment"] = {"context": info["kind"], "selection": "explicit artifact selection",
                "exact_complex_member": prepared.source["complex_id"] in members if members else None,
                "membership_note": "Exact identifier check only; sequence-family novelty is not established",
                "user_basis": scope_note or "Not supplied; training membership/family applicability remains unverified"}
        valid = []
        for i, row in requests.iterrows():
            try:
                sites = _mapped_sites(row.mutations, residues)
                predictions.at[i, "sites"] = json.dumps(sites)
                if i in checked["row_errors"]:
                    raise ValueError(checked["row_errors"][i])
                if target == "pkd" and row.mutations:
                    raise ValueError("The absolute model predicts the supplied structure; use a ddG model for variants")
                if target == "ddg" and not row.mutations:
                    raise ValueError("A ddG model requires an explicit substitution or joint variant")
                if target == "ddg":
                    mutation_features(bundle, row.mutations)
                valid.append(i)
                predictions.at[i, "error"] = ""
            except (ValueError, KeyError) as exc:
                predictions.at[i, "error"] = str(exc)
        if not valid:
            raise ValueError("No valid input variants; see predictions.csv for row-specific errors")
        frame = requests.loc[valid].copy()
        frame["source"], frame["complex_id"] = source, prepared.source["complex_id"]
        frame["structure_path"] = prepared.source["path"]
        frame["chains_a"], frame["chains_b"] = [",".join(p.chains) for p in prepared.partners]
        frame["pdb_id"] = Path(structure).stem
        frame["n_mutations"] = [len(parse_mutations(m)) for m in frame.mutations]
        frame["temperature_K"] = temperature
        frame.to_csv(output / "inputs.csv", index=False)
        report["files"]["inputs"] = "inputs.csv"
        paths, protocol = build_features(prepared, bundle, frame, model_bundle, output / "features",
            sequence_dir=sequence_dir, device=device, log=log, cancel_event=cancel_event,
            shared_embeddings=shared_embeddings, event_callback=event_callback, cache_dir=cache_dir)
        report["feature_protocol"] = protocol
        report["files"]["features"] = {("descriptor_features" if key == "legacy_features" else key):
            str(value.relative_to(output)) for key, value in paths.items()}
        if cancel_event is not None and cancel_event.is_set():
            raise InterruptedError("Analysis cancelled")
        prior, provenance = compute_prior(frame, info, output, prior_csv=prior_csv,
            repository=prior_repository, checkpoint=prior_checkpoint, device=device,
            cancel_event=cancel_event, log=log)
        report["prior"] = provenance
        if prior is not None:
            report["files"]["prior"] = "prior.csv"
        log(f"Predicting with {info['scope']}")
        check_cancel(cancel_event)
        result = predict_features(model, **paths, prior=prior, temperature=temperature, wt_kd=wt_kd)
        aligned = result.set_index("sample_id")
        for column in aligned.columns:
            if column not in predictions:
                predictions[column] = None
            predictions.loc[valid, column] = frame.sample_id.map(aligned[column]).to_numpy()
        report["prediction_status"] = "ok" if len(valid) == len(requests) else "partial"
        report["successful_predictions"] = len(result)
        report["base_only_predictions"] = int(result.status.eq("base_only").sum())
        report["warnings"] = []
        if report["base_only_predictions"]:
            report["warnings"].append("Only the base model was used for some variants; the requested prior was unavailable or failed. See prior status and row errors.")
        if info["kind"] in {"known", "unseen"} and not scope_note:
            report["warnings"].append("Model applicability is unconfirmed: no basis for the known/unseen selection was supplied.")
    except Exception as exc:
        cancelled = isinstance(exc, InterruptedError)
        if report["descriptor_status"] == "pending":
            report["descriptor_status"] = "cancelled" if cancelled else "error"
        report["prediction_status"] = "cancelled" if cancelled else "error"
        report["prediction_error"] = f"{type(exc).__name__}: {exc}"
        if "out of memory" in str(exc).lower():
            report["prediction_error"] += (
                "; available device memory is insufficient. Retry with --device cpu "
                "(GUI Advanced: Device = cpu), or free GPU memory before retrying. The model protocol is unchanged.")
        if isinstance(exc, ImportError):
            report["prediction_error"] += (
                "; install this plugin's [models] dependencies in the active Python environment; "
                "author priors also require their documented source environment")
        failed = predictions.error.isin(["", "Prediction not completed"])
        if cancelled:
            failed &= predictions.predicted.isna()
            predictions.loc[failed, "status"] = "cancelled"
        predictions.loc[failed, "error"] = report["prediction_error"]
        log(report["prediction_error"])
    predictions["analysis_id"] = prepared.analysis_id
    predictions.to_csv(output / "predictions.csv", index=False)
    # Prediction is durable before any optional solver or plotting is attempted.
    _write_report(report, output)
    notify("prediction_ready", report["prediction_status"], "Prediction results are ready; optional layout runs separately",
           analysis_json=str(output / "analysis.json"))

    if layout and not (cancel_event is not None and cancel_event.is_set()):
        _generate_layout(prepared, report, output, optcuts_bin, layout_timeout, log, cancel_event, event_callback)
    elif layout:
        report["layout_status"] = "cancelled"
    _write_report(report, output)
    notify("finished", report["prediction_status"], "Analysis saved", analysis_json=str(output / "analysis.json"))
    return report


def _generate_layout(prepared, report, output, optcuts_bin, timeout, log, cancel_event, event_callback):
    from topoppi.pipeline import layout_interface
    from topoppi.visualization.atlas_io import save_atlas

    report["layout_status"] = report["layout_attempt_status"] = "running"
    previous_warning = report.pop("layout_attempt_warning", None)
    report["warnings"] = [warning for warning in report.get("warnings", []) if warning != previous_warning]
    attempts, failed_sides, retained_sides = [], [], []
    _write_report(report, output)
    for side in (0, 1):
        pending = output / f".atlas_{side}.pending.npz"
        try:
            check_cancel(cancel_event)
            log(f"Generating optional atlas for partner {','.join(prepared.partners[side].chains)}")
            patches, layout_report = layout_interface(prepared, side=side,
                optcuts=replace(DEFAULT_RUN_CONFIG.optcuts, optcuts_bin=optcuts_bin, timeout_sec=timeout),
                cancel_event=cancel_event)
            filename = f"atlas_{side}.npz"
            save_atlas(pending, patches, prepared.visualizer(side),
                style_config={"map_style": "footprints", "footprint_labels": "highlighted"},
                run_metadata={"analysis_id": prepared.analysis_id, "analysis_json": "analysis.json", "side": side,
                    "interface": prepared.metadata(), "layout_report": layout_report})
            check_cancel(cancel_event)
            pending.replace(output / filename)
            report["layout"][str(side)] = {"status": "ok", "atlas": filename}
            status = "ok"
        except Exception as exc:
            status = "cancelled" if isinstance(exc, InterruptedError) or (cancel_event and cancel_event.is_set()) else "error"
            # Keep an existing valid atlas when an explicit regeneration fails.
            if report["layout"].get(str(side), {}).get("status") == "ok":
                report["layout"][str(side)]["last_attempt_error"] = f"{type(exc).__name__}: {exc}"
                retained_sides.append("AB"[side])
            else:
                report["layout"][str(side)] = {"status": status, "error": f"{type(exc).__name__}: {exc}"}
            failed_sides.append("AB"[side])
            log(f"Atlas attempt for partner {'AB'[side]} {status}: {exc}")
        finally:
            pending.unlink(missing_ok=True)
        attempts.append(status)
        _write_report(report, output)
        if event_callback:
            event_callback({"stage": "layout_ready", "status": status,
                            "side": side, "analysis_json": str(output / "analysis.json"), "message": "Atlas status updated"})
    statuses = [item["status"] for item in report["layout"].values()]
    report["layout_status"] = ("ok" if all(s == "ok" for s in statuses) else "partial" if "ok" in statuses
                               else "cancelled" if "cancelled" in statuses else "error")
    report["layout_attempt_status"] = ("cancelled" if "cancelled" in attempts else "error" if "error" in attempts else "ok")
    if failed_sides:
        outcome = "cancelled" if report["layout_attempt_status"] == "cancelled" else "failed"
        warning = f"Latest atlas generation {outcome} for partner(s) {', '.join(failed_sides)}."
        if retained_sides:
            warning += f" Earlier atlases for {', '.join(retained_sides)} remain viewable; they do not reflect this attempt's parameters."
        report["layout_attempt_warning"] = warning
        report["warnings"].append(warning)
        log(warning)
    _write_report(report, output)


def layout_existing(path, *, optcuts_bin="OptCuts_bin", layout_timeout=600., cancel_event=None,
                    log=print, event_callback=None):
    """Generate/retry atlases from a saved interface without recalculating predictions."""
    if not np.isfinite(layout_timeout) or layout_timeout <= 0:
        raise ValueError("Layout timeout must be finite positive seconds")
    document = load_analysis(path)
    if document.interface is None:
        raise ValueError("This analysis has no saved native interface; retry the analysis first")
    _generate_layout(document.interface, document.report, document.path.parent, optcuts_bin, layout_timeout,
                     log, cancel_event, event_callback)
    return document.report


@dataclass
class AnalysisDocument:
    path: Path
    report: dict
    interface: object
    predictions: pd.DataFrame
    residues: pd.DataFrame

    def selection(self, variant_id, side=0):
        matches = self.predictions[self.predictions.variant_id == variant_id]
        if matches.empty:
            available = ", ".join(self.predictions.variant_id.astype(str).head(10))
            raise ValueError(f"Unknown variant ID {variant_id!r}. Available IDs include: {available}")
        row = matches.iloc[0]
        sites = json.loads(row.sites)
        visible = [s["residue_key"] for s in sites if s["side"] == side and s["coverage"] == "footprint"]
        hidden = [s for s in sites if s["residue_key"] not in visible]
        return row, visible, hidden

    def view(self, side=0, mode="surface", variant_id=None, descriptor=None):
        style = {"map_style": "footprints", "view": mode, "footprint_labels": "highlighted",
                 "highlight_residues": (), "show_seams": mode == "atlas"}
        if self.interface is None:
            raise ValueError("This analysis has no prepared interface")
        if mode == "atlas":
            from topoppi.visualization.atlas_io import load_atlas

            item = self.report["layout"].get(str(side), {})
            if item.get("status") != "ok":
                raise ValueError(item.get("error", "No atlas available for this side; use the 3D view"))
            atlas = load_atlas(self.path.parent / item["atlas"])
            if (atlas.metadata.get("analysis_id") != self.report["analysis_id"] or
                    atlas.metadata.get("side") != side or
                    atlas.metadata.get("interface", {}).get("groups") != self.report["interface"]["groups"]):
                raise ValueError("Atlas belongs to a different analysis or partner group")
            patches, viz = atlas.patches, atlas.visualizer
        else:
            patches, viz = self.interface.partners[side].patches, self.interface.visualizer(side)
        if variant_id:
            _, style["highlight_residues"], _ = self.selection(variant_id, side)
        if descriptor:
            if descriptor not in {"contact_degree", "footprint_area_A2", "surface_area_A2"}:
                raise ValueError("Only input descriptors can be projected onto residues")
            selected = self.residues[self.residues.side == side]
            style["annotation_values"] = dict(zip(selected.residue_key, selected[descriptor], strict=True))
            style["annotation_label"] = f"Input descriptor: {descriptor}"
        return patches, viz, style


def load_analysis(path):
    """Reopen portable results without the input structure, model or embedding stack."""
    from topoppi_affinity.model_metadata import public_model_info

    path = Path(path).expanduser().resolve()
    if path.is_dir():
        path /= "analysis.json"
    report = json.loads(path.read_text())
    if report.get("format") != "topoppi-affinity-analysis" or report.get("version") != 1:
        raise ValueError("Unsupported TopoPPI Affinity analysis")
    if report.get("model"):
        report["model"] = public_model_info(report["model"])
    native = report["files"].get("interface")
    interface = load_interface(path.parent / native) if native else None
    if interface is not None and interface.analysis_id != report["analysis_id"]:
        raise ValueError("Native interface belongs to another analysis")
    predictions = pd.read_csv(path.parent / report["files"]["predictions"], keep_default_na=False,
        dtype={"variant_id": str, "sample_id": str, "mutations": str})
    required = {"variant_id", "mutations", "predicted", "status", "error", "sites"}
    missing = required - set(predictions.columns)
    if missing:
        raise ValueError(f"Saved predictions table is missing columns: {', '.join(sorted(missing))}")
    if predictions.empty and report["prediction_status"] in {"ok", "partial"}:
        raise ValueError("Saved predictions table has no rows; restore the complete file")
    residues_path = report["files"].get("residues")
    residues = pd.read_csv(path.parent / residues_path, keep_default_na=False,
        dtype={"chain": str, "insertion": str, "residue_key": str}) if residues_path else pd.DataFrame()
    return AnalysisDocument(path, report, interface, predictions, residues)
