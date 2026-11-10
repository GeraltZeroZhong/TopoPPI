"""Standalone frozen-model ranking of a user-selected WT complex partner.

Feature context is the complete WT interface, including native alanines and
unmeasured residues. Only standard non-alanine residues are ranked candidates.
The default model shares mutation descriptors with the Affinity plugin.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from topoppi.config import SurfaceConfig, TopologyConfig
from topoppi.interface import chain_group, prepare_interface

from .block_ranking import BlockPairwiseRanker
from .features import extract_features
from .graph_operators import build_graph_operators
from .graph_ranking import GraphRanker
from .model_io import predict_frozen
from .partner_features import KEYS, prepare_partner_features
from .rank_fusion import aggregate_ranks, partner_rank_percentiles, standard_nonalanine_candidates
from .representations import (
    all_residue_embeddings,
    center_within_queries,
    relative_new_columns,
)
from .resources import resolve_model_dir


@dataclass
class _ModelBundle:
    directory: Path
    metadata: dict
    feature_protocol: dict
    components: list


def _feature_protocol(value=None, *, mutation_surface=False):
    value = {} if value is None else value
    if not isinstance(value, dict):
        raise ValueError("Model feature_protocol must be a JSON object.")
    try:
        surface = SurfaceConfig(**value.get("surface", {}))
        topology = TopologyConfig(**value.get("topology", {}))
    except TypeError as exc:
        raise ValueError(f"Invalid model feature_protocol: {exc}") from exc
    surface.validate()
    topology.validate()
    contact = value.get("contact_cutoff_angstrom", 6.)
    if mutation_surface:
        radii = value.get("spatial_query_radii_angstrom", [4.0, 12.0])
        if contact != 6. or list(radii) != [4.0, 12.0]:
            raise ValueError("The mutation-derived model uses a 6 A interface and spatial query radii of 4 and 12 A.")
        return {"surface": asdict(surface), "topology": asdict(topology),
                "contact_cutoff_angstrom": contact, "spatial_query_radii_angstrom": list(radii)}
    spatial = value.get("spatial_neighbor_cutoff_angstrom", 10.)
    if contact != 6. or spatial != 10.:
        raise ValueError("Frozen inference currently supports the 6 A interface and 10 A spatial-neighbor protocol.")
    return {"surface": asdict(surface), "topology": asdict(topology),
            "contact_cutoff_angstrom": contact, "spatial_neighbor_cutoff_angstrom": spatial}


def _load_bundle(model_dir, inherited_protocol=None):
    directory = resolve_model_dir(model_dir)
    metadata = json.loads((directory / "model.json").read_text())
    if not isinstance(metadata, dict):
        raise ValueError(f"Model metadata in {directory / 'model.json'} must be a JSON object.")
    kind = metadata.get("kind")
    if kind is None:
        raise ValueError(f"Model metadata in {directory / 'model.json'} must declare a model kind.")
    protocol = _feature_protocol(metadata.get("feature_protocol", inherited_protocol),
                                 mutation_surface=kind == "mutation_surface")
    if kind not in {"linear", "mlp", "xgb", "graph", "block_linear", "fusion", "mutation_surface"}:
        raise ValueError(f"Unsupported frozen model kind: {kind}")
    if kind != "fusion" and "columns" not in metadata:
        raise ValueError(f"Model metadata in {directory / 'model.json'} must declare its feature columns.")
    components = []
    if kind == "fusion":
        for component in metadata["components"]:
            child = _load_bundle(directory / component["path"], protocol)
            if child.feature_protocol != protocol:
                raise ValueError("Fusion components require different structure feature protocols.")
            components.append((component["name"], child))
        if not components:
            raise ValueError("A frozen fusion must contain at least one component.")
    return _ModelBundle(directory, metadata, protocol, components)


def _leaves(bundle):
    if bundle.components:
        return [leaf for _, child in bundle.components for leaf in _leaves(child)]
    return [bundle]


def model_requirements(model_dir=None):
    """Read recursive column/embedding/operator requirements without scoring."""
    bundle = _load_bundle(model_dir)
    leaves = _leaves(bundle)
    return {
        "columns": sorted({column for leaf in leaves for column in leaf.metadata["columns"]}),
        "needs_embeddings": any(leaf.metadata["kind"] == "block_linear"
                                or (leaf.directory / "embedding_pca.npz").exists() for leaf in leaves),
        "needs_graph": any(leaf.metadata["kind"] == "graph"
                           or leaf.metadata.get("incremental_control") == "native_source_shuffle"
                           for leaf in leaves),
        "needs_uv": any(leaf.metadata["kind"] == "mutation_surface" for leaf in leaves),
        "feature_protocol": bundle.feature_protocol,
    }


def _score_bundle(bundle, frame, embeddings, operators, graphs):
    metadata, directory = bundle.metadata, bundle.directory
    if metadata["kind"] == "mutation_surface":
        from .mutation_model import MutationSurfaceRegressor

        model = MutationSurfaceRegressor(directory / "model.json")
        return model.predict(frame[list(model.feature_names)].to_numpy(float), model.feature_names), {}
    if metadata["kind"] == "fusion":
        eligible = standard_nonalanine_candidates(frame)
        queries = frame.query_id.to_numpy()
        scores, contributions = [], {}
        for name, child in bundle.components:
            current, _ = _score_bundle(child, frame, embeddings, operators, graphs)
            scores.append(current)
            contributions[f"component_rank.{name}"] = partner_rank_percentiles(current, queries, eligible)
        return aggregate_ranks(scores, queries, eligible), contributions
    if metadata["kind"] == "linear" and metadata.get("incremental_control", "none") != "none":
        from .feature_controls import control_features

        frame = control_features(frame, graphs, metadata["incremental_control"], metadata["incremental_control_seed"])
    columns = metadata["columns"]
    X = frame[columns].to_numpy(float)
    if not np.isfinite(X).all():
        raise ValueError(f"Required features are nonfinite for frozen model {directory.name}.")
    if metadata["kind"] == "block_linear":
        if embeddings is None:
            raise ValueError("This frozen model requires full WT residue embeddings.")
        centered = center_within_queries(embeddings, frame.query_id.to_numpy())
        model = BlockPairwiseRanker.load(directory / metadata.get("model_directory", "estimator"))
        return model.decision_function(X, centered), {}
    if metadata["kind"] == "graph":
        if (directory / "embedding_pca.npz").exists():
            if embeddings is None:
                raise ValueError("This frozen graph model requires WT residue embeddings.")
            centered = center_within_queries(embeddings, frame.query_id.to_numpy())
            with np.load(directory / "embedding_pca.npz", allow_pickle=False) as pca:
                X = np.column_stack([X, (centered - pca["mean"]) @ pca["components"].T])
        local, context = operators[metadata["operator_kind"]], operators["context"]
        scores = [GraphRanker.load(directory / f"seed_{seed}").decision_function(X, local, context)
                  for seed in metadata["seeds"]]
        return np.mean(scores, axis=0), {}
    return predict_frozen(directory, frame, embeddings), {}


def score_model(model_dir, frame, graphs=None, embeddings=None):
    """Score an already prepared full-interface frame without labels or fitting.

    ``graphs`` is a complex-ID dictionary of the raw full-residue graph arrays
    returned by extract_features. ``embeddings`` contains raw WT embeddings in
    frame row order; query centering happens here. Returns (scores, components),
    where fusion component entries are full-population midrank percentiles
    named ``component_rank.<name>``. Fusion non-candidates have NaN scores.
    """
    bundle = _load_bundle(model_dir)
    leaves = _leaves(bundle)
    needs_operators = any(leaf.metadata["kind"] == "graph" for leaf in leaves)
    needs_graph = needs_operators or any(leaf.metadata.get("incremental_control") == "native_source_shuffle"
                                        for leaf in leaves)
    if embeddings is not None:
        embeddings = np.asarray(embeddings, dtype=float)
        if embeddings.ndim != 2 or len(embeddings) != len(frame) or not np.isfinite(embeddings).all():
            raise ValueError("WT embeddings must be finite and aligned to the full interface frame.")
    if needs_graph and graphs is None:
        raise ValueError("This frozen model requires the WT residue graphs.")
    operators = build_graph_operators(frame, graphs)[0] if needs_operators else None
    return _score_bundle(bundle, frame, embeddings, operators, graphs)


def _join_auxiliary(frame, auxiliary, prefixes):
    columns = [column for column in auxiliary if column.startswith(prefixes)]
    return frame.merge(auxiliary[KEYS + columns], on=KEYS, how="left", validate="one_to_one")


def _prepare_features(structure, groups, complex_id, prepared, requirements):
    """Run only auxiliary extractors needed by the frozen predictor columns."""
    all_features, graph = extract_features(structure, *groups, prepared=prepared)
    all_features["complex_id"] = complex_id
    graphs = {complex_id: graph}
    columns = requirements["columns"]
    relative_sources = sorted({column.split(".", 2)[2] for column in columns
                               if column.startswith(("new_relative.z.", "new_relative.pct."))})
    raw = set(columns) | set(relative_sources)
    reports = {}
    if any(column.startswith("contact_context.") for column in raw):
        from .paired_surface_features import cross_partner_context

        all_features, _ = cross_partner_context(all_features, graphs)
        reports["contact_context"] = {"weights": "inter-partner atom contact counts", "labels_used": False}
    frame = prepare_partner_features(all_features)
    if not len(frame):
        raise ValueError("The selected complete partners have no WT residues within the 6 A interface cutoff.")
    for prefix, module_name, function_name in (
        (("surface_pair.",), "paired_surface_features", "extract_paired_surface"),
        (("mesh_diffusion.", "mesh_wavelet.", "mesh_coverage."), "mesh_diffusion", "extract_mesh_diffusion"),
        (("surface_rim.",), "surface_rim", "extract_surface_rim"),
    ):
        if any(column.startswith(prefix) for column in raw):
            from importlib import import_module

            function = getattr(import_module(f".{module_name}", __package__), function_name)
            auxiliary, report = function(structure, *groups, prepared=prepared)
            auxiliary["complex_id"] = complex_id
            frame = _join_auxiliary(frame, auxiliary, prefix)
            reports[module_name] = report
    graph_kinds = sorted({column.split(".")[1] for column in raw
                          if column.startswith(("diffusion.", "scattering."))})
    if graph_kinds:
        from .footprint_diffusion import add_graph_diffusion

        auxiliary, _, report = add_graph_diffusion(all_features, graphs, kinds=tuple(graph_kinds))
        frame = _join_auxiliary(frame, auxiliary, ("diffusion.", "scattering."))
        reports["footprint_diffusion"] = report
    if any(column.startswith("sasa.") for column in raw):
        from .solvent_accessibility import residue_sasa

        auxiliary = residue_sasa(structure, *groups, prepared=prepared)
        auxiliary["complex_id"] = complex_id
        frame = _join_auxiliary(frame, auxiliary, ("sasa.",))
        reports["sasa"] = {"method": "Shrake-Rupley", "sphere_points": 960, "probe_radius_angstrom": 1.4}
    if any(column.startswith(("typed_contact.", "spatial_baseline.")) for column in raw):
        from .physical_features import SCHEMA, add_physical_context

        frame = add_physical_context(frame, graphs)
        reports["physical_context"] = {"schema": SCHEMA, "labels_used": False}
    if relative_sources:
        frame, _ = relative_new_columns(frame, relative_sources)
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"Frozen model requires unsupported feature columns: {missing}")
    if not np.isfinite(frame[columns].to_numpy(float)).all():
        raise ValueError("Required frozen-model features do not cover every WT interface node.")
    return all_features, frame, graphs, reports


def _registry_embeddings(registry, frame, all_features):
    """Match full observed WT chains before reusing a supplied embedding registry."""
    registry = Path(registry)
    chains = pd.read_csv(registry / "chain_sequences.csv", keep_default_na=False)
    sequences = pd.read_csv(registry / "sequences.csv", keep_default_na=False).set_index("sequence_id").sequence
    observed = {}
    for chain, part in all_features.loc[all_features.ca_resolved].groupby("chain", sort=False):
        observed[str(chain)] = (part.residue_key.tolist(), "".join(part.wt))
    matches = []
    for complex_id, part in chains.groupby("complex_id", sort=True):
        selected = part.loc[part.chain.isin(observed)]
        if set(selected.chain) != set(observed) or selected.chain.duplicated().any():
            continue
        if all((json.loads(row.residue_keys), sequences[row.sequence_id]) == observed[str(row.chain)]
               for row in selected.itertuples(index=False)):
            matches.append(complex_id)
    if not matches:
        raise ValueError("Embedding registry does not match the complete observed WT chain sequences and residue keys.")
    requested = frame.complex_id.iloc[0]
    matched = requested if requested in matches else matches[0]
    aligned = frame.copy()
    aligned["complex_id"] = matched
    return all_residue_embeddings(registry, aligned), {"registry": str(registry.resolve()), "registry_complex_id": matched}


def rank_complex(structure, chains_a, chains_b, partner, model_dir=None, output_dir=None, *,
                 embedding_registry=None, esmc_checkpoint=None, device="cpu", optcuts_bin=None):
    """Rank one complete partner from a supplied WT structure and frozen model."""
    model_dir = resolve_model_dir(model_dir)
    if output_dir is None:
        raise ValueError("Choose an output_dir for the ranking results.")
    structure = Path(structure).expanduser().resolve()
    groups = tuple(sorted((chain_group(chains_a), chain_group(chains_b))))
    selected = chain_group(partner)
    if set(groups[0]) & set(groups[1]):
        raise ValueError("The two complete partner groups must not overlap.")
    if selected not in groups:
        raise ValueError("--partner must exactly match one complete --chains-a or --chains-b group.")
    if device not in {"cpu", "cuda"}:
        raise ValueError("device must be cpu or cuda.")
    if embedding_registry is not None and esmc_checkpoint is not None:
        raise ValueError("Supply either an existing embedding registry or a checkpoint for fresh embeddings.")
    bundle = _load_bundle(model_dir)
    requirements = model_requirements(model_dir)
    if requirements["needs_embeddings"] and embedding_registry is None and esmc_checkpoint is None:
        raise ValueError("This model requires --embedding-registry or an explicit --esmc-checkpoint.")
    protocol = requirements["feature_protocol"]
    complex_id = f"{structure.stem}:{','.join(groups[0])}:{','.join(groups[1])}"
    prepared = prepare_interface(structure, *groups, surface=SurfaceConfig(**protocol["surface"]),
                                 topology=TopologyConfig(**protocol["topology"]), complex_id=complex_id)
    return _rank_prepared(prepared, structure, groups, selected, model_dir, output_dir,
                          bundle, requirements, complex_id, embedding_registry, esmc_checkpoint, device,
                          optcuts_bin=optcuts_bin)


def rank_interface(prepared, partner=None, model_dir=None, output_dir=None, *,
                   embedding_registry=None, esmc_checkpoint=None, device=None, optcuts_bin=None, context=None):
    """Rank from shared Core geometry without reopening the source structure.

    The supplied geometry must match the frozen model's surface/topology protocol
    and contain both partner surfaces. Task features, sequence embeddings and
    ranking remain plugin-owned; no model is trained by this function.
    """
    from topoppi import GeometryRepresentation

    if not isinstance(prepared, GeometryRepresentation):
        raise TypeError("hotspot.rank requires a Core GeometryRepresentation; call topoppi.compute_geometry first.")
    values = dict(context.config_for("hotspot")) if context is not None else {}
    partner = partner if partner is not None else values.get("partner")
    model_dir = model_dir if model_dir is not None else values.get("model_dir")
    output_dir = output_dir if output_dir is not None else values.get(
        "output_dir", context.output_dir if context is not None else None)
    embedding_registry = embedding_registry if embedding_registry is not None else values.get("embedding_registry")
    esmc_checkpoint = esmc_checkpoint if esmc_checkpoint is not None else values.get("esmc_checkpoint")
    device = device if device is not None else values.get("device", "cpu")
    optcuts_bin = optcuts_bin if optcuts_bin is not None else values.get("optcuts_bin")
    model_dir = resolve_model_dir(model_dir)
    if partner is None or output_dir is None:
        raise ValueError("Hotspot ranking requires partner and output_dir arguments or plugin configuration.")
    groups = tuple(sorted(partner.chains for partner in prepared.partners))
    selected = chain_group(partner)
    if selected not in groups:
        raise ValueError("--partner must exactly match one complete --chains-a or --chains-b group.")
    if device not in {"cpu", "cuda"}:
        raise ValueError("device must be cpu or cuda.")
    if embedding_registry is not None and esmc_checkpoint is not None:
        raise ValueError("Supply either an existing embedding registry or a checkpoint for fresh embeddings.")
    if context is not None:
        context.check_cancelled()
    bundle = _load_bundle(model_dir)
    requirements = model_requirements(model_dir)
    if requirements["needs_embeddings"] and embedding_registry is None and esmc_checkpoint is None:
        raise ValueError("This model requires --embedding-registry or an explicit --esmc-checkpoint.")
    protocol = requirements["feature_protocol"]
    if asdict(prepared.surface_config) != protocol["surface"] or asdict(prepared.topology_config) != protocol["topology"]:
        raise ValueError("Supplied Core geometry does not match the frozen model's surface/topology protocol.")
    if any(item.mesh is None for item in prepared.partners):
        raise ValueError("Hotspot ranking requires Core geometry for both partner sides.")
    structure = Path(prepared.source["path"])
    complex_id = f"{structure.stem}:{','.join(groups[0])}:{','.join(groups[1])}"
    return _rank_prepared(prepared, structure, groups, selected, model_dir, output_dir,
                          bundle, requirements, complex_id, embedding_registry, esmc_checkpoint, device,
                          context=context, optcuts_bin=optcuts_bin)


def _rank_prepared(prepared, structure, groups, selected, model_dir, output_dir,
                   bundle, requirements, complex_id, embedding_registry, esmc_checkpoint, device, *,
                   context=None, optcuts_bin=None):
    if bundle.metadata["kind"] == "mutation_surface":
        from .mutation_inference import rank_mutation_surface

        return rank_mutation_surface(prepared, selected, model_dir, output_dir,
                                     optcuts_bin=optcuts_bin, context=context)
    protocol = requirements["feature_protocol"]
    all_features, frame, graphs, reports = _prepare_features(structure, groups, complex_id, prepared, requirements)
    output = Path(output_dir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    embeddings, embedding_report = None, None
    if requirements["needs_embeddings"]:
        if embedding_registry is None:
            from .sequence_embeddings import embed_sequences, prepare_interface_embedding_registry

            embedding_registry = output / "embedding_registry"
            prepare_interface_embedding_registry(prepared, embedding_registry, complex_id=complex_id)
            embed_sequences(embedding_registry, esmc_checkpoint, device=device)
        embeddings, embedding_report = _registry_embeddings(embedding_registry, frame, all_features)
    if context is not None:
        context.check_cancelled()
    scores, components = score_model(model_dir, frame, graphs, embeddings)
    eligible = standard_nonalanine_candidates(frame)
    side = groups.index(selected)
    positions = np.flatnonzero((frame.side.to_numpy() == side) & eligible)
    if not np.isfinite(np.asarray(scores)[positions]).all():
        raise ValueError("Frozen model produced nonfinite candidate ranking scores.")
    result = frame.iloc[positions][KEYS + ["chain", "residue_number", "insertion_code", "side", "wt"]].copy()
    result["rank_score"] = np.asarray(scores)[positions]
    result["rank_min"] = result.rank_score.rank(method="min", ascending=False).astype(int)
    result["rank_max"] = result.rank_score.rank(method="max", ascending=False).astype(int)
    result["full_candidate_percentile"] = partner_rank_percentiles(scores, frame.query_id.to_numpy(), eligible)[positions]
    result["footprint_area_A2"] = frame.iloc[positions]["surface.footprint_area"].to_numpy()
    result["has_CA"] = frame.iloc[positions].ca_resolved.to_numpy(bool)
    result["flags"] = [";".join(flag for flag, present in (
        ("missing_native_footprint", area <= 0), ("missing_CA", not has_ca), ("tied_score", low != high)) if present)
        for area, has_ca, low, high in zip(result.footprint_area_A2, result.has_CA,
                                           result.rank_min, result.rank_max, strict=True)]
    for name, values in components.items():
        result[name] = np.asarray(values)[positions]
    result = result.sort_values(["rank_min", "residue_key"]).reset_index(drop=True)
    partner_frame = frame.loc[frame.side == side]
    coverage = {"partner_chains": ",".join(selected), "n_structure_residues": int((all_features.side == side).sum()),
                "n_WT_interface_context_nodes": len(partner_frame), "n_ranked_candidates": len(result),
                "n_native_alanine_context_nodes": int((partner_frame.wt == "A").sum()),
                "n_other_excluded_interface_nodes": int(len(partner_frame) - len(result) - (partner_frame.wt == "A").sum()),
                "n_candidates_without_native_footprint": int((result.footprint_area_A2 <= 0).sum()),
                "n_candidates_without_CA": int((~result.has_CA).sum())}
    metadata = {
        "model_dir": str(Path(model_dir).resolve()), "model_kind": bundle.metadata["kind"],
        "representation": bundle.metadata.get("representation", Path(model_dir).name),
        "model_id": bundle.metadata.get("model_id"), "model_version": bundle.metadata.get("version"),
        "model_scope": bundle.metadata.get("scope"),
        "structure": str(structure), "complex_id": complex_id, "selected_partner": list(selected),
        "partner_groups": [list(group) for group in groups], "feature_protocol": protocol,
        "feature_protocol_source": "model.json" if "feature_protocol" in bundle.metadata else "current_defaults",
        "prepared_interface": prepared.metadata(), "feature_columns": requirements["columns"],
        "auxiliary_features": reports, "embeddings": embedding_report,
        "n_full_WT_interface_context_nodes": len(frame), "candidate_coverage": coverage,
        "candidate_population": "all standard non-alanine WT residues of the selected partner within 6 A",
        "feature_context_population": "complete WT interface on both partners including native Ala and unmeasured sites",
        "score_type": "uncalibrated within-partner ranking score; higher ranks first",
        "rank_intervals": "descending minimum/maximum ranks within exact score ties",
        "full_candidate_percentile": "ascending midrank percentile (rank-.5)/candidate_count within partner",
        "component_ranks": "equal-weight fusion components expressed as full-candidate midrank percentiles",
        "labels_used": False, "affinity_plugin_used": False,
    }
    result.to_csv(output / "rankings.csv", index=False)
    pd.DataFrame([coverage]).to_csv(output / "candidate_coverage.csv", index=False)
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return result, metadata
