"""Scientific metadata exposed by model reports and saved analyses.

Fitted objects retain their serialization identifiers. User reports describe the
input representation, training membership and numerical protocol independently
of local experiment directories and experiment bookkeeping.
"""
from __future__ import annotations

REPRESENTATIONS = {
    "legacy_v1": "conventional_descriptors",
    "plm_analytic": "sequence_and_analytic_surface",
    "plm_surface": "sequence_and_surface",
    "plm_contact": "sequence_and_contact",
    "plm_only": "sequence",
}

MODEL_SCOPES = {
    "absolute": "Absolute-affinity regression with sequence and interface descriptor kernels",
    "known": "Mutation-effect ensemble for systems represented in training",
    "unseen": "Mutation-effect ensemble for transfer to unrepresented systems",
    "descriptor": "Conventional interface-descriptor regression",
    "fitted": "Custom fitted affinity regressor",
}

# These fields affect applicability, model interpretation or reproducibility.
_FIELDS = {
    "target", "task", "units", "affinity_units", "family", "parameters", "seed",
    "n_train", "n_train_complexes", "n_complexes", "n_families", "n_features",
    "n_failed_prior", "train_sample_ids", "train_complex_ids", "train_family_groups",
    "training_sample_ids", "validation_sample_ids", "test_sample_ids",
    "use_surface", "origin_slope", "kernel_normalization", "weighting",
    "surface_scope", "symmetry", "missing_prior", "input_contract", "ddg_sign",
    "geometry", "curvature", "prior_method", "source_url", "license", "training_data",
    "training_overlap", "overlap_definition", "input_scope", "mc_samples",
    "backbone_noise_A", "antithetic_variates", "batch_tokens", "sign", "n_attempted", "n_success",
    "field", "derivatives", "replaced_graph_columns", "pooled_readout", "parameter_selection",
    "checkpoint", "checkpoint_scope",
}


def training_metadata(report):
    """Keep factual training/protocol fields while omitting research bookkeeping."""
    result = {key: value for key, value in report.items() if key in _FIELDS}
    if result.get("checkpoint"):
        result["checkpoint"] = str(result["checkpoint"]).replace("\\", "/").rsplit("/", 1)[-1]
    if "parameters" in result:
        result["parameters"] = {
            ("descriptor_kernel_weight" if key == "legacy_weight" else key): value
            for key, value in result["parameters"].items()
        }
    if "representation" in report:
        result["representation"] = REPRESENTATIONS.get(report["representation"], report["representation"])
    if "method" in report:
        result["method"] = (
            "ESMC and conventional interface descriptors with normalized two-kernel regression"
            if report["method"] == "ESM/Topo plus legacy normalized two-kernel KRR" else report["method"]
        )
    if "combination" in report:
        result["combination"] = report["combination"].replace("frozen RT-scaled", "RT-scaled")
    if "selection" in report:
        result["parameter_selection"] = (
            "Two inner folds grouped by variant or family according to the task; minimum complex-macro MAE"
            if "two inner folds" in report["selection"].lower()
            else "Fixed ensemble weights" if "equal weights" in report["selection"].lower()
            else "Fixed parameters stored with the fitted estimator"
        )
    if "preprocessing" in report and isinstance(report["preprocessing"], dict):
        result["preprocessing"] = training_metadata(report["preprocessing"])
    for name in ("base", "residual", "prior_provenance"):
        if isinstance(report.get(name), dict):
            result[name] = training_metadata(report[name])
    if "overlap_definition" in result:
        result["overlap_definition"] = result["overlap_definition"].replace("final model flags", "binding-supervised model flags")
    if "prior_provenance" in result:
        prior = result["prior_provenance"]
        checkpoint = report["prior_provenance"].get("checkpoint")
        if checkpoint:
            prior["checkpoint"] = str(checkpoint).replace("\\", "/").rsplit("/", 1)[-1]
        if prior.get("method") == "StaB-ddG-stability-only":
            prior["checkpoint_scope"] = "Folding-stability-supervised author checkpoint; structural pretraining exposure applies"
    return result


def public_model_info(info):
    """Present saved model records through the same scientific vocabulary."""
    result = {key: info[key] for key in ("kind", "target", "scope", "curvature", "training_scope_note", "task", "model_version", "feature_scheme") if key in info}
    result["kind"] = {"pilot": "descriptor", "frozen": "fitted"}.get(result.get("kind"), result.get("kind"))
    if result["kind"] in MODEL_SCOPES:
        result["scope"] = MODEL_SCOPES[result["kind"]]
    result["training"] = training_metadata(info.get("training", {}))
    return result
