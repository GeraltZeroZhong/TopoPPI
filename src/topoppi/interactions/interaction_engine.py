"""Generic sidecar generation and compatible interaction-file APIs.

Chemical implementations are discovered through the default Chemistry plugin;
reading a saved interaction document never requires that plugin or its backend.
"""
from __future__ import annotations

import json
import logging
import os
from datetime import datetime, timezone
from pathlib import Path

from topoppi import __version__
from topoppi.file_utils import sha256_file
from topoppi.interactions.records import (
    INTERACTION_SCHEMA_VERSION,
    PROLIF_TO_STANDARD,
    InteractionResult,
    load_interaction_document,
    load_interaction_partner_map,
    normalize_interaction_name,
    residue_sequence_token,
)
from topoppi.json_utils import dump_json_atomic

# Existing public decoding imports remain available at their historical path.
__all__ = ["INTERACTION_SCHEMA_VERSION", "PROLIF_TO_STANDARD", "InteractionResult",
           "generate_interactions", "generate_prolif_interactions", "load_interaction_document",
           "load_interaction_partner_map", "load_prolif_document", "load_prolif_partner_map",
           "normalize_interaction_name", "residue_sequence_token", "compute_interactions"]
logger = logging.getLogger("InteractionEngine")


def compute_interactions(interface, *, backend="prolif", context=None):
    """Compute typed interactions from Core atoms using the selected backend."""
    from topoppi import run

    return run("chemistry.detect", interface, backend=backend, context=context)


def generate_interactions(pdb_path, chain_a, chain_b, log=None, *, backend="prolif",
                          source_sha256=None, output_dir=None, chain_data=None):
    """Generate or reuse a bound sidecar, reusing caller-selected atoms if given.

    The hash is the established schema-4 source-binding/cache contract. Pure
    object callers should use compute_interactions without hashing or file IO.
    """
    log = log or logger
    target_dir = os.path.abspath(output_dir or os.path.dirname(os.path.abspath(pdb_path)))
    pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]
    chain_tag = f"{str(chain_a).strip()}-{str(chain_b).strip()}".replace(os.sep, "_")
    backend_tag = str(backend).replace(os.sep, "_")
    out_path = os.path.join(target_dir, f"{pdb_name}.{chain_tag}.{backend_tag}.json")
    current_sha = source_sha256
    if os.path.exists(out_path):
        try:
            with open(out_path, "r", encoding="utf-8") as handle:
                existing = json.load(handle)
            source_sha = existing.get("source_sha256") if isinstance(existing, dict) else None
            schema_version = existing.get("interaction_schema_version") if isinstance(existing, dict) else None
            engine = existing.get("engine") if isinstance(existing, dict) else None
        except (OSError, json.JSONDecodeError):
            source_sha, schema_version, engine = None, None, None
        current_sha = current_sha or sha256_file(pdb_path)
        if source_sha == current_sha and schema_version == INTERACTION_SCHEMA_VERSION and engine == backend:
            log.info("Found existing %s output: %s", backend, os.path.basename(out_path))
            return out_path
        log.warning("Existing interaction output has a stale source binding or schema; regenerating.")

    from topoppi import ExecutionContext, get_interface, load_structure
    from topoppi.interface import Interface, PartnerAtoms

    if chain_data is None:
        interface = get_interface(load_structure(pdb_path), chain_a, chain_b)
    else:
        coords_a, atoms_a, coords_b, atoms_b = chain_data
        interface = Interface("chemical-interactions", {"path": os.path.abspath(pdb_path)}, (
            PartnerAtoms((str(chain_a).strip(),), coords_a, atoms_a),
            PartnerAtoms((str(chain_b).strip(),), coords_b, atoms_b),
        ))
    result = compute_interactions(interface, backend=backend, context=ExecutionContext(logger=log))
    if not isinstance(result, InteractionResult):
        raise TypeError("A chemical backend must return an InteractionResult.")
    if (result.engine, result.chain_a, result.chain_b) != (backend, str(chain_a).strip(), str(chain_b).strip()):
        raise ValueError("Chemical backend result does not match the selected backend and chains.")
    current_sha = current_sha or sha256_file(pdb_path)
    # Preserve the existing sidecar schema. Rich in-memory provenance is retained
    # by InteractionResult; the disk record binds to its original source here.
    payload = {
        "engine": result.engine,
        "interaction_schema_version": INTERACTION_SCHEMA_VERSION,
        "topoppi_version": __version__,
        "created_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "source_pdb": os.path.abspath(pdb_path), "source_sha256": current_sha,
        "chain_a": result.chain_a, "chain_b": result.chain_b,
        "interactions": result.interactions,
    }
    Path(target_dir).mkdir(parents=True, exist_ok=True)
    dump_json_atomic(payload, out_path)
    log.info("%s JSON generated successfully: %s", backend, os.path.basename(out_path))
    return out_path


def generate_prolif_interactions(pdb_path, chain_a, chain_b, log=None, *,
                                 source_sha256=None, output_dir=None, chain_data=None):
    """Compatibility entry point for the default ProLIF backend."""
    return generate_interactions(pdb_path, chain_a, chain_b, log, backend="prolif",
                                 source_sha256=source_sha256, output_dir=output_dir, chain_data=chain_data)


def load_prolif_document(json_path):
    """Decode historical ProLIF sidecars without importing a chemistry backend."""
    return load_interaction_document(json_path)


def load_prolif_partner_map(json_path, atoms_a, atoms_b, *, expected_chain_a=None,
                            expected_chain_b=None, expected_source_sha256=None, require_bindings=False):
    """Compatibility reader retaining the historical ProLIF engine constraint."""
    return load_interaction_partner_map(json_path, atoms_a, atoms_b, expected_engine="prolif",
        expected_chain_a=expected_chain_a, expected_chain_b=expected_chain_b,
        expected_source_sha256=expected_source_sha256, require_bindings=require_bindings)
