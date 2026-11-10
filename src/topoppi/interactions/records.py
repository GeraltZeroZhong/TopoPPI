"""Provider-independent interaction records and legacy sidecar decoding.

Historical ProLIF names are file-format aliases only. Chemical perception and
backend execution belong to the default Chemistry plugin.
"""
from __future__ import annotations

import json
import re
from collections import defaultdict
from dataclasses import dataclass, field

from topoppi.interactions.metadata import INTERACTION_TYPES

INTERACTION_SCHEMA_VERSION = 4


@dataclass
class InteractionResult:
    """Canonical typed residue contacts from one chemical backend.

    The current sidecar contract represents a single author chain per partner.
    Residue sequence tokens include insertion codes. Metadata may record source
    provenance without requiring the original file to remain available.
    """

    engine: str
    chain_a: str
    chain_b: str
    interactions: list[dict]
    metadata: dict = field(default_factory=dict)

    def to_document(self):
        return {**self.metadata, "engine": self.engine, "chain_a": self.chain_a,
                "chain_b": self.chain_b, "interactions": self.interactions}


PROLIF_TO_STANDARD = {
    "HydrogenBond": "HydrogenBond",
    "HBAcceptor": "HydrogenBond",
    "HBDonor": "HydrogenBond",
    "Hydrophobic": "Hydrophobic",
    "PiStacking": "PiStacking",
    "FaceToFace": "PiStacking",
    "EdgeToFace": "PiStacking",
    "PiCation": "PiCation",
    "CationPi": "PiCation",
    "Cationic": "Ionic",
    "Anionic": "Ionic",
    "Ionic": "Ionic",
    "HalogenBond": "HalogenBond",
    "XBAcceptor": "HalogenBond",
    "XBDonor": "HalogenBond",
    "MetalAcceptor": "MetalCoordination",
    "MetalDonor": "MetalCoordination",
    "MetalCoordination": "MetalCoordination",
    "VdWContact": "VdWContact",
}

_RESIDUE_TOKEN_RE = re.compile(r"(?<!\d)(-?\d+)([A-Za-z]?)")


def normalize_interaction_name(name):
    if not name:
        return None
    compact = str(name).replace(" ", "").replace("-", "")
    return compact if compact in INTERACTION_TYPES else PROLIF_TO_STANDARD.get(compact, "Other")


def residue_sequence_token(value):
    """Return a canonical PDB residue sequence token such as ``42`` or ``42A``."""
    match = _RESIDUE_TOKEN_RE.search(str(value).strip())
    if match is None:
        return None
    sequence = int(match.group(1))
    insertion = match.group(2).upper()
    return f"{sequence}{insertion}"


def load_interaction_document(json_path):
    """Read canonical records or a saved document without loading a backend."""

    if isinstance(json_path, InteractionResult):
        payload = json_path.to_document()
    elif isinstance(json_path, dict):
        payload = json_path
    else:
        with open(json_path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError("Interaction JSON root must be an object.")
    interactions = payload.get("interactions", [])
    if not isinstance(interactions, list):
        raise ValueError("Interaction JSON interactions must be a list.")
    residue_types = defaultdict(set)
    partners = defaultdict(dict)
    for item in interactions:
        if not isinstance(item, dict):
            continue
        res_a = residue_sequence_token(item.get("res_a_seq"))
        res_b = residue_sequence_token(item.get("res_b_seq"))
        i_type = normalize_interaction_name(item.get("interaction"))
        if res_a is None or not i_type:
            continue
        residue_types[res_a].add(i_type)
        if res_b is not None:
            partners[res_a][res_b] = partners[res_a].get(res_b, 0) + 1
    return payload, dict(residue_types), dict(partners)


def _residue_label_lookup(atoms):
    """Map canonical PDB residue tokens to TopoPPI residue labels."""

    lookup = {}
    for atom in atoms:
        residue = atom.get_parent()
        chain = residue.get_parent()
        token = residue_sequence_token(f"{int(residue.id[1])}{str(residue.id[2]).strip().upper()}")
        label = f"{chain.id}:{residue.get_resname()}:{token}"
        previous = lookup.setdefault(token, label)
        if previous != label:
            raise ValueError(f"Residue token {token!r} is ambiguous within Chain {chain.id!r}.")
    return lookup


def load_interaction_partner_map(
    json_path,
    atoms_a,
    atoms_b,
    *,
    expected_engine=None,
    expected_chain_a=None,
    expected_chain_b=None,
    expected_source_sha256=None,
    require_bindings=False,
):
    """Return the authoritative Interaction Chain-A-to-B residue partner map.

    Only records whose Chain A and Chain B residue tokens resolve against the
    supplied structure are retained. Values count distinct interaction records
    for a residue pair. Residue weights use the number of partner keys.
    """

    payload, _residue_types, token_partners = load_interaction_document(json_path)
    file_engine = str(payload.get("engine") or "").strip().lower()
    file_chain_a = str(payload.get("chain_a") or "").strip()
    file_chain_b = str(payload.get("chain_b") or "").strip()
    file_source_sha256 = str(payload.get("source_sha256") or payload.get("input_sha256") or "").strip().lower()
    requested_a = str(expected_chain_a or "").strip()
    requested_b = str(expected_chain_b or "").strip()
    requested_source_sha256 = str(expected_source_sha256 or "").strip().lower()
    if require_bindings and (not file_chain_a or not file_chain_b or not file_source_sha256):
        raise ValueError("Formal interaction JSON must declare chain_a, chain_b, and source_sha256/input_sha256.")
    if expected_engine and file_engine and file_engine != expected_engine:
        raise ValueError(f"Interaction JSON engine must be {expected_engine!r}, got {file_engine!r}.")
    if file_chain_a and requested_a and file_chain_a != requested_a:
        raise ValueError(f"Interaction chain_a mismatch: {file_chain_a} != {requested_a}")
    if file_chain_b and requested_b and file_chain_b != requested_b:
        raise ValueError(f"Interaction chain_b mismatch: {file_chain_b} != {requested_b}")
    if file_source_sha256 and requested_source_sha256 and file_source_sha256 != requested_source_sha256:
        raise ValueError("Interaction source checksum does not match the selected structure.")

    labels_a = _residue_label_lookup(atoms_a)
    labels_b = _residue_label_lookup(atoms_b)
    partner_map = {}
    for token_a, partners in token_partners.items():
        label_a = labels_a.get(token_a)
        if label_a is None:
            continue
        resolved = {
            labels_b[token_b]: int(count)
            for token_b, count in partners.items()
            if token_b in labels_b and int(count) > 0
        }
        if resolved:
            partner_map[label_a] = dict(sorted(resolved.items()))
    return dict(sorted(partner_map.items()))
