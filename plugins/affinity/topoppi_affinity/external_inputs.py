"""Explicit structure-number mapping for published predictors requiring 1-based chains."""

from __future__ import annotations

import json
import string
from pathlib import Path

from Bio.Data.PDBData import protein_letters_3to1, protein_letters_3to1_extended
from Bio.PDB import PDBIO, MMCIFParser, PDBParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from topoppi_affinity.data import parse_mutations


def chain_aliases(chains):
    """Preserve PDB chain IDs; map longer mmCIF author IDs into free PDB IDs."""
    chains = sorted(set(chains))
    alphabet = string.ascii_uppercase + string.ascii_lowercase + string.digits
    aliases = {chain: chain for chain in chains if len(chain) == 1 and chain in alphabet}
    available = iter(char for char in alphabet if char not in aliases.values())
    for chain in chains:
        if chain not in aliases:
            try:
                aliases[chain] = next(available)
            except StopIteration as exc:
                raise ValueError("This external prior accepts at most 62 partner chains in its PDB input") from exc
    return aliases


def external_partner_groups(chains_a, chains_b):
    groups = [[chain.strip() for chain in group.split(",")] for group in (chains_a, chains_b)]
    aliases = chain_aliases(sum(groups, []))
    return tuple("".join(aliases[chain] for chain in group) for group in groups)


def prepare_numbered_complex(row, output_path):
    """Keep complete partner chains and coordinates; number observed residues from 1.

    Missing positions and insertion codes are explicitly mapped, not interpreted
    as sequential offsets. Only canonical protein residues with CA are accepted.
    """
    output_path = Path(output_path)
    source_path = Path(row["structure_path"])
    parser = MMCIFParser(QUIET=True) if source_path.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    source = parser.get_structure("source", str(source_path))[0]
    structure, model = Structure("external"), Model(0)
    structure.add(model)
    mapping, entries = {}, []
    groups = [[chain.strip() for chain in row[field].split(",")] for field in ("chains_a", "chains_b")]
    aliases = chain_aliases(sum(groups, []))
    sequences = {}
    for chain_id in sorted(set(sum(groups, []))):
        if chain_id not in source:
            raise ValueError(f"Partner chain is absent from the input structure: {chain_id}")
        chain = Chain(aliases[chain_id])
        model.add(chain)
        sequence = ""
        for residue in source[chain_id]:
            if "CA" not in residue:
                continue
            if residue.resname not in protein_letters_3to1_extended:
                continue  # Calcium ions and ligands can also contain an atom named CA.
            aa = protein_letters_3to1.get(residue.resname)
            if aa is None:
                raise ValueError(f"Noncanonical residue unsupported by external input: {chain_id}:{residue.id}:{residue.resname}")
            number = len(sequence) + 1
            sequence += aa
            copied = Residue((" ", number, " "), residue.resname, residue.segid)
            for atom in residue:
                if atom.element not in {"H", "D"}:
                    copied.add(atom.copy())
            chain.add(copied)
            key = (chain_id, residue.id[1], residue.id[2].strip())
            mapping[key] = (number, aa)
            entries.append({"chain": chain_id, "external_chain": aliases[chain_id],
                            "source_number": residue.id[1], "source_insertion": residue.id[2].strip(),
                            "number": number, "aa": aa})
        if not sequence:
            raise ValueError(f"Partner chain contains no canonical protein residues with CA: {chain_id}")
        sequences[chain_id] = sequence
    for label, group in zip(("sequence_a", "sequence_b"), groups, strict=True):
        if label in row and row[label] != "|".join(sequences[c] for c in group):
            raise ValueError("External input sequences do not match the declared partner sequences")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(output_path))
    output_path.with_suffix(".numbering.json").write_text(json.dumps(entries, indent=2) + "\n")
    return mapping


def remap_mutations(mutations, mapping):
    result = []
    aliases = chain_aliases(key[0] for key in mapping)
    for wt, chain, position, insertion, mt in parse_mutations(mutations):
        number, actual = mapping[(chain, position, insertion)]
        if actual != wt:
            raise ValueError(f"WT identity mismatch: {chain}:{position}{insertion}")
        result.append(f"{wt}{aliases[chain]}{number}{mt}")
    return ",".join(result)
