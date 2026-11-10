"""Exact-affinity parsing, mutation notation and structure identities."""
from __future__ import annotations

import re
from functools import lru_cache

import numpy as np
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB import PDBParser

KD_PATTERN = re.compile(r'^Kd=([0-9]+(?:\.[0-9]*)?(?:[eE][+-]?[0-9]+)?)([fpnum]?M)$')
MUTATION_PATTERN = re.compile(r'^([ACDEFGHIKLMNPQRSTVWY])([A-Za-z0-9])(-?\d+)([A-Za-z]?)([ACDEFGHIKLMNPQRSTVWY])$')


def parse_kd(raw: str) -> float:
    """Parse exact PDBbind Kd, rejecting Ki/IC50, inequalities and approximate Kd."""
    match = KD_PATTERN.fullmatch(str(raw).strip())
    if match is None:
        raise ValueError(f'Not an exact Kd with a supported molar unit: {raw}')
    value = float(match[1]) * {'M': 1, 'mM': 1e-3, 'uM': 1e-6, 'nM': 1e-9, 'pM': 1e-12, 'fM': 1e-15}[match[2]]
    if not np.isfinite(value) or value <= 0:
        raise ValueError('Kd must be finite and positive.')
    return value


def exact_positive(raw) -> float:
    text = str(raw).strip()
    if re.fullmatch(r'[+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?', text) is None:
        raise ValueError('Expected an exact positive numerical measurement.')
    value = float(text)
    if not np.isfinite(value) or value <= 0:
        raise ValueError('Expected finite positive measurement.')
    return value


def parse_mutations(raw: str) -> list[tuple[str, str, int, str, str]]:
    result = []
    for token in str(raw).split(','):
        explicit = re.fullmatch(r'([^:\s,]+):([ACDEFGHIKLMNPQRSTVWY])(-?\d+)([A-Za-z]?)([ACDEFGHIKLMNPQRSTVWY])', token.strip())
        match = explicit or MUTATION_PATTERN.fullmatch(token.strip())
        if match is None:
            raise ValueError(f'Unsupported mutation notation: {token}')
        if explicit:
            chain, wt, number, insertion_code, mt = match.groups()
        else:
            wt, chain, number, insertion_code, mt = match.groups()
        result.append((wt, chain, int(number), insertion_code, mt))
    positions = [(chain, number, insertion_code) for _, chain, number, insertion_code, _ in result]
    if len(positions) != len(set(positions)):
        raise ValueError('A mutation set contains the same position more than once.')
    return sorted(result, key=lambda value: (value[1], value[2], value[3]))


def canonical_mutations(raw: str) -> str:
    return ','.join(f'{wt}{chain}{number}{icode}{mt}' if len(chain) == 1 and not chain.isdigit()
                    else f'{chain}:{wt}{number}{icode}{mt}'
                    for wt, chain, number, icode, mt in parse_mutations(raw))


def complex_key(pdb_id: str, chains_a: list[str], chains_b: list[str]) -> str:
    sides = sorted([','.join(sorted(chains_a)), ','.join(sorted(chains_b))])
    return f'{pdb_id.lower()}:{sides[0]}:{sides[1]}'


@lru_cache(maxsize=4096)
def read_structure(path: str) -> tuple[dict[str, str], dict[tuple, str]]:
    model = next(PDBParser(QUIET=True).get_structure('affinity', path).get_models())
    sequences, identities = {}, {}
    for chain in model:
        sequence = []
        for residue in chain:
            letter = protein_letters_3to1_extended.get(residue.resname.upper())
            if letter is None:
                continue
            if 'CA' not in residue:
                continue
            sequence.append(letter)
            identities[(chain.id, residue.id[1], residue.id[2].strip())] = letter
        if sequence:
            sequences[chain.id] = ''.join(sequence)
    return sequences, identities


def attach_structure(row: dict, mutations: str = '') -> dict:
    sequences, identities = read_structure(row['structure_path'])
    chains_a, chains_b = row['chains_a'].split(','), row['chains_b'].split(',')
    if set(chains_a).intersection(chains_b):
        raise ValueError('overlapping_partner_chains')
    if any(chain not in sequences for chain in chains_a + chains_b):
        raise ValueError('missing_partner_chain')
    row['sequence_a'] = '|'.join(sequences[chain] for chain in chains_a)
    row['sequence_b'] = '|'.join(sequences[chain] for chain in chains_b)
    if mutations:
        for wt, chain, number, icode, _mt in parse_mutations(mutations):
            if chain not in chains_a + chains_b:
                raise ValueError('mutation_outside_partner_groups')
            actual = identities.get((chain, number, icode))
            if actual is None:
                raise ValueError('mutation_residue_missing_in_structure')
            if actual != wt:
                raise ValueError(f'mutation_wildtype_mismatch:{chain}{number}{icode}:{wt}!={actual}')
    return row
