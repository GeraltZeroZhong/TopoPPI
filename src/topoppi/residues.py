"""Author residue identities shared by structure and geometry consumers."""

from __future__ import annotations

from typing import Sequence

import numpy as np


def atom_residue_label(atom) -> str:
    """Return the stable chain:name:sequence label, including insertion codes."""
    residue = atom.get_parent()
    chain = residue.get_parent()
    insertion = str(residue.id[2]).strip()
    return f"{chain.id}:{residue.get_resname()}:{residue.id[1]}{insertion}"


def source_atom_residue_labels(atoms: Sequence[object]) -> np.ndarray:
    """Return residue labels aligned with the supplied selected atom array."""
    return np.asarray([atom_residue_label(atom) for atom in atoms], dtype=object)
