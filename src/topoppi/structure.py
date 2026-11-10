"""Loaded protein atoms and author identities, independent of surface geometry."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np


@dataclass
class StructureChain:
    """Selected protein heavy atoms for one author chain, in structure order."""

    chain_id: str
    coords: np.ndarray
    atoms: list
    residue_count: int


@dataclass
class Structure:
    """An in-memory first-model structure with coordinates in angstroms.

    Derived interfaces reuse selected atom objects and assemble their coordinate
    arrays from the cached chains. Consumers editing atoms should copy them.
    ``source`` retains file provenance after the input file is removed.
    """

    source: dict
    chains: tuple[StructureChain, ...]
    numbering: str = "supplied_structure_author_chain_residue_insertion"

    @property
    def chain_ids(self) -> tuple[str, ...]:
        return tuple(chain.chain_id for chain in self.chains)

    def get_protein_chain_ids(self) -> list[str]:
        return list(self.chain_ids)

    def get_chain_data(self, chain_id):
        for chain in self.chains:
            if chain.chain_id == chain_id:
                return chain.coords, chain.atoms
        raise ValueError(f"Chain {chain_id} not found in structure.")

    def get_chain_residue_count(self, chain_id) -> int:
        for chain in self.chains:
            if chain.chain_id == chain_id:
                return chain.residue_count
        raise ValueError(f"Chain {chain_id} not found in structure.")


def _source_metadata(path: Path, *, source="user", complex_id=None) -> dict:
    stat = path.stat()
    return {"path": str(path), "mtime_ns": stat.st_mtime_ns, "size": stat.st_size,
            "source": source, "complex_id": complex_id or path.stem}


def load_structure(path, *, source="user", complex_id=None) -> Structure:
    """Load selected first-model PDB/mmCIF heavy atoms once.

    The established loader retains recognized amino-acid residues and selected
    alternate conformers, excluding hydrogens, water, ions and cofactors.
    Loading performs no surface generation, plotting or model initialization.
    """
    from topoppi.io.io_loader import PDBLoader

    path = Path(path).expanduser().resolve()
    loader = PDBLoader(str(path))
    chains = []
    for chain_id in loader.get_protein_chain_ids():
        coords, atoms = loader.get_chain_data(chain_id)
        chains.append(StructureChain(chain_id, coords, atoms, loader.get_chain_residue_count(chain_id)))
    return Structure(_source_metadata(path, source=source, complex_id=complex_id), tuple(chains))
