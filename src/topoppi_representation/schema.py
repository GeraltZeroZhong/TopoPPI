"""Scientific carrier shared by all complete codecs (length Å, mass Å²)."""

from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import csr_matrix


@dataclass
class Patch:
    side: int
    patch_id: int
    xyz: np.ndarray  # native, identity-welded vertices [N,3]
    faces: np.ndarray  # [F,3] indexing xyz
    residue_index: np.ndarray  # [N] into ComplexSource residue table
    uv: dict = field(default_factory=dict)  # ordinary/aware -> [F,3,2], raw corners


@dataclass
class ComplexSource:
    sample_id: str = ""
    family_id: str = ""
    residue_keys: list[str] = field(default_factory=list)  # side-prefixed unique keys
    amino_acids: list[str] = field(default_factory=list)  # one letter, X unknown
    patches: list[Patch] = field(default_factory=list)
    residue_fields: dict = field(default_factory=dict)  # per-residue arrays, e.g. SASA Å²
    metadata: dict = field(default_factory=dict)


@dataclass
class Surface:
    xyz: np.ndarray  # decoded quadrature/cell positions, common 3D frame
    edges: np.ndarray  # undirected [E,2], within physical patch
    mass: csr_matrix  # [N,R] source mass per cell and residue (Å²)
    side: np.ndarray  # [N], 0 or 1
    patch: np.ndarray  # [N], physical patch ids; distinct per side
    residue_keys: list[str]
    amino_acids: list[str]
    residue_fields: dict = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)
    faces: np.ndarray | None = None  # available for decoded native meshes, not UV cells

    @property
    def weights(self):
        return np.asarray(self.mass.sum(axis=1)).ravel()


@dataclass
class Encoded:
    blob: bytes
    metadata: dict = field(default_factory=dict)  # diagnostic only; decoder cannot read it

    @property
    def nbytes(self):
        return len(self.blob)
