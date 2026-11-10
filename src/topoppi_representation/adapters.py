"""Adapt existing Core geometry and face-corner layouts to representation inputs."""

from __future__ import annotations

import numpy as np

from .chemistry import aa_from_residue_key
from .schema import ComplexSource, Patch


def from_geometry(geometry, *, layouts=None, residue_fields=None, sample_id=None):
    """Reuse native geometry and explicitly supplied UV arrays.

    ``layouts[(side, patch_index)]`` maps layout names to arrays of shape
    ``(n_faces, 3, 2)`` in the corresponding Core patch's face-corner order.
    Names are ``ordinary``, ``aware`` (contact-weighted), or ``aware_uniform``.
    This operation preserves caller partner order and rounds native XYZ to the
    common float32 coordinate contract. It performs no UV optimization.

    ``residue_fields`` maps a field name to side-qualified residue keys, such as
    ``side0:A:ALA:12``. Entries without a surface footprint are kept in the residue
    table. Every supplied field must cover that complete table.
    """
    from topoppi.atlas.footprints import mesh_vertex_residue_labels
    from topoppi.mesh.provenance import OPTCUTS_GEOMETRY_VERTEX_IDS
    from topoppi.residues import source_atom_residue_labels

    layouts = layouts or {}
    keys, amino_acids, lookup, patches = [], [], {}, []

    def register(key):
        if key not in lookup:
            lookup[key] = len(keys)
            keys.append(key)
            amino_acids.append(aa_from_residue_key(key))
        return lookup[key]

    max_rounding = 0.0
    for side, partner in enumerate(geometry.partners):
        atom_labels = source_atom_residue_labels(partner.atoms)
        for index, mesh in enumerate(partner.patches):
            labels = np.asarray(mesh_vertex_residue_labels(mesh, atom_labels), dtype=str)
            identity = np.asarray(mesh.metadata.get(OPTCUTS_GEOMETRY_VERTEX_IDS, np.arange(len(mesh.vertices))))
            _, first, inverse = np.unique(identity, return_index=True, return_inverse=True)
            if not np.array_equal(labels, labels[first][inverse]):
                raise ValueError("Geometry copies disagree in residue identity")
            xyz = np.asarray(mesh.vertices)[first].astype(np.float32).astype(float)
            if len(xyz):
                max_rounding = max(max_rounding, float(np.abs(np.asarray(mesh.vertices)[first] - xyz).max()))
            uv = {}
            for name, array in layouts.get((side, index), {}).items():
                array = np.asarray(array, dtype=float)
                if array.shape != (len(mesh.faces), 3, 2):
                    raise ValueError("UV arrays must follow the native patch face-corner order")
                uv[name] = array.copy()
            patches.append(
                Patch(
                    side,
                    index,
                    xyz,
                    inverse[np.asarray(mesh.faces, dtype=int)],
                    np.asarray([register(f"side{side}:{key}") for key in labels[first]], dtype=int),
                    uv,
                )
            )
    for field in (residue_fields or {}).values():
        for key in field:
            register(key)
    fields = {
        name: np.asarray([values[key] for key in keys], dtype=float) for name, values in (residue_fields or {}).items()
    }
    return ComplexSource(
        sample_id=sample_id or geometry.analysis_id,
        residue_keys=keys,
        amino_acids=amino_acids,
        patches=patches,
        residue_fields=fields,
        metadata={
            "coordinate_units": "angstrom",
            "native_mass_units": "angstrom_squared",
            "common_coordinate_precision": "float32",
            "max_coordinate_rounding_A": max_rounding,
            "chain_groups": [list(partner.chains) for partner in geometry.partners],
        },
    )
