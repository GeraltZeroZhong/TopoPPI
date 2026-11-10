"""Heavy-atom residue contact geometry, independent of surfaces and layouts."""

from __future__ import annotations

from collections import defaultdict
from typing import Sequence

import numpy as np
from scipy.spatial import cKDTree

from topoppi.residues import source_atom_residue_labels


def geometric_contact_partner_map(
    coords_a: np.ndarray,
    atoms_a: Sequence[object],
    coords_b: np.ndarray,
    atoms_b: Sequence[object],
    *,
    distance_cutoff: float,
) -> dict[str, dict[str, int]]:
    """Count selected heavy-atom pairs within the cutoff for each residue pair."""
    points_a = np.asarray(coords_a, dtype=np.float64)
    points_b = np.asarray(coords_b, dtype=np.float64)
    if len(points_a) == 0 or len(points_b) == 0:
        return {}
    labels_a = source_atom_residue_labels(atoms_a)
    labels_b = source_atom_residue_labels(atoms_b)
    neighborhoods = cKDTree(points_b).query_ball_point(points_a, r=float(distance_cutoff))
    partners: dict[str, dict[str, int]] = defaultdict(dict)
    for index_a, indices_b in enumerate(neighborhoods):
        label_a = str(labels_a[index_a])
        counts = partners[label_a]
        for index_b in indices_b:
            label_b = str(labels_b[index_b])
            counts[label_b] = counts.get(label_b, 0) + 1
    return {label_a: dict(sorted(counts.items())) for label_a, counts in sorted(partners.items()) if counts}
