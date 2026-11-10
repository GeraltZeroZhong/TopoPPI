"""Selected partner interfaces and compatible native geometry records."""
from __future__ import annotations

import copy
import math
import uuid
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

import numpy as np

from topoppi.config import SurfaceConfig, TopologyConfig, VisualizationConfig
from topoppi.contacts import geometric_contact_partner_map
from topoppi.errors import InputDataError
from topoppi.residues import source_atom_residue_labels
from topoppi.structure import Structure, _source_metadata, load_structure


def chain_group(value):
    chains = tuple(sorted(set(str(v).strip() for v in (value.split(",") if isinstance(value, str) else value))))
    if not chains or any(not chain for chain in chains):
        raise ValueError("Partner chain groups must be nonempty (comma-separated chain IDs).")
    return chains


def load_partner_data(loader, chains):
    available = loader.get_protein_chain_ids()
    missing = [chain for chain in chains if chain not in available]
    if missing:
        raise InputDataError(f"Selected chain(s) {', '.join(missing)} were not found. "
                             f"Available protein chains: {', '.join(available) or 'none'}.")
    arrays, atoms = [], []
    for chain in chains:
        coords, current = loader.get_chain_data(chain)
        if not len(coords):
            raise InputDataError(f"Chain {chain} has no standard protein atoms.")
        arrays.append(coords)
        atoms.extend(current)
    return np.concatenate(arrays), atoms


@dataclass
class PartnerAtoms:
    """A complete selected partner group and aligned protein heavy atoms."""

    chains: tuple
    coords: np.ndarray
    atoms: list

    @property
    def residue_keys(self):
        return source_atom_residue_labels(self.atoms)


@dataclass
class Interface:
    """Two oriented partner groups and a declared heavy-atom contact rule.

    This representation contains no meshes. It remains usable after the source
    file is removed because selected atom data and provenance are held in memory.
    """

    analysis_id: str
    source: dict
    partners: tuple[PartnerAtoms, PartnerAtoms]
    contact_distance_angstrom: float = 6.0
    numbering: str = "supplied_structure_author_chain_residue_insertion"
    structure: Structure | None = field(default=None, repr=False)
    _contact_map: dict | None = field(default=None, repr=False)

    def partner(self, chains):
        group = chain_group(chains)
        return next(partner for partner in self.partners if partner.chains == group)

    def chain_data(self, side=0):
        left, right = self.partners[side], self.partners[1 - side]
        return left.coords, left.atoms, right.coords, right.atoms

    def contacts(self, side=0, cutoff=None):
        cutoff = self.contact_distance_angstrom if cutoff is None else float(cutoff)
        if not math.isfinite(cutoff) or cutoff <= 0:
            raise ValueError("Contact distance must be finite and positive.")
        if cutoff != self.contact_distance_angstrom:
            return geometric_contact_partner_map(*self.chain_data(side), distance_cutoff=cutoff)
        if self._contact_map is None:
            self._contact_map = geometric_contact_partner_map(*self.chain_data(), distance_cutoff=cutoff)
        if side == 0:
            return {key: dict(values) for key, values in self._contact_map.items()}
        if side != 1:
            raise IndexError("Interface side must be 0 or 1.")
        reverse = {}
        for left, partners in self._contact_map.items():
            for right, count in partners.items():
                reverse.setdefault(right, {})[left] = count
        return {key: dict(sorted(values.items())) for key, values in sorted(reverse.items())}

    def metadata(self):
        return {"analysis_id": self.analysis_id, "source": self.source,
                "groups": [list(partner.chains) for partner in self.partners],
                "numbering": self.numbering, "coordinates": "angstrom_original_structure",
                "contact_distance_angstrom": self.contact_distance_angstrom,
                "contact_definition": "selected_protein_heavy_atom_pair_distance"}


def get_interface(structure: Structure, chains_a, chains_b, *, contact_distance=6.0, analysis_id=None) -> Interface:
    """Select complete partner groups and contacts from an already loaded structure.

    Caller A/B orientation is preserved. This cutoff governs heavy-atom contacts;
    surface-face selection is configured separately during geometry computation.
    """
    if not isinstance(structure, Structure):
        raise TypeError("get_interface requires a Structure returned by load_structure.")
    group_a, group_b = chain_group(chains_a), chain_group(chains_b)
    if set(group_a) & set(group_b):
        raise ValueError("Partner chain groups overlap.")
    cutoff = float(contact_distance)
    if not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("Contact distance must be finite and positive.")
    partners = tuple(PartnerAtoms(group, *load_partner_data(structure, group)) for group in (group_a, group_b))
    interface = Interface(analysis_id or uuid.uuid4().hex, dict(structure.source), partners,
                          cutoff, structure.numbering, structure)
    interface.contacts()
    return interface


@dataclass
class PartnerInterface(PartnerAtoms):
    """Selected partner atoms with their computed native surface and patches."""

    mesh: object = None
    patches: list = field(default_factory=list)
    surface_report: dict = field(default_factory=dict)
    topology_report: dict = field(default_factory=dict)

    def layout_patches(self):
        """Layout may edit topology and metadata, so give it independent meshes."""
        return copy.deepcopy(self.patches)


@dataclass
class PreparedInterface:
    """Native geometry representation, retained under its established type name."""

    analysis_id: str
    source: dict
    partners: tuple[PartnerInterface, PartnerInterface]
    surface_config: SurfaceConfig
    topology_config: TopologyConfig
    numbering: str = "supplied_structure_author_chain_residue_insertion"
    _interface: Interface | None = field(default=None, repr=False, compare=False)

    @property
    def interface(self) -> Interface:
        if self._interface is None:
            partners = tuple(PartnerAtoms(p.chains, p.coords, p.atoms) for p in self.partners)
            self._interface = Interface(self.analysis_id, self.source, partners, numbering=self.numbering)
        else:
            # Retried analyses assign a new analysis ID to reused native geometry.
            self._interface.analysis_id = self.analysis_id
            self._interface.source = self.source
        return self._interface

    def partner(self, chains):
        group = chain_group(chains)
        return next(p for p in self.partners if p.chains == group)

    def chain_data(self, side=0):
        left, right = self.partners[side], self.partners[1 - side]
        return left.coords, left.atoms, right.coords, right.atoms

    def contacts(self, side=0, cutoff=None):
        return self.interface.contacts(side, cutoff)

    def metadata(self):
        metadata = {"analysis_id": self.analysis_id, "source": self.source,
                "groups": [list(p.chains) for p in self.partners], "numbering": self.numbering,
                "surface": asdict(self.surface_config), "topology": asdict(self.topology_config),
                "coordinates": "angstrom_original_structure", "uv_required": False}
        if self.interface.contact_distance_angstrom != 6.0:
            metadata["contact_distance_angstrom"] = self.interface.contact_distance_angstrom
        return metadata

    def visualizer(self, side=0, **kwargs):
        from topoppi.visualization.visualizer import InterfaceVisualizer

        coords_a, atoms_a, coords_b, atoms_b = self.chain_data(side)
        config = replace(kwargs.pop("config", None) or VisualizationConfig(), use_geometric_interaction_fallback=True)
        visualizer = InterfaceVisualizer(atoms_a, coords_a, coords_b, atoms_b,
            chain_a_id=",".join(self.partners[side].chains), chain_b_id=",".join(self.partners[1-side].chains),
            structure_label=self.source.get("complex_id", Path(self.source["path"]).stem),
            interaction_partner_map=self.contacts(side),
            contact_distance_angstrom=self.interface.contact_distance_angstrom, config=config, **kwargs)
        visualizer.interaction_residue_source = "geometric"
        return visualizer


def prepare_partner_surface(partner, other_coords, *, surface=None, topology=None):
    """Compatibility entry point for native surface generation."""
    from topoppi.geometry import prepare_partner_surface as compute_partner_surface

    return compute_partner_surface(partner, other_coords, surface=surface, topology=topology)


def prepare_interface(structure_path, chains_a, chains_b, *, surface=None, topology=None,
                      analysis_id=None, source="user", complex_id=None, sides=(0, 1), chain_data=None,
                      log=None, cancel_event=None):
    """Load, select and compute native geometry using the public Core API.

    Empty interface domains are retained as coverage information. A consumer
    requiring a visible patch reports that restriction in its own branch.
    ``sides=(0,)`` preserves the single-surface cost of the mapping workflow.
    """
    from topoppi.geometry import compute_geometry

    group_a, group_b = chain_group(chains_a), chain_group(chains_b)
    if set(group_a) & set(group_b):
        raise ValueError("Partner chain groups overlap.")
    if chain_data is None:
        structure = load_structure(structure_path, source=source, complex_id=complex_id)
        interface = get_interface(structure, group_a, group_b, analysis_id=analysis_id)
    else:
        path = Path(structure_path).expanduser().resolve()
        coords_a, atoms_a, coords_b, atoms_b = chain_data
        interface = Interface(analysis_id or uuid.uuid4().hex,
            _source_metadata(path, source=source, complex_id=complex_id),
            (PartnerAtoms(group_a, coords_a, atoms_a), PartnerAtoms(group_b, coords_b, atoms_b)))
    return compute_geometry(interface, surface=surface, topology=topology, sides=sides,
                            log=log, cancel_event=cancel_event)
