"""Native Gaussian surface geometry computed from a selected interface."""

from __future__ import annotations

from topoppi.config import SurfaceConfig, TopologyConfig
from topoppi.errors import PipelineError
from topoppi.interface import Interface, PartnerInterface, PreparedInterface

# Existing imports and native record consumers retain the same concrete class.
GeometryRepresentation = PreparedInterface


def prepare_partner_surface(partner, other_coords, *, surface=None, topology=None):
    """Generate the established Gaussian surface and retained face domains."""
    from topoppi.mesh.surface import SurfaceGenerator
    from topoppi.mesh.topology import TopologyManager

    generator = SurfaceGenerator(partner.coords, config=surface or SurfaceConfig())
    mesh = generator.generate_mesh()
    if mesh is None or not len(mesh.vertices):
        raise PipelineError(f"Failed to generate surface mesh: {generator.last_report}")
    manager = TopologyManager(mesh, other_coords, config=topology or TopologyConfig())
    partner.mesh = mesh
    partner.patches = list(manager.get_interface_patches())
    partner.surface_report = dict(generator.last_report)
    partner.topology_report = dict(manager.last_report)
    return partner


def compute_geometry(interface: Interface, *, surface=None, topology=None, sides=(0, 1),
                     log=None, cancel_event=None) -> GeometryRepresentation:
    """Compute native geometry without rereading the source or invoking OptCuts.

    ``sides=(0,)`` computes only the first partner surface. Empty retained patch
    sets remain coverage information for downstream consumers. The input
    interface's partner orientation and contact rule are preserved.
    """
    if not isinstance(interface, Interface):
        raise TypeError("compute_geometry requires an Interface returned by get_interface.")
    partners = tuple(PartnerInterface(p.chains, p.coords, p.atoms) for p in interface.partners)
    geometry = PreparedInterface(interface.analysis_id, interface.source, partners,
        surface or SurfaceConfig(), topology or TopologyConfig(), interface.numbering, interface)
    for side in sides:
        if cancel_event is not None and cancel_event.is_set():
            raise InterruptedError("Interface preparation cancelled")
        if log:
            log(f"Preparing native interface for {','.join(geometry.partners[side].chains)}")
        prepare_partner_surface(geometry.partners[side], geometry.partners[1 - side].coords,
                                surface=geometry.surface_config, topology=geometry.topology_config)
    return geometry
