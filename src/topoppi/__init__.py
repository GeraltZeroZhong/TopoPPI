"""TopoPPI scientific platform and optional plugin execution API.

Structure loading, contact interfaces, native geometry and plugin discovery are
available independently. Numerical and optional task modules load on first use.
"""

from importlib import import_module

from topoppi._version import __version__
from topoppi.config import (
    BenchmarkConfig,
    OptCutsConfig,
    ParameterizationConfig,
    SurfaceConfig,
    TopologyConfig,
    TopoPPIRunConfig,
    VisualizationConfig,
)

__all__ = [
    "__version__",
    "BenchmarkConfig",
    "OptCutsConfig",
    "ParameterizationConfig",
    "SurfaceConfig",
    "TopologyConfig",
    "TopoPPIRunConfig",
    "VisualizationConfig",
    "Structure",
    "StructureChain",
    "Interface",
    "PartnerAtoms",
    "PartnerInterface",
    "PreparedInterface",
    "GeometryRepresentation",
    "load_structure",
    "get_interface",
    "compute_geometry",
    "compute_interactions",
    "InteractionResult",
    "prepare_interface",
    "save_interface",
    "load_interface",
    "geometric_contact_partner_map",
    "source_atom_residue_labels",
    "vertex_area_masses",
    "residue_dual_boundaries",
    "run",
    "plugins",
    "PluginMetadata",
    "PluginRegistry",
    "ExecutionContext",
]

_LAZY_EXPORTS = {
    "Structure": ("topoppi.structure", "Structure"),
    "StructureChain": ("topoppi.structure", "StructureChain"),
    "load_structure": ("topoppi.structure", "load_structure"),
    "Interface": ("topoppi.interface", "Interface"),
    "PartnerAtoms": ("topoppi.interface", "PartnerAtoms"),
    "PartnerInterface": ("topoppi.interface", "PartnerInterface"),
    "PreparedInterface": ("topoppi.interface", "PreparedInterface"),
    "get_interface": ("topoppi.interface", "get_interface"),
    "prepare_interface": ("topoppi.interface", "prepare_interface"),
    "GeometryRepresentation": ("topoppi.geometry", "GeometryRepresentation"),
    "compute_geometry": ("topoppi.geometry", "compute_geometry"),
    "compute_interactions": ("topoppi.interactions.interaction_engine", "compute_interactions"),
    "InteractionResult": ("topoppi.interactions.records", "InteractionResult"),
    "save_interface": ("topoppi.interface_io", "save_interface"),
    "load_interface": ("topoppi.interface_io", "load_interface"),
    "geometric_contact_partner_map": ("topoppi.contacts", "geometric_contact_partner_map"),
    "source_atom_residue_labels": ("topoppi.residues", "source_atom_residue_labels"),
    "vertex_area_masses": ("topoppi.mesh.metrics", "vertex_area_masses"),
    "residue_dual_boundaries": ("topoppi.mesh.metrics", "residue_dual_boundaries"),
    "run": ("topoppi.plugins", "run"),
    "plugins": ("topoppi.plugins", None),
    "PluginMetadata": ("topoppi.plugins", "PluginMetadata"),
    "PluginRegistry": ("topoppi.plugins", "PluginRegistry"),
    "ExecutionContext": ("topoppi.execution", "ExecutionContext"),
}


def __getattr__(name):
    if name not in _LAZY_EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attribute = _LAZY_EXPORTS[name]
    module = import_module(module_name)
    value = module if attribute is None else getattr(module, attribute)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__))
