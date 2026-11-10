"""Interface representations, finite encodings and physical descriptors.

Importing this package leaves codec and machine-learning dependencies unloaded.
"""

from importlib import import_module

from topoppi._version import __version__

_EXPORTS = {
    "Patch": ("schema", "Patch"),
    "ComplexSource": ("schema", "ComplexSource"),
    "Surface": ("schema", "Surface"),
    "Encoded": ("schema", "Encoded"),
    "from_geometry": ("adapters", "from_geometry"),
    "physical_descriptors": ("descriptors.physical", "physical_descriptors"),
    "residue_physics": ("descriptors.physical", "residue_physics"),
    "physical_targets": ("descriptors.surface", "physical_targets"),
    "global_spatial_descriptors": ("descriptors.spatial", "global_spatial_descriptors"),
    "SpatialDescriptorCache": ("descriptors.spatial", "SpatialDescriptorCache"),
    "prepare_spatial_source": ("spatial", "prepare_spatial_source"),
    "encode_spatial_source": ("spatial", "encode_spatial_source"),
    "prepare_spatial_descriptors": ("spatial", "prepare_spatial_descriptors"),
    "encode_uv": ("codecs.uv", "encode_uv"),
    "decode_uv": ("codecs.uv", "decode_uv"),
    "encode_native": ("codecs.native", "encode_native"),
    "decode_native": ("codecs.native", "decode_native"),
    "encode_direct": ("codecs.direct", "encode_direct"),
    "decode_direct": ("codecs.direct", "decode_direct"),
    "save_source": ("source", "save_source"),
    "load_source": ("source", "load_source"),
    "native_surface": ("source", "native_surface"),
    "reference_quadrature": ("source", "reference_quadrature"),
    "canonical_source_bytes": ("source", "canonical_source_bytes"),
}
__all__ = ["__version__", *_EXPORTS]


def __getattr__(name):
    if name not in _EXPORTS:
        raise AttributeError(name)
    module, attribute = _EXPORTS[name]
    value = getattr(import_module(f"{__name__}.{module}"), attribute)
    globals()[name] = value
    return value
