"""Lightweight registration; feature and model imports wait for execution."""

from topoppi.plugins import PluginMetadata

from . import __version__


class HotspotPlugin:
    metadata = PluginMetadata(
        name="hotspot", version=__version__,
        description="Frozen-model ranking of residues on a selected interface partner.",
        capabilities=("hotspot.rank",),
    )

    def register(self, registry):
        registry.register_capability("hotspot.rank", "topoppi_hotspot.inference:rank_interface")
        registry.register_command("hotspot", "topoppi_hotspot.cli:main")
        registry.register_gui("hotspot", "topoppi_hotspot.gui:create_controller", label="Hotspot ranking")


plugin = HotspotPlugin()
