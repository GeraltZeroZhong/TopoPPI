"""Lightweight plugin declaration; importing it loads no model dependencies."""
from topoppi.plugins import PluginMetadata

from topoppi_affinity import __version__


class AffinityPlugin:
    metadata = PluginMetadata(
        name="affinity", version=__version__,
        description="Affinity prediction and linked native-interface review",
        capabilities=("affinity.predict",),
    )

    def register(self, registry):
        registry.register_capability("affinity.predict", "topoppi_affinity.capabilities:predict")
        registry.register_command("affinity", "topoppi_affinity.cli:main")
        registry.register_gui("affinity", "topoppi_affinity.integration:create_controller", label="Affinity prediction")


plugin = AffinityPlugin()
