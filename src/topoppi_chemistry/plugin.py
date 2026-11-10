"""Default lightweight chemistry provider declaration."""

from topoppi.plugins import PluginMetadata
from topoppi_chemistry import __version__


class ChemistryPlugin:
    metadata = PluginMetadata(
        name="chemistry", version=__version__,
        description="Chemical interaction perception with discoverable backends",
        capabilities=("chemistry.detect", "chemistry.backends", "chemistry.classify_contact"),
    )

    def register(self, registry):
        registry.register_capability("chemistry.detect", "topoppi_chemistry.api:detect")
        registry.register_capability("chemistry.backends", "topoppi_chemistry.api:backends")
        registry.register_capability("chemistry.classify_contact", "topoppi_chemistry.annotations:classify_contact")
        registry.register_command("chemistry", "topoppi_chemistry.cli:main")
        registry.register_gui("chemistry", "topoppi_chemistry.gui:create_controller", label="Chemical interactions")


plugin = ChemistryPlugin()
