"""Independent fixtures and explicit local plugin registration for source tests."""
import pytest


@pytest.fixture(autouse=True)
def affinity_registration(monkeypatch):
    from topoppi import plugins
    from topoppi_chemistry.plugin import plugin as chemistry_plugin

    from topoppi_affinity.plugin import plugin

    registry = plugins.PluginRegistry(entry_points=())
    registry.register_plugin(chemistry_plugin)
    registry.register_plugin(plugin)
    monkeypatch.setattr(plugins, "_default_registry", registry)
    return registry
