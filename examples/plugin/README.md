# Minimal TopoPPI plugin

This example plugin reports counts from a shared geometry object.

The example depends only on `topoppi`, whose default distribution includes Core
and Chemistry. From a checkout with TopoPPI installed:

```bash
python -m pip install ./examples/plugin
topoppi plugins list
topoppi plugins info stats
python -m topoppi.examples demo --synthetic
topoppi stats demo/synthetic_complex.pdb --chains-a A --chains-b B
```

```python
import topoppi

geometry = topoppi.compute_geometry(
    topoppi.get_interface(topoppi.load_structure("demo/synthetic_complex.pdb"), "A", "B")
)
print(topoppi.run("stats.describe", geometry))
```

The package metadata declares the `topoppi.plugins` entry point. Its provider
registers lazy capability and command references; the registry supplies an
execution context when running the capability. Local development may use
`PluginRegistry.register_plugin(plugin)` with the same provider object. Start a new
Python process after installing another plugin, or refresh registry discovery.

See [the plugin contract](../../docs/architecture.md#plugin-contract) for
configuration, errors and compatibility boundaries.
