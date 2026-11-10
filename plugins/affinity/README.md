# TopoPPI Affinity

TopoPPI Affinity adds affinity modelling and linked interface review to [TopoPPI](../../README.md). It extracts native interface descriptors and includes CPU models for binding affinity (pKd) and mutation effects (ΔΔG). Users can fit regression models to supplied data or load custom models, then review saved predictions alongside the corresponding residues and contacts.

## Installation

Use Python 3.10 and install the core and plugin from the repository root:

```sh
python -m pip install -e . -e './plugins/affinity[dev]'
```

The optional `models` extra supplies neural representation dependencies. Install it when using a model that requires them. Tk is required for the optional GUI. OptCuts is required for atlas generation. Native geometry and prediction run independently of atlas generation.

## Platform capability

Installation registers `affinity.predict`, `topoppi affinity ...`, and the
optional GUI integration through `topoppi.plugins`. Model dependencies load when
needed for prediction. New analyses can share a Core geometry object:

```python
import topoppi

geometry = topoppi.compute_geometry(
    topoppi.get_interface(topoppi.load_structure("complex.pdb"), "A", "B")
)
result = topoppi.run("affinity.predict", geometry, task="absolute", output="results")
```

The capability reuses selected atoms and native surfaces without reopening the
input. `workflow.analyze` and the standalone CLI also accept structure paths.
Namespaced defaults use `ExecutionContext.plugin_config["affinity"]`; explicit
arguments take precedence.

The Core GUI lists **Binding affinity** and **Mutation effects** in the shared
**Binding analysis** navigation, alongside Hotspot ranking when that plugin is
installed. The prediction panel has two scientific tasks: whole-complex pKd and
mutation ΔΔG. Keep **Prepared task model** to use the bundled physical model.
Advanced known/unseen ensembles retain their separately documented training
scope; **Local model file** selects a compatible custom estimator.

Selecting a local model reads its target and declared scope. A later task or
scope mismatch is reported before computation. Variant drafts remain available
when switching tasks; absolute-affinity prediction uses the unmodified complex.
Selected interface residues can be sent to the mutation editor as separate
single-substitution variants or one variant containing multiple substitutions.
Adding variants appends rows, retains existing edits and resolves duplicate row
identifiers. An existing
draft for another structure or partner grouping must be switched explicitly.

Task UI and resource setup live in `topoppi_affinity.integration` and
`topoppi_affinity.installer`. The old `topoppi-install-affinity` command is a
deprecated forwarding adapter and requires the plugin to be installed first.

## Synthetic quick start

This example converts a dissociation constant to pKd and back:

```python
from topoppi_affinity.thermodynamics import kd_to_pkd, pkd_to_kd

assert float(kd_to_pkd(1e-8)) == 8.0
assert float(pkd_to_kd(8.0)) == 1e-8
```

The `examples/tiny_complex.pdb` fixture contains independently constructed toy coordinates. Its variant CSV exercises alternative substitutions, a variant containing multiple substitutions and an absent residue. Use it to check input handling and row-error reporting.

## Structure example and prediction

The installed example contains the documented [PDB 1BVK structure](https://www.rcsb.org/structure/1BVK) and example variants. It illustrates the antibody–lysozyme interface. Export it with:

```sh
topoppi-affinity example --output affinity-demo
```

Analyze it with the bundled affinity model:

```sh
topoppi-affinity analyze --structure affinity-demo/1bvk.pdb \
  --chains-a A,B --chains-b C --task absolute \
  --output affinity-demo/results
topoppi-affinity inspect affinity-demo/results
topoppi-affinity gui affinity-demo/results
```

The plugin contains one default model per task, including fitted preprocessing
and feature definitions. For mutation effects, use `--task mutation` and provide
`--mutations` or `--variants`. Default prediction runs on CPU without model
downloads. See [models](docs/models.md) and the [workflow](docs/workflow.md).

## Inputs and outputs

Supply PDB or mmCIF structures and complete partner chain groups. Mutation input uses `variant_id,mutations` CSV records or structured residue fields. Affinity targets are pKd or mutation ΔΔG in kcal/mol; positive ΔΔG means weaker binding. Temperature is in kelvin and the dissociation constant Kd is in molar units.

Saved analyses contain a JSON report, prediction tables, the native interface and any requested layout results. Prediction rows retain their identifiers, status and error messages. The report records model scope and separate prediction/layout status. Models determine their feature requirements and applicability; review their training metadata before interpreting predictions.

Run `topoppi-affinity --help` for CLI commands. See the [workflow guide](docs/workflow.md), [model interfaces](docs/models.md), [evaluation primitives](docs/evaluation.md), [feature tools](tools/README.md) and [test instructions](tests/README.md).

## Contributing

Use independently constructed fixtures for unit tests, preserve target units and row identities, and add coverage for changed public behavior. Run the affinity test suite together with affected core tests. Public interfaces are listed in the [API guide](docs/api.md).
