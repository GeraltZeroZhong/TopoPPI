# Example inputs

## Antibody–lysozyme example

The installed plugin includes PDB structure 1BVK and a CSV with two mutation sets. Follow the [installation guide](../docs/workflow.md#install-and-prepare), then export the example:

```sh
topoppi-affinity example --output affinity-demo
```

Run the bundled affinity model and open the result:

```sh
topoppi-affinity example --output affinity-demo --run
topoppi-affinity gui affinity-demo/results
```

The example uses antibody chains A,B against lysozyme chain C. Repeated runs save results in a new directory and print its path. Native-surface review is available immediately; OptCuts adds optional atlases. The [model guide](../docs/models.md#applicability) describes model scope, and the [workflow guide](../docs/workflow.md#variant-input) explains mutation input.

For a desktop start, select **Try 1BVK example...**, choose a folder and select **Predict** to use the bundled model. Advanced settings allow a compatible local model or a task model prepared from an explicit download catalog. The panel fills the absolute task and both partner groups. Reuse preserves previous results; edited example files can be analyzed directly or kept while a fresh copy is exported elsewhere.

## Synthetic row-error fixture

`tiny_complex.pdb` is a small synthetic structure shared with the core tests. Its `variants.csv` contains two alternative substitutions at one site, a variant with substitutions in both partners and an absent position. Use it to exercise valid and invalid rows with the bundled mutation model:

```sh
topoppi-affinity analyze --structure tiny_complex.pdb --chains-a A --chains-b B \
  --task mutation --variants variants.csv --output software-fixture
topoppi-affinity gui software-fixture
```

Run these commands from this directory. The absent position produces a row error and exit code 2; valid predictions remain in the result. Add `--model` to use a compatible [custom model](../docs/workflow.md#custom-models).
