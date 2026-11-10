# Antibody–lysozyme example

This example contains [PDB structure 1BVK](https://www.rcsb.org/structure/1BVK), with antibody chains A,B as one partner and lysozyme chain C as the other. The mutation CSV uses residues present in the structure. A row containing multiple substitutions specifies one variant.

Export the example from an installed plugin:

```sh
topoppi-affinity example --output affinity-example
```

Run the bundled CPU affinity model:

```sh
topoppi-affinity analyze --structure affinity-example/1bvk.pdb \
  --chains-a A,B --chains-b C --task absolute --output affinity-result
```

For mutation effects, use `--task mutation` and add
`--variants affinity-example/variants.csv`. Each default model includes its
fitted preprocessing. Use `--model` for a compatible custom estimator. The
[model guide](https://github.com/GeraltZeroZhong/TopoPPI/blob/HEAD/plugins/affinity/docs/models.md#applicability)
describes applicability and advanced neural-model resources.
