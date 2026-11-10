# Migration to TopoPPI 3.0

Affinity is an optional distribution named `topoppi-affinity`. Install it alongside TopoPPI 3.0 to use prediction and linked result review. Scientific tasks are available through the shared **Binding analysis** navigation.

## Shared geometry and saved results

`topoppi.run("affinity.predict", geometry, ...)` accepts an existing Core `GeometryRepresentation`. Path-based callers can continue to use `topoppi_affinity.workflow.analyze`. Core continues to read saved atlas geometry and its annotations independently of the prediction plugin.

Fitted estimator classes, feature names, units and serialized fields retain their numerical contracts. Existing model keys such as `legacy_feature_names` continue to identify their stored feature block.

## Installation and launch commands

`topoppi_affinity.installer` manages Affinity resources. The `topoppi-install-affinity` command forwards to the installed plugin; explicit imports of `topoppi.install_affinity` emit a deprecation warning. The old `topoppi-gui --affinity` and `--affinity-setup` options remain deprecated aliases of plugin launch options.

## Explicit fitting inputs

The feature and model tools accept the files required by each operation:

- `fit_absolute_model`: `--training-table`, `--sequence-features`, `--descriptor-features` and `--output`.
- `refit_frozen_model`: `--training-table`, `--features`, `--output` and the model configuration. `--antisymmetric` makes predictions change sign when the signs of the mutation features are reversed.
- `refit_residual_model`: `--training-table`, `--features`, `--prior`, `--output` and the grouping task.
- `pool_manifest_embeddings`: `--graphs`, `--sequences`, `--manifest`, `--output` and `--target`.
- `extract_analytic_features`: `--manifest`, `--graphs`, `--features`, `--output` and `--target`.

The former `--root`, `--development-only` and `--auxiliary-features` options for these tools are replaced by explicit input files. Concatenate intended training pools before fitting, retain unique sample IDs, and construct evaluation splits independently of the refit.

## Preparation metadata

`cath_protocol.prepare_cath_protocol` reports target differences as `target_mean_absolute_difference_from_reported_temperature`. The temperature-derived values remain in `ddg_at_reported_temperature`.
