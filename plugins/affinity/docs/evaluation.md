# Evaluation primitives

The plugin provides grouped regression metrics, deterministic partitioning, feature ablations and reusable fitting APIs. Supply observations, groups and model parameters explicitly. Fit preprocessing within each training partition. Compare predictors on matching successful rows and retain prediction coverage alongside the scores.

Absolute affinity uses pKd, the negative base-10 logarithm of molar Kd. Mutation effects use the change in binding free energy, ΔΔG, in kcal/mol, with positive values denoting weaker binding. Complex-level summaries give each complex equal weight; row-level summaries describe all observations.

## Grouped fitting and prediction

`topoppi_affinity.grouped_fitting.fit_nested` selects parameters using grouped inner folds. `learning_subset` constructs nested subsets of whole families. `evaluation.regression_metrics` computes regression summaries; `splits` and `benchmark_data` provide grouping and partitioning utilities for supplied tables.

`artifacts.fit_frozen` fits an explicit configuration to a training table and a feature NPZ. The table contains `sample_id`, `complex_id`, `family_group` and the selected target. Feature files contain `sample_ids`, `X` and `feature_names`. Rows are matched by sample ID, and the model stores its fitted preprocessing and ordered feature names.

`residual_artifacts.fit_residual_model` additionally accepts a prior prediction table containing `sample_id,predicted,target,method,status,error`. It selects radial basis function (RBF) kernel parameters within two training folds grouped by mutation or family. The saved model applies the same feature definition and prior method during inference.

Use held-out predictions to assess predictive performance. A model refitted for application carries its training metadata separately from the evaluation results. The synthetic test suite checks numerical invariants, training-only preprocessing, row alignment and serialization.
