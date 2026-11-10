# Feature extraction and model tools

Install the core and plugin in Python 3.10. Run a tool as `python -m topoppi_affinity_tools.<group>.<name> --help`; paths, model parameters and output destinations are explicit inputs.

| Package | Function | Commands |
| --- | --- | --- |
| `features` | Native descriptors, residue graphs, analytic curvature and sequence representations | `extract_features`, `extract_graph_features`, `extract_sequence_embeddings`, `extract_analytic_features`, `prepare_analytic_mutation_inputs`, `refine_graph_curvature`, `pool_manifest_embeddings`, `extract_prospective_sasa` |
| `models` | User-selected fitting and named-feature prediction | `fit_absolute_model`, `refit_frozen_model`, `refit_residual_model`, `bundle_mutation_ensemble`, `predict_absolute_model`, `predict_frozen_affinity`, `predict_mutation_ensemble`, `predict_residual_affinity`, `predict_attribution` |

Inspect required inputs with `--help` before executing a command. Neural feature extraction requires its optional dependencies and local checkpoints or explicit resource preparation. Everyday structure analysis and review use `topoppi-affinity analyze`, `inspect` and `review`; see the [workflow](../docs/workflow.md).
