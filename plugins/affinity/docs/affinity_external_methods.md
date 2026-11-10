# External prediction adapters

The optional adapters run upstream methods using explicit input manifests, source repositories and model resources. Each upstream package retains its installation requirements, model license and input conventions. Run adapters in a compatible environment and keep the source revision and checkpoint identity with the output.

## Absolute affinity

| Module and interface | Method | Required inputs |
| --- | --- | --- |
| `external_absolute.run_prodigy` | PRODIGY | Structure manifest, output directory and upstream resource directory |
| `external_absolute.run_mirage` | MIRAGE | Structure manifest, upstream source/model directory and sequence-embedding resources |
| `author_references.run_bindpred_features` | BindPred | Input manifest, named embedding features and fitted checkpoint |
| `author_references.run_proaffinity_graphs` | ProAffinity-GNN | Structure manifest, prepared graphs and their partner index, upstream repository and checkpoint |

The structure manifest identifies `sample_id`, `pdb_id`, `structure_path`, `chains_a` and `chains_b`. Supply complete partner chain groups. Energy-based outputs are converted to pKd using the method's stated sign convention and gas constant. `external_absolute.native_dg_to_pkd` exposes this conversion directly.

PRODIGY and MIRAGE process all supplied rows by default (`subset="all"`).
To select a held-out partition, set `subset="heldout"` and provide
`split_complex` and `split_sequence` columns. This selects rows marked `test`
in either column.

For ProAffinity-GNN, pass `graph_root` and `partner_index` for your prepared graphs.
The index records the two partner chain groups for each graph; they must match
the manifest. Each tab-separated line contains the graph key and two
semicolon-separated partner groups, with commas separating chains within a group.

## Mutation effects

| Module and interface | Method | Required inputs |
| --- | --- | --- |
| `external_mutation.run_rde` | RDE-Net | Mutation input table and upstream repository/checkpoints |
| `external_mutation.run_stabddg` | StaB-ddG | Mutation input table, upstream repository and chosen checkpoint |
| `expanded_references.run_mutation_reference` | PPIformer, Pythia-PPI, SAAMBE-3D | Mutation manifest, repository and method identifier |
| `cached_author_references.run_cached_reference` | DDAffinity, GearBind | Mutation manifest, compatible author cache, repository and checkpoint |
| `external_energy.run_evoef` | EvoEF1 | Mutation manifest and the upstream executable |

Mutation adapters preserve complete partner groups, wild-type residue identities and chain-renumbering maps. Outputs use ΔΔG in kcal/mol with positive values indicating weaker mutant binding. Supported substitution scope and missing structural inputs are recorded per row. A checkpoint's training provenance remains relevant when interpreting its predictions.

## Input and output records

Prediction tables expose sample identity, prediction, target, method and execution status. Mutation tables use the common fields `sample_id,predicted,target,method,status,error`; additional columns retain method-specific information. The output metadata records the selected software and model resources.

Pass explicit input tables to run the intended samples. The RDE and StaB-ddG adapters accept `records_frame` for in-memory records. Select checkpoints according to their documented applicability, and use matched successful sample IDs when comparing predictions from different methods.

For normal TopoPPI prediction, see the [workflow guide](workflow.md). For metric definitions and grouped fitting, see [evaluation](evaluation.md).
