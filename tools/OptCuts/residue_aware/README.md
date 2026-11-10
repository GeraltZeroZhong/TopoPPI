# Residue-aware OptCuts integration

This directory contains the dependency-free C++ state engine and patch set that add residue-footprint fragmentation energy to the pinned TopoPPI OptCuts build. OptCuts supplies symmetric Dirichlet energy optimization and topology feasibility checks.

## How the objective works

The engine fixes residue footprints on the original face-dual graph. Cutting a seam removes an original dual adjacency. A candidate split or merge is evaluated as one transaction, and the engine recomputes only the residues touched by its edges. Transactional evaluation captures the combined effect of multiple cuts through a cyclic footprint.

Footprint mass is the integral of the residue indicator over the source surface, in Å². For an original footprint component with total area `M` and post-cut piece areas `m_k`, fragmentation is:

```text
1 - sum_k (m_k / M)^2
```

Component scores are weighted by their share of the footprint area. The sidecar then assigns a non-negative objective weight to each residue. TopoPPI uses `1 + contact degree`, where contact degree counts distinct partner residues. ProLIF supplies these partners by default; single-run `--interaction-source geometric` uses the configured heavy-atom contact distance.

The pinned patch maps every OptCuts split or merge path to original edge IDs, adds `alpha * candidateDelta(...)` to the seam-energy change, and commits the selected transaction after the topology operation succeeds. The residue-aware objective and the weight-zero geometry ablation use the same candidate search.

## Sidecar format, version 2

TopoPPI passes the objective state to OptCuts in a deterministic whitespace-delimited file:

```text
TOPOPPI_FOOTPRINT_V2
COUNTS <face_count> <residue_count> <internal_edge_count> <input_vertex_count>
SOURCES <source_vertex_0> ... <source_vertex_V-1>
WEIGHTS <w_0> ... <w_R-1>
FACE <face_id> <entry_count> <residue_id> <mass> ...
EDGE <edge_id> <source_vertex_0> <source_vertex_1> <face_0> <face_1> <initial_cut:0|1> <entry_count> <residue_id> ...
```

The records have these roles:

| Record | Meaning |
| --- | --- |
| `COUNTS` | Expected faces, residues, internal source edges, and OBJ input vertices |
| `SOURCES` | Root surface vertex for each OBJ input vertex, including repaired vertex-fan copies |
| `WEIGHTS` | Dense zero-based residue objective weights |
| `FACE` | Integrated piecewise-linear residue mass on one face |
| `EDGE` | Original source edge, adjacent faces, initial seam state, and continuous residue support |

Every face and internal source edge has one record. Diskification copies share their geometry vertex, while OBJ texture indices carry initial UV seams. The Python measurements use the same integrated corner-indicator mass stored in each `FACE` record.

## Build the executable

From the repository root, run:

```bash
bash tools/OptCuts/build_residue_aware_optcuts.sh
```

Pass an output path as the first argument when a separate build artifact is useful:

```bash
bash tools/OptCuts/build_residue_aware_optcuts.sh /tmp/OptCuts_bin
```

The script checks out the pinned upstream commit, applies the residue and reproducibility patches, builds in a temporary directory, and installs the stripped executable at the requested path. Set `OPTCUTS_BUILD_JOBS` to control build parallelism.

## Mesh and cycle handling

TopoPPI's OBJ writer collapses diskification copies, retains repaired vertex-fan copies, and records the root source ID for every geometry vertex. Supplied UV seams use texture indices. Sparse local solves keep candidate memory bounded.

An A-B-A-B topology cycle is identified from the complete face-corner connectivity partition and cohesive-edge state. OptCuts uses its critical-lambda query after confirming two distinct states with matching seam energy and distortion on each recurrence.

See the [benchmark evidence schema](../../../docs/benchmark_schema.md#residue-footprint-fragmentation) for the formal objective definition and exported measurements. Build provenance and platform distribution details are in the [OptCuts notice](../NOTICE.md).
