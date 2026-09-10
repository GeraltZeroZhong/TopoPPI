# TopoPPI benchmark evidence schema (v2.0)

_A guide to the benchmark reports and evidence artifacts written by TopoPPI 2.1._

---

## Start here

Each benchmark writes one evidence bundle under `BenchmarkConfig.output_root`.
The bundle contains a human-readable summary, machine-readable details, row-level
audit records, resume state, and file checksums.

This guide covers the TopoPPI `2.1` working-tree version. Within
`benchmark_report.json`, `topoppi_version` identifies the application and
`schema_version` identifies the serialized layout. Their current values are
`2.1` and `2.0`, respectively; each is versioned independently for readers and analysis scripts.

All JSON files use strict JSON syntax. Undefined and non-finite scientific
values are serialized as `null`, which allows standard JSON parsers to read the
files directly.

### Choose the right file

| Question | Open this artifact |
| --- | --- |
| What was run and what completed? | `benchmark_report.json` |
| Which structures succeeded or failed? | `benchmark_summary.csv` |
| Which inputs and chain pairs entered the run? | `benchmark_manifest.csv` |
| Where did a structure or method fail? | `benchmark_failures.csv` |
| How much geometry and biology survived preparation? | `benchmark_per_patch.csv` |
| How did individual residue footprints behave? | `benchmark_per_residue.csv.gz` |
| Which source elements produced the final mesh? | `benchmark_provenance.csv.gz` |
| Which OptCuts command and binary produced a patch? | `benchmark_optcuts_executions.jsonl.gz` |
| Can this bundle be checked for file changes? | `benchmark_artifact_checksums.json` |

### Bundle layout

```text
output_root/
|-- benchmark_report.json
|-- benchmark_summary.csv
|-- benchmark_manifest.csv
|-- benchmark_failures.csv
|-- benchmark_per_patch.csv
|-- benchmark_per_face_sample.csv
|-- benchmark_per_residue.csv.gz
|-- benchmark_provenance.csv.gz
|-- benchmark_optcuts_executions.jsonl.gz
|-- benchmark_checkpoint.json
|-- benchmark_artifact_checksums.json
`-- worker_logs/
```

## Execution profiles and comparison sets

TopoPPI separates comparative evaluation from operational timing. The selected
profile is recorded in `config.execution_profile` and in each structure row.

| Profile | Methods executed | Main use |
| --- | --- | --- |
| `comparative` | Five parameterizations plus selected OptCuts arms | Shared-domain quality and method comparisons |
| `operational_optcuts` | One automatic OptCuts arm | End-to-end runtime and usable-output rate |

The comparative parameterizations are `lscm`, `harmonic`, `slim`, `spherical`,
and `cylindrical`. Available OptCuts arms are `optcuts_automatic`,
`optcuts_lscm_initialized`, and the complete TopoPPI arm identified internally
as `residue_aware_optcuts`.

### Structure status terms

| Term | Meaning |
| --- | --- |
| Attempted | The input received a benchmark job or a recorded preprocessing outcome |
| Valid | At least one common prepared patch reached the comparative domain |
| Complete comparison | Every configured method returned finite metrics on the exact common source faces |
| Complete TopoPPI pair | Geometry-only OptCuts and TopoPPI returned an exact finite pair on the matched domain |
| Failed | A preprocessing, structure, topology, solver, timeout, or resource error was recorded |
| Right censored | Runtime reached a method budget, worker timeout, or memory limit before a usable completion |

Failures remain in attempted-structure counts and failure-rate summaries.
Complete-case quality summaries use the relevant exact comparison domain.
Operational summaries use every observed non-warm-up run.

### Comparison domains

Every prepared patch carries stable source face, vertex, and atom identifiers
plus a source-face hash. Comparative metrics use the same prepared 3D faces for
all methods in a comparison. A complete domain requires:

- every expected patch exactly once;
- matching source-face identifiers and hashes;
- finite primary metric values;
- method-specific feasibility checks recorded in the method status.

Three paired analyses have their own domain records:

- `comparison_domain` for the standard same-domain comparison;
- `initialization_comparison_domain` for automatic and LSCM-initialized OptCuts;
- `residue_aware_comparison_domain` for geometry-only OptCuts and TopoPPI.

The topology-preparation ablation has a separate paired record. Its results stay
outside the standard same-domain method comparison.

## `benchmark_report.json`

### Top-level fields

| Field | Contents |
| --- | --- |
| `schema_version` | Machine-readable report layout, currently `2.0` |
| `topoppi_version` | TopoPPI application version, currently `2.1` |
| `created_at` | UTC creation time |
| `config` | Complete serialized `BenchmarkConfig` |
| `runtime` | Execution model, environment, Git state, resources, hashes, and time budgets |
| `metric_protocol` | Machine-readable metric definitions and evidence filenames |
| `preprocessing` | Accepted and excluded inputs, chain selection, integrity checks, and grid estimates |
| `files` | One detailed record per attempted structure |
| `summary` | Counts, distributions, paired analyses, retention, runtime, and atlas summaries |

### Per-structure records

Each item in `files` brings together the evidence needed to interpret one input:

| Record group | Contents |
| --- | --- |
| Input identity | Filename, SHA-256, chain pair, record ID, split, family, and sequence clusters |
| Source metadata | Dataset source, accession, license, structure type, method, resolution, and confidence fields |
| Preprocessing | Surface generation, interface extraction, topology components, and prepared-domain signatures |
| Method execution | Status, time, hashes, initialization, diagnostics, and quality for each configured method |
| Pair quality | Exact geometry-only/TopoPPI, initialization, and topology-ablation domains |
| Retention | Patch-level geometry, residue, contact, hotspot, interaction, and confidence retention |
| Resources | Stage wall time, CPU time, repeated measurements, and process-tree peak RSS |
| Artifact references | Per-face, per-residue, provenance, and OptCuts-execution evidence files |

Method failures include a stage and reason. Successful methods include their
domain signature, quality metrics, injectivity status, and OptCuts feasibility
certificate when applicable.

### Summary groups

The `summary` object includes:

- attempted, valid, complete, incomplete, and failed structure counts;
- method quality distributions on their declared domains;
- all-attempted method execution and reliability counts;
- standard paired comparisons;
- the LSCM-initialization diagnostic;
- `residue_aware_optcuts_comparisons` for the matched TopoPPI analysis;
- the topology-preparation ablation;
- biological-retention and mesh-cardinality summaries;
- multi-patch atlas summaries;
- isolated runtime, peak RSS, termination, and censoring summaries.

Each paired comparison reports its baseline, treatment, metric path, structure
count, cluster count, cluster source, effect estimates, confidence intervals,
and inferential diagnostics.

## Artifact reference

### `benchmark_summary.csv`

This file provides one row for every attempted structure. It is the quickest
place to join benchmark outcomes with dataset metadata.

Generic method columns include the source domain and signature used for each
value. The exact geometry-only/TopoPPI comparison uses the
`residue_aware_pair_*` columns. These paired columns carry the matched domain,
method statuses, quality values, and fragmentation values required for the
TopoPPI comparison.

### `benchmark_manifest.csv`

This is the realized dataset inventory. It records accepted and excluded files,
input hashes, selected chains, chain-selection mode, record and cluster IDs,
dataset provenance, structure/confidence metadata, `prolif_file`,
`prolif_sha256`, and surface-grid estimates.

The realized manifest makes automatic chain selection and preprocessing
decisions visible even when the input configuration used folder discovery.

### `benchmark_failures.csv`

Each row identifies a failure scope, stage, structure, optional patch or method,
and reason. The file covers preprocessing, structure loading, topology
components, shared domains, parameterizations, OptCuts, timeouts, memory limits,
and output validation.

### `benchmark_per_patch.csv`

Each row follows one raw interface component through topology sanitation and
parameterization. It records:

- raw, topology-stage, and final face/vertex/area counts;
- geometry, source-vertex, source-atom, and residue retention;
- geometric contact-pair retention;
- declared hotspot and declared interaction retention;
- pLDDT retention when confidence evidence is available;
- exact denominator and component provenance.

The raw interface component before sanitation supplies the overall retention
denominator. Stage-specific columns isolate topology and parameterization
effects.

### `benchmark_per_face_sample.csv`

This table contains a deterministic sample from every prepared source-face
domain. Available method values include log-stretch distortion, angle error,
log-area error, local orientation, and face-level identifiers. The sampling
seed and sample size are stored in the configuration.

### `benchmark_per_residue.csv.gz`

Each row describes one interface-residue footprint for one method and evidence
domain. Fields cover:

- domain and source-face signatures;
- original and final footprint-component counts;
- footprint mass and fragmentation;
- seam crossings and nonseparating seam crossings;
- interaction and objective weights;
- exact-pair and independent-arm provenance.

### `benchmark_provenance.csv.gz`

This table maps final face and vertex indices back to source faces, source
vertices, and nearest source atoms. It supports residue tracing, retention
audits, and reconstruction of the geometry lineage.

### `benchmark_optcuts_executions.jsonl.gz`

Each JSON line stores one patch-level OptCuts execution. Records include the
command, resolved executable and SHA-256, input/output hashes, initialization,
effective settings, timeout budget, numeric sidecar hash, topology diagnostics,
and the independently recomputed distortion certificate.

### Checkpoint, checksums, and worker logs

`benchmark_checkpoint.json` stores terminal structure attempts under the full
configuration fingerprint. Resume validates this fingerprint together with
input hashes, chain pairs, interaction hashes, Git revision, and executable.

`benchmark_artifact_checksums.json` records SHA-256 and byte count for the core
evidence files. `worker_logs/` stores compact job/result JSON, stdout, stderr,
and compressed row-level detail spools for every warm-up and measured process.
Worker logs are listed separately because their number depends on the execution
plan.

## Metric definitions

### UV representation and alignment

UV coordinates use shape `(F, 3, 2)`, with one coordinate for each face corner.
This representation preserves distinct coordinates on opposite sides of a seam.

Shape, angle, and log-area metrics use one uniform similarity scale per patch to
align total absolute UV area with total 3D area. Symmetric Dirichlet reporting
adds one analytically optimized global scale. Translation, rotation, one global
reflection, and uniform input scaling leave the reported scores unchanged.

Face means and percentiles use original 3D triangle area as weight. A positive-
area face with an invalid Jacobian marks the corresponding structure comparison
as incomplete.

### Geometric quality metrics

For face `f`, let `J_f` be the aligned 2D Jacobian and let `s1` and `s2` be its
singular values.

| Metric | Definition | Reading direction |
| --- | --- | --- |
| Log stretch | `(\|log(s1)\| + \|log(s2)\|) / 2` | Lower values indicate less stretch |
| Angle distortion | Mean absolute error across the three corner angles, in radians | Lower values indicate better angle preservation |
| Log-area distortion | `\|log(A_uv / A_3d)\|` after similarity alignment | Lower values indicate better area preservation |
| Symmetric Dirichlet | `(s1^2 + s2^2 + s1^-2 + s2^-2) / 2` | Identity value `2`; lower values indicate less distortion |
| Flip rate | Area-weighted fraction of faces with non-positive signed UV area after one patch-level orientation choice | Local orientation diagnostic |
| Global injectivity | Polygonal overlap and self-contact checks over the full patch | Whole-map validity diagnostic |

Flip rate describes local triangle orientation. Global injectivity, overlap,
extreme distortion, and residue fragmentation have dedicated fields.

### OptCuts feasibility certificate

For every returned OptCuts patch, TopoPPI recomputes the raw-scale native
symmetric-Dirichlet constraint:

```text
area-weighted mean of ||J||_F^2 + ||J^-1||_F^2
```

Its identity value is `4`. A returned patch passes when the value is finite and
no larger than the configured bound plus the recorded numerical tolerance.
The reporting score above has identity value `2` after global scale
optimization. Both values and their conventions are stored in the execution
evidence.

### Seams

An internal seam is an original internal 3D edge whose incident face-corner UV
assignments are discontinuous. Boundary edges are counted separately.

| Level | Normalized seam length |
| --- | --- |
| Patch | 3D internal seam length divided by `sqrt(patch area)` |
| Structure | Total 3D internal seam length divided by `sqrt(total retained patch area)` |

The structure-level definition remains stable when identical geometry is
partitioned into a different number of charts.

### Residue-footprint fragmentation

Residue footprints are defined on the original face-dual graph. A residue
receives mass `face_area * labelled_corner_count / 3` on each face. Two incident
faces are adjacent for that residue when their shared edge has a labelled
endpoint.

Original connected components provide the zero-fragmentation baseline. For one
component with total mass `M` and seam-induced pieces of mass `m_k`, the score
is:

```text
1 - sum((m_k / M)^2)
```

Reports include mean, footprint-area-weighted, interaction-weighted, and
objective-weighted fragmentation. They also record cycle rank, seam crossings,
and nonseparating seam crossings.

The interaction weight equals the number of distinct partner residues paired
with each mapped residue in the bound ProLIF record. Formal runs require
`prolif_file`, `prolif_sha256`, chain bindings, and a source-structure digest for
every included structure. A non-formal run without declared interaction
evidence records `geometric_fallback` and derives partner degrees from
heavy-atom pairs inside `contact_distance_angstrom`. TopoPPI assigns objective
weight `1 + interaction weight`. The matched geometry-only arm is
`optcuts_automatic`; the complete TopoPPI arm is serialized as
`residue_aware_optcuts`. Both arms use the same automatic initialization and
candidate-search policy.

`tools/publication/prepare_manifest_prolif.py` generates these records in bulk and writes the two manifest fields. Relative `prolif_file` paths are resolved from the manifest directory.

### Atlas packing

Each chart receives an area-matching scale followed by an optional 90-degree
rotation and a translation. The configured gap supplies chart padding.
Triangle and polygon unions determine covered area, within-chart overlap,
between-chart overlap, utilization, waste, minimum gap, and padding violations.

Multi-patch atlas summaries use structures with more than one common retained
patch. The report records the reference method and every applied chart
transform.

## Biological retention and structure metadata

### Retention stages

Retention is reported at three scopes:

| Scope | Numerator | Denominator |
| --- | --- | --- |
| Overall | Final prepared patch | Raw interface component |
| Topology | Patch after topology sanitation | Raw interface component |
| Parameterization | Final prepared patch | Patch after topology sanitation |

Geometry, atoms, residues, contacts, declared hotspots, declared interactions,
and optional confidence values share these stage definitions. Summary blocks
provide both per-component ratio distributions and pooled component-incidence
totals. Rejected components remain represented in all-attempted retention.

### Experimental and predicted structure fields

| Field family | Interpretation |
| --- | --- |
| `resolution_angstrom` | One official scalar `resolution_combined` value from the experimental record |
| `confidence_metric` | Residue-level confidence metric for predicted coordinates |
| `confidence_source` | Coordinate field or source that supplied confidence |
| `confidence_threshold` | Frozen inclusion threshold, when configured |
| `confidence_stratum` | Mean-pLDDT group: `[0,70)`, `[70,90)`, or `[90,100]` |
| `afdb_ipsae` | AFDB complex-level ipSAE value |
| `afdb_ipsae_stratum` | ipSAE group: `[0,0.50)`, `[0.50,0.70)`, or `[0.70,1]` |
| `paired_geometry_*` | Contact, ligand RMSD, clash, mapping, and paired-geometry evidence |

Formal predicted-structure inputs use `plddt_bfactor` and require one finite
value in `[0,100]` for every retained heavy atom in both selected chains.
Residue-level pLDDT and complex-level ipSAE retain separate fields and strata.

`resolution_angstrom` is populated when the source provides one official scalar
resolution. NMR and records with another resolution representation retain the
original source metadata and leave this scalar field empty.

Complex-level ipTM, ipSAE, pDockQ, and LIS fields are populated for actual AFDB
complex coordinates. Independently aligned AFDB monomer replacements and their
experimental references leave these complex-level fields empty.

### Sequence mapping and fixed-pose replacements

When a chain-specific SIFTS interval is unavailable, sequence-based mapping uses
at least `0.90` selected-pair consensus across every optimal correspondence.
Ambiguous repeat positions are recorded as mapping attrition.

Fixed-pose monomer replacement requires:

- at least ten sequence-matched C-alpha atoms;
- sequence identity of at least `0.70`;
- non-collinear matched coordinates for a unique rigid fit.

Whole-chain alignment coverage is recorded as a descriptor. Mapping-completeness
and paired-geometry strata capture partial AFDB domain availability. High and
moderate paired-geometry strata also require selected-pair consensus of at
least `0.90` across every optimal sequence correspondence.

## Statistical summaries

### Analysis units

Face observations contribute to area-weighted structure metrics. Structures
supply the descriptive unit. Paired inference first computes one difference per
structure and then averages within the resolved dependency cluster.

| Cohort | Primary clustering key |
| --- | --- |
| Experimental structures | `family_id` |
| Predicted structures | `inference_family_id` |
| Mixed experimental and predicted structures | `analysis_split_component_id` |

Outside formal mode, missing family IDs fall back to `cluster_id` and then to a
single-structure cluster. Reports expose `cluster_source_counts` so readers can
see which key supplied each analysis unit.

### Intervals and tests

The primary interval bootstraps resolved cluster means with the configured seed
and iteration count. Each paired block also reports:

- structure-level and cluster-level means and differences;
- relative improvement and its declared reference value;
- a cluster-level standardized paired effect;
- a two-sided Wilcoxon signed-rank test on resolved cluster means;
- a shared-protein dyadic-robust sensitivity when partner-cluster metadata is
  complete;
- an exact discordant-pair diagnostic for binary reliability endpoints.

The dyadic sensitivity covers heterotypic interaction families with two distinct
partner nodes. Homotypic families remain in the primary family-cluster analysis.
The sensitivity reports asymptotic and finite-degree intervals, node degrees,
and the largest leave-one-protein-cluster-out mean shift.

When a numerical p-value underflows to zero, the stored value becomes the
smallest positive `float64` and `wilcoxon_p_value_censored_from_zero` becomes
`true`.

### TopoPPI endpoint family

Within `summary.residue_aware_optcuts_comparisons`:

- `automatic_objective_weighted_fragmentation` is the prespecified primary
  efficacy comparison;
- `automatic_distortion_mean`, `automatic_symmetric_dirichlet_mean`, and
  `automatic_normalized_seam_length` form the supporting trade-off family;
- `automatic_unusable_output` supplies the paired reliability endpoint.

The primary endpoint has one prespecified comparison and receives no
multiplicity adjustment. Supporting comparisons receive Benjamini-Hochberg
adjustment as one family. Standard-method quality, symmetric-Dirichlet,
unusable-output, jointly-injective, initialization, and topology-ablation
comparisons retain their own declared families.

## Runtime and resource evidence

Each structure and repetition runs in a fresh subprocess. The parent process
sets thread variables, assigns a disjoint CPU-affinity block, samples the worker
and its descendants for peak RSS, and records wall and CPU time.

Formal quality runs use one measured repetition and zero warm-ups. Formal
performance runs use one worker, at least one warm-up, and at least three
measured repetitions. Warm-up observations stay outside quality and timing
estimates.

For right-censored runs:

- `runtime_observation_sec` is the elapsed lower bound at the censoring event;
- supervisor wall time includes later process shutdown and evidence writing;
- `termination_reason` records method budget, worker timeout, memory limit, or
  another termination state.

The `operational_optcuts` profile emits an `execution_certificate`. A usable
structure returns every prepared patch on the exact source-face domain, passes
global injectivity, and satisfies the independently recomputed raw-scale
OptCuts bound.

## Failure, resume, and integrity behavior

Preprocessing, structure, topology, method, timeout, memory, and output-validation
failures remain visible in `benchmark_failures.csv` and the all-attempted report
counts.

Formal manifest-integrity errors stop execution before workers start. When
preprocessing accepts zero structures, TopoPPI writes the report, realized
manifest, failure log, empty tabular artifacts, checkpoint, provenance header,
and checksums before returning an error.

Resume preserves terminal successes, failures, and incomplete comparisons. The
checkpoint fingerprint binds the complete configuration, manifest, Git commit,
resolved OptCuts executable, and relevant input evidence. Per-structure hashes,
chain pairs, and interaction hashes are revalidated during resume.

For archival use, retain:

- the benchmark configuration and original manifest;
- the coordinate-audit file and preflight report used for a formal run;
- the code commit and environment export;
- the OptCuts executable hash and build provenance;
- the complete evidence bundle and artifact checksums.

Tables and figures should use one configuration fingerprint and one report
schema. Report attempted, complete, and failed sample counts together with the
analysis unit.

## Sensitivity plans

Each sensitivity scenario writes a standard benchmark evidence bundle with the
same schema described above. A plan also writes `sensitivity_results.json` and
`sensitivity_summary.csv` across scenarios.

Supported axes are interface cutoff, surface grid spacing, Gaussian sigma,
isosurface level, OptCuts initial lambda, and OptCuts distortion bound.
The baseline configuration includes `optcuts_automatic` in `optcuts_variants`.
Add `residue_aware_optcuts` when the study will measure stability of the
TopoPPI treatment effect.
One-factor plans contain one baseline plus each non-baseline perturbation.
Factorial plans use the Cartesian product of configured values. Dataset, chain
selection, seed, metrics, and unvaried settings remain fixed across the plan.
