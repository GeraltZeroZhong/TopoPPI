# Affinity Python API

TopoPPI Affinity provides structure analysis, prediction from named features, model fitting and review of predictions alongside interface residues and contacts. Structural identifiers retain author chain IDs, residue numbers and insertion codes.

## Predict from a structure

The bundled CPU models predict binding affinity and mutation effects. Choose the task and provide the structure and partner chain groups:

```python
from topoppi_affinity.workflow import analyze, load_analysis

report = analyze(
    "complex.pdb", "A", "B",
    task="absolute", output="affinity_result", device="cpu",
)
affinity = load_analysis("affinity_result").predictions
print(affinity)
```

For a substitution, supply the wild-type complex. Prediction uses its fixed coordinates:

```python
report = analyze(
    "complex.pdb", "A", "B",
    task="mutation", mutations=["A:G42V"],
    output="mutation_result", device="cpu",
)
mutation_effects = load_analysis("mutation_result").predictions
print(mutation_effects)
```

Replace `A:G42V` with a substitution matching your structure's chain ID, author residue number, insertion code if present and wild-type amino acid. Each list entry describes one variant; comma-separated substitutions in one entry form a variant with multiple substitutions. Results include prediction units and row-specific status. The returned report records the saved files and overall prediction status. Supply `model=` to use a custom model.

## Reuse Core geometry

Pass a prepared Core `GeometryRepresentation` to the registered capability:

```python
import topoppi

structure = topoppi.load_structure("complex.pdb")
interface = topoppi.get_interface(structure, "A,B", "C")
geometry = topoppi.compute_geometry(interface)
result = topoppi.run(
    "affinity.predict", geometry,
    task="absolute", output="shared_geometry_result",
)
```

The capability reuses the selected atoms and compatible geometry. Use `task="mutation"` with `mutations=["A:G42V"]` for mutation effects, or supply `model="fitted_model.joblib"` for a custom model. The model determines its required features. See the [workflow guide](workflow.md) for saved results and linked review.

## Scientific interfaces

| Module | Interface | Purpose |
| --- | --- | --- |
| `features` | `extract_complex_features`, `mutation_features` | Conventional interface descriptors and mutation features |
| `data` | `parse_mutations`, `canonical_mutations`, `read_structure` | Mutation notation and structural identities |
| `thermodynamics` | `kd_to_pkd`, `pkd_to_kd`, `ddg_from_kd` | Affinity units and temperature-dependent conversion |
| `inference` | `model_info`, `predict_features` | Inspect a loaded model bundle or predict from named feature files |
| `artifacts` | `fit_frozen`, `predict_frozen` | Fit and use a fixed-configuration named-feature regressor |
| `legacy_kernel` | `fit_legacy_bundle`, `predict_legacy_kernel` | Two-block kernel fitting and prediction |
| `residual_artifacts` | `fit_residual_model`, `predict_residual` | Fit or apply a scalar-prior residual model |
| `grouped_fitting` | `fit_nested` | Select model parameters within supplied training groups |
| `evaluation` | `regression_metrics` | Row-weighted and equal-complex regression summaries |

Absolute-affinity predictions use pKd = −log10(Kd/M). Mutation effects use ΔΔG in kcal/mol, with positive values indicating weaker mutant binding. Conversion from pKd to binding free energy uses the supplied temperature, which defaults to 298.15 K.

## Feature and model contracts

Named feature NPZ files contain `sample_ids`, `X` and `feature_names`. Fitting aligns rows by sample ID and stores the ordered feature schema with the fitted model. The two-block kernel additionally requires descriptor `feature_groups`. Prediction checks the feature order against that stored schema.

Fitting functions take explicit observations and feature inputs. Models retain their fitted preprocessing, units and training metadata. Use model files from a trusted source because joblib deserialization reconstructs Python objects. See [model interfaces](models.md), [evaluation](evaluation.md) and [migration notes](migration.md).
