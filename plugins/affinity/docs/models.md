# Model interfaces

TopoPPI Affinity 3.0 includes two official CPU inference models: whole-complex
binding affinity (pKd) and mutation effects (ΔΔG in kcal/mol). The models,
fitted preprocessing and ordered feature definitions are installed with the
plugin. Each task uses one final fitted estimator. Model details are in the
[default model card](../topoppi_affinity/data/production/MODEL_CARD.md).

## Default models

The affinity model uses linear ridge regression on 477 physicochemical descriptors.
The mutation model uses kernel ridge regression with a radial basis function (RBF)
kernel. Its physicochemical inputs comprise 263 descriptors of the native
interface and substitutions and 167 atomic
solvent-accessible surface area (SASA) descriptors, plus temperature divided by
298.15 K. Both use median imputation, missingness indicators and standardization.
Model parameters were selected by group-macro mean absolute error across
inner-validation splits.
Final fits use all eligible training samples with inverse similarity-group
weights. The affinity model stores regression coefficients and an intercept;
the mutation model stores standardized training vectors and dual coefficients.
Both include fitted preprocessing parameters. Cross-validation performance is
estimated from separate held-out predictions.

```sh
topoppi-affinity analyze --structure complex.pdb --chains-a A --chains-b B \
  --task absolute --output affinity-result
topoppi-affinity analyze --structure complex.pdb --chains-a A --chains-b B \
  --task mutation --mutations A:G42V --output mutation-result
```

`models status` uses local files to report model availability. The default
models require CPU computation and the installed scientific dependencies. Optional
2D atlas generation additionally requires OptCuts.

## Applicability

Provide a prepared protein complex and the complete interacting chain groups.
The affinity output is pKd = −log10(Kd/M). Positive mutation ΔΔG indicates weaker
mutant binding. Mutation predictions describe single or multiple substitutions on
the fixed wild-type structure; the atomic SASA context is calculated
from that structure. Review sequence identity, missing residues and conformation
when interpreting predictions. Per-prediction confidence intervals are unavailable.

## Custom models and upstream resources

Use `--model` with a compatible `model.json` or `.joblib` artifact. Load joblib
models from trusted sources. `AffinityRegressor` stores an estimator, named feature
schema, target and training metadata. Dictionary bundles support named-feature
regressors, sequence/descriptor kernels and mutation ensembles with declared priors.
`inference.model_info` describes the required representation and target.

Custom neural models use the plugin's `models` extra. ESMC and upstream priors
retain their source locations, revisions and license information in the resource
catalog. For advanced known/unseen ensembles, provide an explicit catalog using
`models download --catalog CATALOG.json`, or select a compatible local model.
These identifiers describe training-system applicability. Select them according
to the estimator's documented training scope. `--cache-dir` chooses the resource
cache. The default physical models use descriptors calculated from the supplied
structure.
