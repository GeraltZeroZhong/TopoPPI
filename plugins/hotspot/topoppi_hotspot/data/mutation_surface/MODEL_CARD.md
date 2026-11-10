# Mutation-derived hotspot ranking model

Model ID: `topoppi-hotspot-mutation-surface`  
Model version: `3.0.0`  
Feature scheme: `mutation_physical_430_temperature_surface_121`

## Intended use

The model predicts the binding free-energy change for alanine substitution at
each standard non-alanine interface residue and ranks these substitutions within
one selected binding partner. Inputs are a protein complex structure and the
complete chain groups forming the two partners. Candidates have at least one
heavy atom within 6 Å of the other partner. Each prediction uses the supplied
wild-type coordinates and a reference temperature of 298.15 K.

Scores are predicted ΔΔG in kcal/mol. Positive values indicate weaker binding
after substitution; larger values rank first. The saved table retains residue
identities, tied rank intervals and surface coverage for review and export.

## Inputs and representation

The model takes 552 inputs: 430 physicochemical descriptors, temperature divided
by 298.15 K, and 121 descriptors of the decoded residue-aware UV surface.
Physicochemical inputs share Affinity's mutation feature calculation, including
amino-acid composition, contacts, wild-type geometry and mutation-conditioned atomic
solvent-accessible surface area (SASA). SASA uses the Shrake–Rupley algorithm
with a 1.4 Å probe and 960 sphere points per atom.

The spatial block contains 60 whole-interface descriptors and 61
mutation-conditioned descriptors. These describe partner area and extent,
chemical composition, hydrophobic organization and cross-partner chemical
proximity. Local descriptors use Gaussian kernels with 4 and 12 Å length scales.
Distances follow shortest paths along the surface graph from each residue's
footprint; kernels are truncated at three times their length scale. The ordered
feature names appear in `model.json`.

Residue-aware UV parameterization uses OptCuts. The complete encoded file has a
byte budget of one quarter of the uncompressed native source size. This reference
size includes geometry, residue identities and isolated-partner, bound-complex
and buried SASA. Target cell counts are 1, 1/4 and 1/16 of the native vertex
count; coordinate quantization steps are 0.001, 0.02 and 0.1 Å. Selection
prioritizes cell resolution, then quantization precision, within the byte budget.
Spatial descriptors are calculated after decoding the selected representation.

## Training and prediction

The model is fitted on 4,192 mutation observations from 243 complexes in 35
similarity groups. These observations have paired wild-type and mutant
dissociation constants from SKEMPI 2.0 with exact values and reported
temperatures, together with complete physicochemical and encoded-surface
inputs. The prediction target is ΔΔG = RT ln(Kd_mutant/Kd_wild-type), where
R is the gas constant and T is the reported temperature.

The regressor uses a radial basis function (RBF) kernel with α = 10 and
γ = 1/552. These settings give the lowest mean absolute error among the
evaluated linear ridge and RBF kernel ridge models, averaging errors equally across
similarity groups and then across five inner-validation partitions.
Final fitting uses all eligible observations, with inverse-group-size weights
normalized to mean one. Held-out evaluations use separately fitted models.

Preprocessing uses fitted median imputation, missing-value indicators and
standardization. The bundled files contain these parameters, standardized
training vectors and dual coefficients. Prediction runs on CPU.

## Applicability

Predictions depend on the supplied bound conformation, selected partner groups
and resolved heavy atoms. The calculation retains wild-type geometry throughout
the alanine scan. Native alanines remain in the structural context and are
omitted from the substitution candidates. Missing local surface coverage is
recorded with each result and handled by the fitted preprocessing.

Training coverage varies across interaction families. Use the predicted effects
and rankings to prioritize experimental follow-up.

## Files and license

- `model.json`: model identity, feature order, numerical settings and training scope.
- `weights.npz`: preprocessing arrays, standardized training vectors and dual coefficients.
- `MODEL_CARD.md`: intended use, representation and applicability.

The distributed model parameters and accompanying TopoPPI code use the MIT
license. SKEMPI 2.0 is available at
https://life.bsc.es/pid/skempi2/database/index.
