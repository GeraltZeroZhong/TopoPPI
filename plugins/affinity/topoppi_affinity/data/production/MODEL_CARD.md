# Default binding models, version 3.0

TopoPPI Affinity includes two CPU inference models. Each model
comprises a JSON feature and training specification and an NPZ containing median
imputation values, missingness-indicator positions and scaling parameters. The
linear ridge affinity model stores regression coefficients and an intercept.
The mutation model uses kernel ridge regression with a radial basis function (RBF)
kernel and stores standardized training vectors and dual coefficients.
Inference uses NumPy and calculates all requested
physical inputs from the supplied structure.

| Model | Target | Inputs | Training source |
| --- | --- | --- | --- |
| Absolute affinity | pKd = −log10(Kd/M) | 477 exposure, burial, chemical-composition and atomic-contact descriptors | PDBbind 2020 release 1 exact dissociation constants with BIORAD curation |
| Mutation effect | ΔΔG in kcal/mol; positive weakens binding | 263 native-interface and substitution descriptors; 167 WT SASA descriptors; temperature/(298.15 K) | SKEMPI 2.0 paired dissociation constants at reported temperatures |

Model selection compares six candidates: linear ridge regression with
regularization parameter α ∈ {0.1, 10, 1000} and RBF kernel ridge regression with
α ∈ {0.1, 1, 10}. The RBF kernel uses γ = 1/d, where d is the number of input
features before adding missingness indicators. Selection minimizes group-macro
mean absolute error (MAE), averaged across five inner-validation partitions.
The affinity model uses linear ridge regression with α = 1000; the mutation
model uses the RBF kernel with α = 1 and γ = 1/431.
Each selected model is refitted on all eligible training samples with inverse
similarity-group weights normalized to mean one. Imputation and scaling are fitted with the model.
The JSON records the training-set counts, parameters and exact feature ordering.
Generalization performance is estimated from separate held-out fits.

Use prepared protein complexes and complete interacting partners. Mutation
predictions retain wild-type (WT) coordinates and support single and multiple
substitutions. Solvent-accessible surface area (SASA) is computed on atoms with a
1.4 Å probe and 960 sampling points per atom. Geometric descriptors use the
Gaussian molecular surface. User-selected chain names and author residue numbering, including insertion codes, identify
substitutions. Structural views show input contacts and residue locations.
Individual predictions have no calibrated confidence interval. Their accuracy
depends on the supplied conformation and similarity to training systems.

Sources: [PDBbind](https://www.pdbbind-plus.org.cn/download),
[BIORAD](https://github.com/lanl/BIORAD/tree/main/data),
[SKEMPI 2.0](https://life.bsc.es/pid/skempi2/database/index).
