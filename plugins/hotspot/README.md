# TopoPPI Hotspot

TopoPPI Hotspot ranks interface residues within a selected binding partner by predicting the binding free-energy change of each alanine substitution. Its bundled model combines physicochemical mutation descriptors with spatial descriptors recovered from a residue-aware UV representation. Inputs are a structure and the complete chain groups forming the two partners. You can also supply a compatible fitted model.

## Installation

From the repository root, with the core installed:

```bash
python -m pip install -e ./plugins/affinity -e ./plugins/hotspot
```

Hotspot installs the Affinity package for shared mutation features and the Core `representation` extra for surface encoding. The bundled model also requires native OptCuts; follow the [platform installation instructions](../../README.md#optcuts-cannot-be-found). Published pip installations can use `topoppi-install-optcuts` on supported platforms. Set `TOPOPPI_OPTCUTS_BIN` to an existing executable, or pass `--optcuts-bin` to the ranking command. Once OptCuts and the Python dependencies are installed, prediction runs offline on CPU using the bundled weights.

The package uses NumPy, SciPy, pandas, Biopython, and scikit-learn through its dependencies. Optional extras enable PyTorch neural rankers (`neural`), XGBoost (`xgboost`), and ESMC embeddings (`sequence`). Embedding inference requires an explicit local checkpoint.

## Desktop workflow

Start `topoppi-gui`, open **Binding analysis**, and choose **Rank hotspots**. The
Analysis menu and **Hotspot ranking** in **Plugins** also open the workspace directly. The structure and
partner chain groups are filled from the main window; **Use main selection**
refreshes them after a change.

1. Choose the complete partner to rank. Leave **Custom model (optional)** empty to use the bundled model; select a folder containing `model.json` to use your own model.
2. Choose a new results folder. Sequence models also use a prepared embedding registry under **Sequence model options**.
3. Select **Rank residues**. The calculation runs on CPU, displays progress, and can be cancelled.
4. Review the ranked table, select one or more residues, and choose **Show ranked surface** to color the native 3D interface in the main canvas. All selected residues with a native surface footprint are highlighted together; a residue without a footprint remains available in the table for variant creation.
5. Export the full ranking table, or use **Export residue values...** for a `residue,value` annotation CSV that can color an existing Core atlas.
6. Select one or more residues with Ctrl/Cmd-click or Shift-click and choose **Create variants...** to continue in Affinity. Choose separate single substitutions or one joint variant. The selected residues retain their author chain IDs, residue numbers, insertion codes, and wild-type identities. Review the prepared variants, then choose **Predict** in Affinity to run the prediction.

The shared task navigation connects hotspot ranking, binding affinity, and mutation
effects using the current structure and partner input fields. **Create variants...**
continues from the displayed ranking, using its saved structure and partner groups
even when the next ranking's input fields have been edited.
Variant creation requires Affinity to be enabled. Saved rankings can still be
reviewed and exported without it.

For the bundled model, higher predicted ΔΔG values rank first. Positive values
indicate weaker binding after alanine substitution. Predictions use a reference
temperature of 298.15 K. Native alanines remain part of the feature context.
The table shows tied ranks and coverage details where relevant; custom scores
follow the selected model's documented scope.

Compatible geometry from the main workspace is reused. The model's surface and
topology settings determine any additional geometry computation. Desktop runs
save `interface.npz` alongside the ranking files, so **Open saved ranking...**
can recover the table and surface after the original structure or model moves.
Earlier CLI output folders can still be opened as tables. Existing result files
are preserved; **New results folder** creates a fresh destination for another run.
Closing and reopening the Hotspot window in the same TopoPPI session retains
unfinished model, output and rank-partner settings, alongside the displayed
result and selected residues. The next ranking's settings remain separate from
the source and selected partner of the displayed result.

Custom sequence models use local fitted models and prepared embedding registries.
The Python and CLI interfaces also support explicit local ESMC checkpoints.

## Synthetic quick start

```bash
python plugins/hotspot/examples/synthetic_ranking.py
python -m topoppi_hotspot --help
```

The example constructs eight one-dimensional feature rows and two invented queries. It fits a pairwise ranker and prints within-query scores and average precision. Its labels follow the sign of the feature by construction; the output demonstrates the API.

```python
import numpy as np
from topoppi_hotspot.ranking_models import PairwiseHotspotRanker

features = np.array([[-2.0], [-1.0], [1.0], [2.0]])
labels = np.array([0, 0, 1, 1])
queries = np.array(["toy"] * 4)
ranker = PairwiseHotspotRanker(alpha=0.1).fit(features, labels, queries)
scores = ranker.decision_function(features)
```

## Structure ranking

```bash
python -m topoppi_hotspot rank-complex \
  --structure input.pdb --chains-a A --chains-b B --partner A \
  --output-dir ranking_output
```

The unified command is `topoppi hotspot rank-complex`. `--output` is an alias for `--output-dir`. Add `--json` to receive the ranked residues, model information and saved-file locations as JSON on standard output; progress messages go to standard error.

Choose a new output directory for each ranking. The CLI preserves existing ranking files; an empty directory can be reused after a failed calculation.

`--partner` must match a complete input partner group. Comma-separated groups support multichain partners. The optional `--model-dir fitted_model` selects a custom model. Its directory contains `model.json` and the coefficient or estimator files described by that metadata. Supported kinds are `mutation_surface`, `linear`, `mlp`, `xgb`, `graph`, `block_linear`, and `fusion`. See the [official model card](topoppi_hotspot/data/mutation_surface/MODEL_CARD.md) for the bundled model's training source, feature definitions, units and applicability.

Models that use sequence features additionally require `--embedding-registry` or `--esmc-checkpoint`. Fresh embedding inference defaults to CPU; `--device cuda` explicitly selects CUDA. The ESMC implementation uses the 600M architecture and limits observed chains to 1,000 residues.

The bundled model saves:

| File | Content |
| --- | --- |
| `rankings.csv` | Residue identity, alanine substitution, predicted ΔΔG in kcal/mol, ranking score, tied rank interval, within-partner percentile, footprint area, and coverage flags |
| `candidate_coverage.csv` | Structural, interface-context, candidate, and missing-feature counts |
| `metadata.json` | Input locations, feature protocol, model requirements, partner groups, and scoring semantics |
| `features.csv` | The bundled model's 552 inputs for each candidate |
| `representation.bin` | The encoded residue-aware UV representation used by the bundled model |
| `interface.npz` | Native geometry for reopening results and linked surface views |

Feature calculation uses both complete wild-type (WT) partners, including native alanines. Ranked candidates are standard non-alanine residues on the selected partner with a heavy atom within 6 Å of the other partner. The bundled radial basis function (RBF) kernel ridge regressor receives 430 physicochemical mutation descriptors, temperature divided by 298.15 K, and 121 spatial descriptors. Spatial features are computed after decoding a residue-aware UV representation with a byte budget equal to one quarter of the uncompressed native source size. `rank_score` equals `predicted_ddg_kcal_mol` for this model.

Mesh surface areas describe the Gaussian surface; solvent-accessible surface area (SASA) is calculated separately with the Shrake–Rupley algorithm. Custom models use the feature protocol and score meaning recorded in their metadata.

## Python structure ranking

Use the bundled model directly from Python:

```python
from topoppi_hotspot.inference import rank_complex

ranking, metadata = rank_complex(
    "input.pdb", "A", "B", partner="A", output_dir="ranking_output",
)
print(ranking[["residue_key", "predicted_ddg_kcal_mol", "rank_min", "rank_max"]])
```

`ranking` is a pandas DataFrame; `metadata` records the model, units and feature protocol. `output_dir` selects the folder for the saved results. Supply `optcuts_bin` to use a specific OptCuts executable, or `model_dir` to select a custom model.

## Shared Core geometry

Installation registers `hotspot.rank`, the Hotspot desktop workspace, and the `topoppi hotspot` command through
`topoppi.plugins`. Existing `topoppi-hotspot` and `rank_complex` calls remain valid.
The object API reuses Core atoms, native surfaces, SASA input and sequence keys:

```python
import topoppi
from topoppi_hotspot.inference import model_requirements

protocol = model_requirements()["feature_protocol"]
geometry = topoppi.compute_geometry(
    topoppi.get_interface(topoppi.load_structure("input.cif"), "A", "B"),
    surface=topoppi.SurfaceConfig(**protocol["surface"]),
    topology=topoppi.TopologyConfig(**protocol["topology"]),
)
ranking, metadata = topoppi.run(
    "hotspot.rank", geometry, partner="A", output_dir="shared_geometry_ranking",
)
```

For a custom model, pass its directory to both `model_requirements` and `topoppi.run` using `model_dir=`.

The capability accepts `partner`, `model_dir`, `output_dir`, `embedding_registry`,
`esmc_checkpoint`, `optcuts_bin` and `device` from `ExecutionContext.plugin_config["hotspot"]`;
explicit non-None arguments take precedence. Output may also use the context's
`output_dir`. Geometry must contain both surfaces and match the fitted model's
protocol. Incompatible geometry raises an error; prepare geometry with the model's
settings before retrying. Fresh sequence registries use Core-selected atoms for
both PDB and mmCIF input.

## Python modules

| Modules | Purpose |
| --- | --- |
| `features`, `partner_features`, `relative_features` | Residue chemistry, contacts, surface geometry, and interface-relative normalization |
| `mesh_diffusion`, `footprint_diffusion`, `surface_rim`, `paired_surface_features` | Mesh and graph diffusion, rim distances, and cross-partner context |
| `solvent_accessibility` | Isolated-partner, bound-complex, and buried residue SASA |
| `ranking_models`, `block_ranking`, `nonlinear_ranking`, `graph_ranking` | Trainable ranking estimators |
| `ranking_metrics`, `evaluation`, `organization` | Ranking metrics, grouped evaluation, and spatial organization |
| `graph_io`, `graph_operators`, `representations`, `feature_controls` | Array loading, graph alignment, normalization, and feature controls |
| `model_io`, `rank_fusion`, `inference`, `sequence_embeddings` | Model loading, rank aggregation, structure scoring, and explicit checkpoint inference |

The amino-acid tables in `features` use the Kyte-Doolittle hydropathy scale and Zamyatnin residue-volume values. Feature extraction uses these fixed reference descriptors. Pairwise training pairs positive and negative examples within each query and gives queries equal weight.

## Development

Run the independent synthetic checks from the repository root:

```bash
PYTHONPATH=src:plugins/affinity:plugins/hotspot python -m pytest plugins/hotspot/tests -q
```

Tests construct small numerical arrays and use the bundled synthetic two-chain PDB fixture. They check graph alignment, conservation, invariance, train-only scaling, serialization, normalization, and inference dispatch. Optional neural tests skip when PyTorch is unavailable. Contributions should preserve feature names, units, seeds, portable model formats, and meaningful synthetic assertions.

Import normalization from `representations`, partner-relative features from `partner_features`, graph mapping from `graph_operators`, and SASA from `solvent_accessibility`. Existing portable model kinds, metadata keys, feature names, and `rank-complex` arguments remain supported.

The distribution includes the MIT license, the bundled model and synthetic tests.
