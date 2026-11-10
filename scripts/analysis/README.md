# Representation analysis scripts

These scripts use `topoppi_representation` and accept
user-selected inputs and output directories. Install TopoPPI 3.0 with its
representation dependencies first. Run each script with `--help` for parameters.

| Script | Purpose | Additional dependencies |
| --- | --- | --- |
| `compare_representations.py` | Compare one encoding with its native source and measure CPU cost | Representation extra |
| `physical_recovery.py` | Fit ridge regression or kernel ridge regression with a radial basis function (RBF) kernel for physical targets | pandas, scikit-learn, joblib |
| `finite_readout.py` | Fit an area-weighted graph readout for physical targets | pandas, PyTorch |

## Representation and resource comparison

```bash
python scripts/analysis/compare_representations.py interface_source.npz \
  --codec uv --layout ordinary --cells 1024 --quant-step 0.001 \
  --output comparison
```

The input is produced by `topoppi_representation.save_source`; keep its companion
JSON file. Other codec selections are `native` and `direct`. The native options
are `--fraction` and `--bits`. UV inputs need the selected layout on every patch.

Outputs are `encoded.bin` and `comparison.json`. Measurements include complete
encoded bytes, canonical source bytes, and wall/process CPU times for encoding,
decoding and physical queries. `python_traced_peak_bytes` is the peak memory
tracked by Python's allocator instrumentation. Total process resident set size
(RSS) and GPU memory require separate measurements. Timing starts after source
loading and initial reference-query construction.

The direct table stores the two hydrophobic-organization fractions and eleven
chemical-proximity components. Its recovery error applies to those named queries.
Native and UV comparisons report organization and proximity errors along with
carried-area agreement. Native proximity uses the requested triangle quadrature;
UV proximity uses its decoded cell positions.

## Physical recovery

```bash
python scripts/analysis/physical_recovery.py features.npz targets.csv \
  --kind ridge --alpha 1 --output recovery
```

The feature NPZ contains `X` (samples by features), `sample_ids`, and
`feature_names`. The target CSV contains:

- `sample_id`: the same unique identifiers as the feature NPZ;
- `group`: the identifier used to keep related samples in the same partition;
- `split`: `train`, `validation`, or `test`;
- the thirteen columns in `analysis.graphs.TARGET_NAMES`.

Those columns are `O_side0`, `O_side1`, and `X_PP`, `X_AA`, `X_CC`, `X_PA`,
`X_PC`, `X_AC`, `X_PX`, `X_AX`, `X_CX`, `X_XX`, `X_none`. Missing target values
are represented by empty CSV cells. Chemical-proximity distributions are supplied
as complete eleven-component targets. Groups must remain within one partition.

The script fits the selected configuration on training rows, including imputation
and standardization, and reports predictions on test rows. Validation rows remain
available for configuration selection outside this command. Outputs are
`reader.joblib`, `predictions.csv`, and `report.json`.
Load joblib models from trusted sources.

## Finite readout

```bash
python scripts/analysis/finite_readout.py graphs.csv \
  --view native --depth 1 --width 32 --epochs 200 --seed 17 \
  --output graph_reader
```

Build samples with `analysis.graphs.make_probe_sample` from a `ComplexSource` and
explicit physical targets, then write each sample using `save_probe_sample`.
The manifest CSV contains `path` and `split`; relative paths resolve beside the
manifest. Splits are `train`, `validation`, and `test`. Sample groups come from
the sample's `family_id` and stay disjoint across partitions.

The available views are self-only, native adjacency, ordinary UV adjacency,
contact-weighted residue-aware UV adjacency and uniformly weighted residue-aware UV
adjacency. View names are `self`, `native`, `ordinary`, `aware`, and
`aware_uniform`. Uniform views use the layout metadata emitted by
`layouts.augment_uniform_layout`.

All views use the same native nodes, feature channels and surface-area-weighted
pooling. Early stopping uses validation loss. Test rows are used for prediction after
the selected checkpoint is restored. Outputs are `reader.pt`,
`predictions.csv`, and `fit.json`. CPU is the default device; GPU use requires an
explicit `--device` selection.

