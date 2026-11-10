# Prediction and interface review

## Install and prepare

Install the core and affinity plugin in Python 3.10 as described in the [README](../README.md). The plugin includes default CPU models for binding affinity and mutation effects. Use `--model` for a custom estimator; neural models require the `models` extra and their declared resources. Inspect current resource availability with `topoppi-affinity models status` or `doctor`. An explicit model catalog can be supplied to `models download --catalog`.

## Run the included example

Export the documented 1BVK antibody–lysozyme example with
`topoppi-affinity example --output affinity-demo`. The complete partner groups
are A,B and C. In the Affinity window, select **Try 1BVK example...**, retain the
default task model and select **Predict**. The `example --run` shortcut performs
an affinity prediction with the bundled CPU model and saves its outputs.

**Local model file** selects a compatible custom model. Models with an explicit
download source can use **Prepare and predict**. **Open saved result...** reads
saved predictions and geometry independently of the model resources.

## Choose a scientific task

The Core **Binding analysis** navigation opens **Binding affinity**, **Mutation
effects**, or **Rank hotspots**. Affinity and mutation effects share a
prediction editor. Its task selector offers whole-complex pKd and
ΔΔG for specified substitutions; model applicability and model source have separate
controls.

The default mutation model supports substitutions on the supplied wild-type (WT)
complex. Advanced **known** and **unseen** ensembles describe training-system applicability.
A local model's declared target and scope are read when it is selected. Changing
tasks retains the variant draft and checks model compatibility before prediction.
Closing and reopening the prediction window also retains the draft and its
model, applicability and thermodynamic settings within the current TopoPPI
session. To continue after exiting the application, choose **Save setup...**
above the setup tabs and save a prediction setup JSON file. **Open setup...**
replaces the editable form with that file. Saving preserves incomplete variant
rows and model settings as entered. Setup files contain paths and settings only;
keep the referenced structure and model files available. Saving a setup is
independent of model loading and the displayed result. Review paths and finish
the draft before predicting.

Selected interface residues can open the shared variant composer. Choose
separate single substitutions for a scan, or combine multiple substitutions
in one variant row. Appending keeps existing rows, adds suffixes to duplicate identifiers and preserves the selected model settings.
If the mutation draft belongs to another structure or partner grouping, save
any rows you need and switch those inputs before appending. An
active prediction must finish or be cancelled before the draft can change.

## Analyze your structure

```sh
topoppi-affinity analyze --structure complex.pdb --chains-a A --chains-b B \
  --task absolute --output results --device cpu
```

The unified command is `topoppi affinity analyze`. `--output-dir` is an alias for `--output`. Add `--json` to receive a JSON summary on standard output; progress messages go to standard error.

For scripted `analyze --layout` runs, exit code 0 confirms prediction success; check `layout_status` in the JSON summary for atlas completion, or use the separate `layout` command and its exit code.

Use `--check-only` to validate inputs before geometry calculation or prediction, `--layout` to request atlases and `--optcuts-bin` to select OptCuts. Each run creates its own result directory; `--retry` reuses a compatible saved interface in a new sibling directory.

## Variant input

A CSV can contain `variant_id,mutations`, with comma-separated substitutions quoted in one field, or structured `variant_id,chain,position,wt,mt` columns and optional insertion codes. A substitution such as `A:G1W` changes the glycine at chain A, residue 1 to tryptophan. Each repeated `--mutations` flag describes one variant; commas join changes in the same variant. Invalid rows retain their identifiers and error messages.

In the prediction editor, CSV import defaults to **Append rows**. Choose
**Replace draft** when the file should become the whole draft. Imported sample
IDs are preserved; resolve conflicting IDs before appending. Ctrl+Z reverses
the most recent import in one step.

For the default mutation model, use `--task mutation` and provide substitutions. Advanced known/unseen ensembles use `--mutation-context known` or `unseen` with the basis recorded in `--scope-note`. An optional measured `--wt-kd` in molar units anchors a mutation prediction to an absolute affinity using `--temperature` in kelvin.

## Review, export and retry

```sh
topoppi-affinity inspect results --json
topoppi-affinity review results --view contacts
topoppi-affinity render results --view surface --output surface.png
topoppi-affinity layout results --optcuts-bin /path/to/OptCuts_bin
```

Review commands use the predictions and geometry stored in saved results. Structure, contact and atlas selections refer to the same author residue identifiers. Layout errors remain separate from prediction results.

In the Affinity window, select **Open saved result...** and choose `analysis.json`. Select a prediction row and binding partner to inspect the 3D interface. Use **Show in main window** to place that selection on the main canvas, **Export predictions...** to share the table, or **Reuse inputs** to start another run. The **Activity** tab and **Show result details** control provide progress messages, units, row errors and saved model context. Failed or missing predictions display **Unavailable** in the table.

## Custom models

Named-feature inference uses NPZ arrays `X`, `sample_ids` and `feature_names`. Select `predict --features`, `--plm-features`, or `--descriptor-features` according to the model schema. A prior CSV identifies rows by `sample_id`. The [model guide](models.md) describes supported fitted objects and the [tools guide](../tools/README.md) lists explicit feature extraction and fitting commands.
