# Residue footprints

TopoPPI draws filled residue regions on the calculated UV atlas and, from version 2.1, on the native 3D interface surface. Each triangle is partitioned into three equal-area quadrilaterals using its corners, edge midpoints and barycentre. Each quadrilateral inherits its corner's source residue. The 2D atlas includes separated residue pieces, original surface boundaries and both occurrences of each optimization seam. The [3D interface view](3d_interface.md) shows these regions on the curved surface with the optimized seams in their structural locations.

Changing colors, labels or external values reuses the UV coordinates. Footprint mode displays every retained optimized patch, including patches below the marker view's interaction-count threshold. `residue_scope` controls labels, highlights and numerical coloring while the complete region geometry stays visible.

## Create a map

```bash
topoppi complex.cif -A A -B B \
  --map-style footprints \
  --highlight A:GLU:37 A:TYR:40 \
  --export-atlas interface.atlas.npz \
  -o interface.svg
```

The CLI selects full patch annotation scope when `--map-style footprints` is requested, unless `--residue-scope interaction` is supplied. Pale blue denotes ordinary regions and magenta denotes selected residues. Labels identify the surface chain's author residue numbers, including insertion codes.

The standard interaction source remains ProLIF. For a run using distance-defined optimization weights, add `--interaction-source geometric`. This uses the distinct partner residues within `--geometric-fallback-distance` (default 6 Å). The surface-domain cutoff is a separate setting, `--cutoff` (default 4 Å). External annotations and highlights are applied during rendering after optimization.

## Add numerical annotations

Create a CSV with `residue,value` columns:

```csv
residue,value
A:GLU:37,0.98
A:TYR:40,2.32
A:ILE:24,NA
```

These rows illustrate the file format. Supply the values and units for your measurement and save the file as UTF-8 CSV; files exported with a UTF-8 byte-order mark are accepted. Full residue keys have the form `chain:three-letter-name:author-number`. Unambiguous `A:37` and `37` aliases are also accepted. Empty values and `NA`, `N/A` or `NaN` denote unavailable measurements. Missing interface values appear grey. Duplicate or unknown keys are reported with their residue identity. A whole-source-chain table is accepted; the report counts entries outside the mapped interface.

```bash
topoppi render interface.atlas.npz \
  --annotation-file effects.csv \
  --annotation-label 'Effect (kcal/mol)' \
  --vmin -2.5 --vmax 2.5 \
  --export-atlas annotated.atlas.npz \
  -o interface_effects.pdf
```

The default blue–white–brown value scale is symmetric about zero and spans the largest absolute interface value. Explicit lower and upper limits can be used to match several figures. Ranges crossing zero retain white at zero; a range entirely above or below zero uses a linear scale. All pieces of a residue receive its value. Numeric coloring takes precedence over manual region colors; clearing annotations restores the saved manual colors. Highlights still prioritize labels in a numeric view.

Arrowheads on the colorbar indicate displayed values below or above the selected range. Those residues use the endpoint color and retain their original values in the saved atlas. The rendering report gives `below_scale_residue_count`, `above_scale_residue_count` and `colorbar_extend` for the active annotation scope.

## Edit in the desktop app

1. Under **Map Display**, select **Residue footprints** under **Map style**. Choose **2D atlas** or **3D interface** under **View**. Footprints select **Full patch context** by default. **Residue scope** in Advanced settings can narrow annotations to interaction residues.
2. Enter selected residues in **Highlight**, separated by commas or spaces. Choose all, highlighted, or no footprint labels.
3. Use the border and seam controls and the region, highlight and missing-value colors to style the map.
4. Load a CSV to color by value. Use **Current map** to edit the displayed atlas or **Next run** to prepare values for the next selected structure. Set the colorbar label and optional numerical limits; **Clear** removes the selected annotation layer.
5. Apply the style to redraw the current atlas. In 2D, click a region to recolor that residue across its pieces, or drag a label. In 3D, drag to rotate the surface and double-click a visible residue to recolor it. When numerical values are active, clear them before assigning manual region colors.
6. Use **Save Figure** for an image or vector file. Use **Save Atlas** to preserve the editable result, then **Open Atlas** to continue in a later session.

Mode changes and style updates use the cached optimization result. The displayed map identifies its structure and chains, so edits remain associated with the current atlas as another input is prepared. A loaded atlas includes the atom identities, resolved interaction data, UV geometry, external values and plotting style needed for offline editing.

When computation completes and a display setting needs correction, the GUI retains the calculated atlas. Correct the setting and choose **Apply Style**, or use **Save Atlas** to continue editing later.

## Render a saved atlas

```bash
topoppi render annotated.atlas.npz -o another_copy.svg

topoppi render annotated.atlas.npz --clear-annotations \
  --highlight A:37 A:40 --labels highlighted \
  --highlight-color '#A64D79' -o selected_residues.png
```

Rendering reads the geometry, interactions and embedded annotations from the saved atlas. The original structure, input JSON/CSV files and native solver can be stored separately. `--export-atlas` saves the updated style into another atlas file. The compressed NPZ contains numeric/string arrays and JSON.

Saved atlases keep all optimized patches. The 2D marker view applies its interaction-count threshold; footprints and the 3D interface view display every retained patch. View changes preserve residue annotations, and saved 3D cameras can be reopened or reused by `topoppi render`.

PNG output is 300 dpi. TIFF uses 600 dpi and LZW compression. SVG keeps text editable and PDF embeds TrueType fonts. All four formats share the same renderer. Atlas width is 178 mm; height follows the mapped geometry, with a compact fixed-height footer when a value scale is present.

## Options

| Option | Purpose |
| --- | --- |
| `--map-style markers\|footprints` | Choose the residue representation |
| `--view atlas\|surface` | Choose the 2D atlas or 3D interface view |
| `--highlight RESIDUE ...` | Highlight selected residues; commas are also accepted |
| `--labels all\|highlighted\|none` | Choose footprint label coverage |
| `--annotation-file CSV` | Read external residue values |
| `--annotation-label TEXT` | Label the value scale and its units |
| `--vmin NUMBER`, `--vmax NUMBER` | Set common value limits |
| `--footprint-color COLOR` | Set the neutral region color |
| `--highlight-color COLOR` | Set the selected-residue color |
| `--missing-color COLOR` | Set unavailable-value color |
| `--hide-residue-borders`, `--show-residue-borders` | Hide or show internal residue borders |
| `--hide-seams`, `--show-seams` | Hide or show optimization seam lines |
| `--export-atlas FILE.npz` | Save complete geometry and editable annotations |
| `render --clear-annotations` | Remove embedded external values for this rendering |

Original patch boundaries remain visible. Label anchors lie inside actual regions. When dense labels cannot be placed without overlap, their regions remain visible and the rendering report lists omitted labels. Explicitly dragged labels retain their offsets across redraws and saving.

## Python API

```python
from topoppi.config import TopoPPIRunConfig, VisualizationConfig
from topoppi.pipeline import run_interface_mapping
from topoppi.visualization.atlas_io import load_atlas, save_atlas

result = run_interface_mapping(TopoPPIRunConfig(
    pdb_file="complex.cif", chain_a="A", chain_b="B",
    output_file="interface.pdf", atlas_output="interface.atlas.npz",
    interaction_source="geometric",
    visualization=VisualizationConfig(
        map_style="footprints", residue_scope="patch",
        highlight_residues=("A:GLU:37", "A:TYR:40"),
    ),
))

atlas = load_atlas("interface.atlas.npz")
style = {**atlas.style, "highlight_color": "#7A4F8C"}
figure = atlas.visualizer.plot_patches(
    atlas.patches, style_config=style, output_file="restyled.svg", show=False,
)
save_atlas("restyled.atlas.npz", atlas.patches, atlas.visualizer,
           run_metadata=atlas.metadata)
```

Use a saved atlas to reproduce a particular geometry exactly. A fresh calculation additionally depends on the input coordinates, chain direction, surface settings, interaction-weight source, solver version and optimization settings. Two independent interfaces have independently optimized coordinates; corresponding residues are identified by their source identities.
