# Explore an interface in the desktop app

This guide describes TopoPPI 3.0. Documentation for the version used in the TopoPPI paper is available at [v2.1](https://github.com/GeraltZeroZhong/TopoPPI/tree/v2.1/docs).

TopoPPI's main window is a workspace for viewing protein interfaces. **Setup**
selects the next structure and chain pair. **Appearance** edits the displayed
map. **Advanced** contains calculation and backend settings. The main canvas
stays available while you prepare another input or open an analysis tool.

## Create your first map

1. Choose **Try 1BVK Example...** and select an output folder, or choose your own
   PDB/mmCIF structure. The example selects antibody chain B and lysozyme chain C.
2. Check the surface and partner chains. The surface chain appears on the map;
   the partner chain defines its contact region.
3. Choose **Create Interface Map**. Stage messages show structure preparation,
   interaction detection, surface construction and map layout. **Cancel** stops
   the current task while preserving the previous successful map.
4. Use **2D atlas** and **3D interface** above the canvas to inspect the result.

Chemistry is included and uses ProLIF by default. Native OptCuts supplies the
2D map layout. The [installation guide](../README.md#install-topoppi) describes
the dependencies for each platform.

## Style and save the displayed map

Choose **Appearance** to change residue colors, labels, highlights and numerical
annotations. Footprints show residue regions; markers emphasize selected
residues. The header names the displayed structure, chains and actual view.
If a setting needs correction, an inline message explains the next action and
the existing map remains visible. Correct the setting and choose **Apply Style**.

Load a CSV with `residue,value` columns to color residues by your own measurements
or analysis scores. Choose **Current map** to edit the displayed result or
**Next run** to prepare values for the selected input. See the
[residue footprints guide](residue_footprints.md) for formatting and examples.
The two targets retain separate unit labels and color limits. Preparing the
next dataset leaves the displayed map and its exported units unchanged.

- **Save Figure...** exports the displayed view as an image or vector graphic.
- **Save Atlas...** preserves a map's geometry, labels, colors and annotations
  for later editing.
- **Open Atlas...** restores that saved map directly. This action is available
  in the empty workspace and the main toolbar.

A plugin may provide a native 3D surface before a 2D layout has been computed.
Save its figure or plugin result folder; **Save Atlas** becomes available for a
map with a 2D layout. Affinity and Hotspot result folders retain their native
geometry for further review.

## Open or add a plugin

Choose **Plugins → Manage Plugins...** or the **Plugins...** button beside the
main mapping action. The manager shows each tool's purpose and installation
state. Select an installed tool and choose **Open**. **Chemical interactions**
opens the backend setting for the next map.

Affinity and Hotspot are optional. Choose **Install** to add a package to the
running TopoPPI environment, or **Install from file...** to select a local wheel
or source archive. Installation progress and details appear in the same window.
**Refresh** discovers packages installed while the application is open.

The Affinity plugin includes one official model each for binding affinity and
mutation effects; Hotspot includes one official ranking model. All three run
on CPU. Affinity is ready after plugin installation; Hotspot also requires
native OptCuts. Once these dependencies are installed, prediction runs offline
from local structures. Custom models are available through the advanced model controls.

## Binding analysis

Choose **Binding analysis...** beside the main mapping action or from the
**Analysis** menu. Choose whether to start from the current input fields or the
displayed result. The three tasks share this navigation and carry the structure
and partner groups between workspaces:

- **Binding affinity** predicts whole-complex pKd and the corresponding dissociation constant Kd.
- **Mutation effects** predicts ΔΔG for explicit variants.
- **Rank hotspots** orders residues within one selected partner.

Affinity supplies the first two tasks; Hotspot supplies the third. They remain
separately installable. Missing tasks explain which plugin to install.

### Predict affinity or mutation effects

Open **Binding affinity** or **Mutation effects**. Reuse the structure and complete
partner chain groups, retain the official default model, and enter explicit
substitutions for a mutation task. Review the inputs and results folder, then run
the prediction. Custom models can declare known-system or transfer applicability
and may require additional resources; inspecting the model shows those requirements. The results table links predictions to interface views and
offers CSV export. **Open saved result** restores an earlier analysis for review.
Closing the tool window keeps its unfinished setup during the current application session.
Use **Save setup...** to keep that setup across application restarts;
**Open setup...** restores its inputs, settings and candidate text, including
incomplete rows. Prediction results are saved separately in analysis folders.
The [Affinity workflow guide](../plugins/affinity/docs/workflow.md) explains
whole-complex and mutation-effect tasks.

### Rank interface residues

Open **Rank hotspots** in Binding analysis. Check the partner
chains, choose which partner to rank, and choose a results folder. Leave
**Custom model (optional)** empty to use the official model. Choose **Rank residues**. The table orders residues by model
score; **Show ranked surface** displays the scores on the main canvas. Export
rankings or a residue-value CSV for another map. **Open saved ranking...**
restores the result folder and its native geometry. See the
[Hotspot guide](../plugins/hotspot/README.md) for score interpretation and model
requirements. The default score is predicted alanine-substitution ΔΔG in
kcal/mol at 298.15 K. Higher scores indicate greater predicted loss of binding
and rank first within the selected partner. Native alanines remain in the interface
context and are excluded from the alanine-substitution candidate list.

### Turn selected residues into variants

Select one or more rows in a hotspot table (Ctrl/Cmd-click or Shift-click), then
choose **Create variants...**. On a displayed map, enable **Select residues**:
click residues in 2D or double-click them in 3D. Orange rings mark the selection;
click again to deselect. Dragging in 3D still rotates the structure. Pan and zoom
temporarily take precedence over selection. **Clear** removes the selection.

The variant editor retains each residue's author chain, number, insertion code
and wild-type amino acid. Choose a target amino acid (initially alanine), then
choose **Separate single variants** or **One joint variant**. Inspect the exact
changes before choosing **Add to mutation analysis**, then choose **Predict** in
Mutation effects to run them. Native alanines require another target or removal from the selection.
Scroll horizontally to inspect long joint constructs. **Copy variants** copies
selected preview rows, or all rows when none are selected, as a tab-separated
table with variant IDs and exact changes. Ctrl/Cmd+C works in the preview.

Variants use the displayed result's original structure and partner groups even
if the next-input fields have changed. Existing variants remain in the Affinity
setup; a conflicting structure or partner grouping requires explicitly changing
that setup first. A saved map without its original PDB/mmCIF remains viewable,
and **Locate structure...** in the variant editor lets you supply its matching
source. The editor checks recorded partner chains and selected residue identities
before continuing. These checks cover residue identities; coordinate correspondence
depends on the file you select. The saved atlas retains its original provenance.
Selections survive style/view changes and clear when a new result is displayed.

## Work with a smaller window or keyboard

Resize the sidebar to balance controls and the canvas. Settings scroll as
keyboard focus moves through them. **Tab** moves between controls and
**Page Up / Page Down** scroll settings. Use **Ctrl+O / Ctrl+S** on Windows and
Linux, or **Command+O / Command+S** on macOS, to open a structure and save a figure.
**Show activity** expands detailed progress when you need it.
