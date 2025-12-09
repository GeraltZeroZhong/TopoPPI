# ProtSurf2D

## Overview
ProtSurf2D is a toolkit for extracting, flattening, and visualizing protein–protein interface patches as 2D maps. It loads PDB/CIF structures, generates a solvent-excluded surface for a target chain, detects interface regions against a binding partner, parameterizes each patch to 2D with Least Squares Conformal Mapping (LSCM), and renders publication-ready figures via CLI or GUI workflows.

## Features
- **Interface-aware surface generation**: Builds smoothed solvent-excluded surfaces from atomic coordinates using Gaussian density fields and marching cubes meshing.
- **Patch extraction**: Identifies interacting surface regions between two chains based on a configurable distance cutoff.
- **2D parameterization**: Flattens each patch with LSCM so contacts can be inspected and annotated on a planar map.
- **Rich visualization**: Exports PNG images from the CLI or interactively explores styles (interaction coloring, fonts, custom colors) in the Tkinter GUI.
- **Biopython-powered parsing**: Reads both PDB and mmCIF inputs while filtering heteroatoms and water.

## Installation
1. Install [conda](https://docs.conda.io/en/latest/miniconda.html).
2. Create the environment:
   ```bash
   conda env create -f environment.yml
   conda activate bio3d
   ```
3. If you prefer pip, ensure Python 3.10+ and install the dependencies listed in `environment.yml` (numpy, scipy, biopython, scikit-image, matplotlib, tk, igl, trimesh, networkx, pillow, rtree, shapely).

## Usage
### Command line
Run the full pipeline and export a figure:
```bash
python main.py <structure.pdb> \
  --chain_a A \
  --chain_b B \
  --cutoff 5.0 \
  --res 1.0 \
  --sigma 1.5 \
  --output interface_map.png \
  --verbose
```
- `--chain_a` is the receptor chain for surface generation.
- `--chain_b` is the ligand chain used to find contacts.
- `--cutoff` controls the interaction distance (Å).
- `--res` and `--sigma` tune surface detail and smoothing.

### Graphical interface
Launch the Tkinter GUI to browse a structure, configure parameters, and adjust visualization styles interactively:
```bash
python gui.py
```
Use **Run Full Analysis** to execute the pipeline, then **Update Style Only** to recolor by interaction type or custom residue color.

## Configuration
Key runtime parameters are passed via CLI flags or GUI fields:
- `--cutoff`: distance threshold for interface detection (Å).
- `--res`: surface grid resolution (Å) affecting mesh detail.
- `--sigma`: Gaussian smoothing sigma for surface density.
- `--output`: filename for the generated image.
- GUI style options: color by interaction type, enable/disable interaction categories, font family/size, custom residue color, and pick-to-recolor on the plot.

## Folder Structure
```
ProtSurf2D/
├── main.py            # CLI pipeline entrypoint
├── gui.py             # Tkinter GUI for interactive analysis
├── src/
│   ├── io_loader.py   # PDB/mmCIF parsing and chain extraction
│   ├── surface.py     # Solvent-excluded surface generation
│   ├── topology.py    # Interface patch detection
│   ├── parameterization.py  # LSCM flattening of surface patches
│   └── visualizer.py  # Plotting and styling for interface maps
├── environment.yml    # Conda environment with runtime dependencies
└── LICENSE            # MIT License
```

## License
This project is distributed under the terms of the MIT License. See `LICENSE` for details.
