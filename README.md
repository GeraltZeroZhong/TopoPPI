# TopoPPI

## Overview
TopoPPI maps protein–protein interfaces onto 2D surfaces so you can quickly inspect the residue contacts that drive binding. The pipeline loads a PDB/CIF structure, builds a solvent-excluded surface for a selected receptor chain, finds interface patches against a ligand chain, flattens those patches with least-squares conformal mapping (LSCM) or Harmonic Parameterization, and renders publication-ready interaction maps. A Tkinter GUI and a CLI are both available for interactive exploration or scripted runs.

## Features
- **Multi-format input**: Parse PDB or mmCIF structures via Biopython.
- **Surface generation**: Convert receptor atoms into a Gaussian-smoothed solvent-excluded surface mesh.
- **Interface detection**: Identify surface patches within a configurable distance cutoff to the ligand chain.
- **2D parameterization**: Flatten interface patches using LSCM for distortion-minimized maps. If the singular matrix proves intractable to remedy, use Harmonic Parameterization instead.
- **Visualization**: Export interaction maps that can incorporate Arpeggio interaction annotations.
- **Two front ends**: Run via command line (`main.py`) or a Tkinter desktop app (`gui.py`).

## Installation
1. Install [conda](https://docs.conda.io/en/latest/miniconda.html).
2. Create the environment:
   ```bash
   conda env create -f environment.yml
   conda activate bio3d
   ```
3. Ensure that `pdbe-arpeggio` is available (installed automatically via `pip` in the environment file) when you want Arpeggio interaction overlays.

## Usage
### Command line
Run the full pipeline on a structure and export a PNG interface map:
```bash
python main.py <structure.pdb> --chain_a A --chain_b B --output interface_map.png \
  --cutoff 5.0 --res 1.0 --sigma 1.5 [--arpeggio path/to/arpeggio.json] [--verbose]
```
Key arguments:
- `--chain_a` / `--chain_b`: Receptor and ligand chain IDs.
- `--cutoff`: Distance cutoff (Å) for defining interface residues.
- `--res`: Grid resolution (Å) for surface voxelization.
- `--sigma`: Gaussian smoothing sigma for the density field.
- `--arpeggio`: Optional precomputed Arpeggio JSON file; if omitted, the GUI can generate one automatically.
- `--output`: Destination image file for the interface map.

### Graphical UI
Launch the desktop application:
```bash
python gui.py
```
Steps:
1. Choose a PDB/CIF file and (optionally) an Arpeggio JSON file.
2. Enter receptor/ligand chain IDs and adjust cutoff, grid resolution, and smoothing sigma.
3. Click **Run Full Analysis** to generate and visualize interface patches. Use the style controls to recolor interactions or filter interaction types.

## Configuration
- **Interface cutoff** (`--cutoff` or GUI field): Controls which surface vertices are considered part of the interface.
- **Grid resolution** (`--res`): Smaller values capture more surface detail at the cost of runtime and memory.
- **Gaussian sigma** (`--sigma`): Higher values smooth the surface more aggressively.
- **Arpeggio interactions**: Provide an Arpeggio JSON file to color residues by interaction type; otherwise, only geometric proximity is used.
- **Logging**: Pass `--verbose` to the CLI for debug-level logging.

## Folder Structure
- `main.py` — CLI entry point for running the full analysis pipeline.
- `gui.py` — Tkinter desktop application for interactive analysis.
- `src/` — Core library code:
  - `io_loader.py` — PDB/mmCIF loading and chain extraction.
  - `surface.py` — Solvent-excluded surface generation from atomic coordinates.
  - `topology.py` — Interface patch detection and topology utilities.
  - `parameterization.py` — Patch flattening via LSCM and Harmonic Parameterization.
  - `visualizer.py` — Plotting routines for interface maps.
- `data/getHuman2ChainsPDB.py` — Downloads Homo sapiens dimers with two chains and normalizes chain IDs to A/B.
- `sensitivity_analysis.py` — Runs cutoff/sigma sensitivity sweeps and writes interface-area summaries to CSV.
- `environment.yml` — Conda environment specification with Python and scientific stack dependencies.
- `LICENSE` — MIT License details.

## License
This project is licensed under the terms of the MIT License. See the [LICENSE](LICENSE) file for details.
