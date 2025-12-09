ProtSurf2D

Overview

ProtSurf2D is a toolkit for extracting, flattening, and visualizing protein–protein interface patches as 2D maps. It loads PDB/CIF structures, generates a solvent-excluded surface for a target chain, detects interface regions, parameterizes each patch to 2D using LSCM, and visualizes residue interactions using Arpeggio or geometric heuristics.

Features

Accurate Interactions: Integrated support for PDBe-Arpeggio to detect detailed interactions (salt bridges, cation-pi, pi-stacking, etc.).

Interface-aware surface generation: Builds smoothed solvent-excluded surfaces using Gaussian density fields.

2D parameterization: Flattens curved protein surfaces into 2D maps using LSCM.

Rich visualization: Interactive GUI to explore contacts and export publication-quality figures.

Installation

Single Environment Setup

Create the unified environment (Python 3.10) which includes all dependencies (GUI, Surface algorithms, and Arpeggio):

conda env create -f environment.yml
conda activate bio3d


Note: This installs openbabel via conda-forge, which is required for Arpeggio to work.

Usage

Graphical Interface (Recommended)

Launch the GUI:

python gui.py


Click Browse PDB... to select your structure.

(Optional) Provide an existing Arpeggio JSON file. If left empty, ProtSurf2D will automatically run Arpeggio in the background.

Click Run Full Analysis.

Command Line

Run the pipeline manually:

python main.py structure.pdb \
  --chain_a A \
  --chain_b B \
  --arpeggio structure.json \
  --output result.png


Configuration

Interaction Types: Supports salt_bridge, cation_pi, pi_stack, hbond, sulfur_pi, hydrophobic, halogen_bond, metal_complex.

Parameters: Adjustable cutoff distance, surface resolution, and smoothing sigma.

License

MIT License
