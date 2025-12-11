import argparse
import sys
import os
import time
import logging
import subprocess
import gemmi

# Ensure src is in python path
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

from src.io_loader import PDBLoader
from src.surface import SurfaceGenerator
from src.topology import TopologyManager
from src.parameterization import Parameterizer
from src.visualizer import InterfaceVisualizer

def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%H:%M:%S'
    )

def generate_arpeggio_interactions(pdb_path, logger):
    """
    Generates an Arpeggio-compatible mmCIF file from a PDB and runs pdbe-arpeggio.
    Matches logic in gui.py.
    """
    logger.info("Checking Arpeggio requirements...")
    
    pdb_dir = os.path.dirname(os.path.abspath(pdb_path))
    pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]
    
    # Expected output JSON
    json_name = f"{pdb_name}.json"
    expected_json = os.path.join(pdb_dir, json_name)
    
    if os.path.exists(expected_json):
        logger.info(f"Found existing Arpeggio output: {json_name}")
        return expected_json

    # Target CIF file for Arpeggio input
    cif_path = os.path.join(pdb_dir, f"{pdb_name}.cif")
    target_file = cif_path

    # --- Conversion / Injection Logic ---
    try:
        # Always process the file using Gemmi to ensure _chem_comp exists
        logger.info("Preparing Arpeggio-compatible mmCIF...")
        
        # 1. Read Structure (handles PDB or CIF)
        if pdb_path.lower().endswith('.cif'):
            doc = gemmi.cif.read(pdb_path)
            block = doc.sole_block()
            # Check if we need to add chem_comp
            if not block.find_loop("_chem_comp.id"):
                structure = gemmi.read_structure(pdb_path) 
            else:
                target_file = pdb_path # Use original if valid
                structure = None
        else:
            structure = gemmi.read_structure(pdb_path)

        # 2. Re-write as MMCIF with _chem_comp if we have a structure object
        if structure:
            doc = structure.make_mmcif_document()
            block = doc.sole_block()

            # Gather unique residue names present in the structure
            res_names = set()
            for model in structure:
                for chain in model:
                    for res in chain:
                        res_names.add(res.name)

            # Define standard amino acids (Peptide linking)
            STANDARD_AAS = {
                'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 
                'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 
                'TYR', 'VAL', 'UNK'
            }

            # 3. Inject the _chem_comp loop
            loop = block.init_loop('_chem_comp.', ['id', 'type', 'name'])
            
            for r in sorted(res_names):
                if r in STANDARD_AAS:
                    ctype = 'L-peptide linking'
                    name = 'AMINO ACID'
                elif r == 'HOH':
                    ctype = 'water'
                    name = 'WATER'
                else:
                    ctype = 'non-polymer'
                    name = 'LIGAND'
                
                # Explicitly quote strings containing spaces
                row = [
                    gemmi.cif.quote(r),
                    gemmi.cif.quote(ctype),
                    gemmi.cif.quote(name)
                ]
                loop.add_row(row)

            doc.write_file(cif_path)
            target_file = cif_path
            logger.info("Created Arpeggio-compatible CIF with _chem_comp dictionary.")

    except Exception as e:
        logger.warning(f"CIF Preparation failed: {e}")
        return None

    # --- Run Arpeggio ---
    cmd = ["pdbe-arpeggio", target_file]
    
    try:
        logger.info(f"Running: pdbe-arpeggio on {os.path.basename(target_file)}...")
        # Capture output to avoid cluttering CLI unless verbose
        subprocess.run(cmd, check=True, cwd=pdb_dir, capture_output=True)
        
        if os.path.exists(expected_json):
            logger.info("Arpeggio JSON generated successfully.")
            return expected_json
        else:
            logger.warning("Arpeggio ran but output JSON not found.")
            return None
            
    except subprocess.CalledProcessError as e:
        logger.warning(f"Arpeggio generation failed: {e}")
        return None
    except FileNotFoundError:
        logger.warning("Error: 'pdbe-arpeggio' command not found. Skipping interaction analysis.")
        return None

def main():
    parser = argparse.ArgumentParser(description="ProtSurf2D: Protein Interface 2D Map Generator")
    parser.add_argument("pdb_file", help="Path to input PDB/CIF file")
    parser.add_argument("-A", "--chain_a", required=True, help="Chain ID for the surface (Receptor)")
    parser.add_argument("-B", "--chain_b", required=True, help="Chain ID for the ligand")
    parser.add_argument("--arpeggio", help="Path to Arpeggio JSON interactions file (Optional)", default=None)
    parser.add_argument("--cutoff", type=float, default=5.0, help="Interface distance cutoff (Angstroms)")
    parser.add_argument("--res", type=float, default=1.0, help="Grid resolution for surface generation (Angstroms)")
    parser.add_argument("--sigma", type=float, default=1.5, help="Gaussian smoothing sigma")
    parser.add_argument("--output", "-o", default="interface_map.png", help="Output image filename")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose debug logging")
    
    args = parser.parse_args()
    setup_logging(args.verbose)
    logger = logging.getLogger("Main")

    start_time = time.time()

    # --- Step 0: Arpeggio Handling ---
    arpeggio_file = args.arpeggio
    
    # If path is empty or invalid, try to generate it
    if not arpeggio_file or not os.path.exists(arpeggio_file):
        generated_json = generate_arpeggio_interactions(args.pdb_file, logger)
        if generated_json:
            arpeggio_file = generated_json
        else:
            arpeggio_file = None # Explicitly None to trigger heuristic fallback in Visualizer

    # --- Step 1: Data Loading ---
    logger.info(f"Loading {args.pdb_file}...")
    try:
        loader = PDBLoader(args.pdb_file)
        coords_A, atoms_A = loader.get_chain_data(args.chain_a)
        coords_B, atoms_B = loader.get_chain_data(args.chain_b)
        logger.info(f"Loaded Chain {args.chain_a}: {len(coords_A)} atoms")
        logger.info(f"Loaded Chain {args.chain_b}: {len(coords_B)} atoms")
    except Exception as e:
        logger.error(f"Failed to load PDB data: {e}")
        sys.exit(1)

    # --- Step 2: Surface Generation ---
    logger.info("Generating SES Surface for Chain A...")
    surf_gen = SurfaceGenerator(coords_A)
    mesh_A = surf_gen.generate_mesh(grid_resolution=args.res, sigma=args.sigma)
    
    if mesh_A is None or len(mesh_A.vertices) == 0:
        logger.error("Failed to generate surface mesh.")
        sys.exit(1)

    # --- Step 3: Topology Extraction ---
    logger.info("Extracting interface patches...")
    topo_mgr = TopologyManager(mesh_A, coords_B)
    patches = topo_mgr.get_interface_patches(distance_cutoff=args.cutoff)
    
    if not patches:
        logger.error("No interface patches found. Try increasing --cutoff.")
        sys.exit(0)

    # --- Step 4: Parameterization (LSCM) ---
    logger.info(f"Parameterizing {len(patches)} patches...")
    valid_patches = []
    param = Parameterizer()
    
    for i, patch in enumerate(patches):
        logger.info(f"  Flattening Patch {i+1} ({len(patch.vertices)} verts)...")
        uv = param.flatten_patch(patch)
        
        if uv is not None:
            # Attach UV data to the mesh object for the visualizer
            patch.metadata['uv'] = uv
            valid_patches.append(patch)
        else:
            logger.warning(f"  Skipping Patch {i+1} due to parameterization failure.")

    if not valid_patches:
        logger.error("All patches failed to parameterize.")
        sys.exit(1)

    # --- Step 5: Visualization ---
    logger.info("Visualizing results...")
    
    # Updated signature to support Arpeggio
    viz = InterfaceVisualizer(
        chain_A_atoms=atoms_A, 
        chain_A_coords=coords_A, 
        chain_B_coords=coords_B, 
        chain_B_atoms=atoms_B,
        chain_a_id=args.chain_a,
        chain_b_id=args.chain_b,
        arpeggio_file=arpeggio_file
    )
    
    viz.plot_patches(valid_patches, output_file=args.output)

    elapsed = time.time() - start_time
    logger.info(f"Done! Pipeline finished in {elapsed:.2f}s. Saved to {args.output}")

if __name__ == "__main__":
    main()
