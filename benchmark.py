import os
import sys
import numpy as np
import trimesh
import csv
import time
import logging

# Ensure src is in python path
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

from src.io_loader import PDBLoader
from src.surface import SurfaceGenerator
from src.topology import TopologyManager
from src.parameterization import Parameterizer

# Setup logging
logging.basicConfig(level=logging.INFO, format='[Benchmark] %(message)s')
logger = logging.getLogger("Benchmark")

class BenchmarkAlgorithms:
    """
    Implements comparison algorithms (Planar, Spherical, Cylindrical)
    to serve as baselines for benchmarking TopoPPI (LSCM).
    """
    
    @staticmethod
    def pca_align(mesh):
        """
        Helper: Aligns the mesh's principal axes to the Cartesian axes using PCA.
        - X axis: Smallest variation (Normal)
        - Y axis: Middle variation
        - Z axis: Largest variation (Longest axis)
        """
        if len(mesh.vertices) < 3: return mesh
        
        points = mesh.vertices - mesh.centroid
        cov = np.dot(points.T, points) / points.shape[0]
        evals, evecs = np.linalg.eigh(cov)
        
        sort_idx = np.argsort(evals)
        
        rotation = np.eye(4)
        rotation[:3, 0] = evecs[:, sort_idx[0]] # X
        rotation[:3, 1] = evecs[:, sort_idx[1]] # Y
        rotation[:3, 2] = evecs[:, sort_idx[2]] # Z
        
        transform = rotation.copy()
        transform[:3, :3] = rotation[:3, :3].T
        
        aligned_mesh = mesh.copy()
        try:
            aligned_mesh.apply_transform(transform)
        except:
            pass # Ignore transform errors on degenerate meshes
        return aligned_mesh

    @staticmethod
    def planar_projection(mesh):
        aligned = BenchmarkAlgorithms.pca_align(mesh)
        v = aligned.vertices - aligned.centroid
        
        # Project onto YZ plane (drop X)
        y = v[:, 1]
        z = v[:, 2]
        
        def normalize(arr):
            rng = arr.max() - arr.min()
            if rng < 1e-9: return np.zeros_like(arr)
            return (arr - arr.min()) / rng

        u = normalize(y)
        v_coord = normalize(z)
        
        return np.column_stack([u, v_coord])

    @staticmethod
    def spherical_projection(mesh):
        aligned = BenchmarkAlgorithms.pca_align(mesh)
        v = aligned.vertices - aligned.centroid
        
        r = np.linalg.norm(v, axis=1)
        r[r == 0] = 1e-9 
        
        theta = np.arctan2(v[:, 1], v[:, 0])
        u = (theta + np.pi) / (2 * np.pi)
        
        phi = np.arccos(v[:, 2] / r)
        v_coord = phi / np.pi
        
        return np.column_stack([u, v_coord])

    @staticmethod
    def cylindrical_projection(mesh):
        aligned = BenchmarkAlgorithms.pca_align(mesh)
        v = aligned.vertices - aligned.centroid
        
        theta = np.arctan2(v[:, 1], v[:, 0])
        u = (theta + np.pi) / (2 * np.pi)
        
        z = v[:, 2]
        
        def normalize(arr):
            rng = arr.max() - arr.min()
            if rng < 1e-9: return np.zeros_like(arr)
            return (arr - arr.min()) / rng

        v_coord = normalize(z)
            
        return np.column_stack([u, v_coord])

class MetricsCalculator:
    """
    Calculates quantitative distortion metrics.
    """
    @staticmethod
    def compute_metrics(mesh, uv):
        if uv is None or len(uv) != len(mesh.vertices):
            return {'area_std': 0.0, 'angle_mae': 0.0}

        try:
            uv_3d = np.column_stack([uv, np.zeros(len(uv))])
            mesh_2d = trimesh.Trimesh(vertices=uv_3d, faces=mesh.faces, process=False)
            
            # 1. Area Distortion
            areas_3d = mesh.area_faces
            areas_2d = mesh_2d.area_faces
            
            valid_mask = areas_3d > 1e-12
            a3 = areas_3d[valid_mask]
            a2 = areas_2d[valid_mask]
            
            if len(a3) == 0:
                return {'area_std': 0.0, 'angle_mae': 0.0}

            ratios = a2 / a3
            global_scale = np.mean(ratios)
            if global_scale == 0: global_scale = 1.0
            
            normalized_ratios = ratios / global_scale
            area_std = np.std(normalized_ratios)
            
            # 2. Angle Distortion
            angles_3d = np.degrees(mesh.face_angles)
            angles_2d = np.degrees(mesh_2d.face_angles)
            
            angle_diff = np.abs(angles_3d - angles_2d)
            angle_mae = np.mean(angle_diff[valid_mask])
            
            return {
                'area_std': round(float(area_std), 4),
                'angle_mae': round(float(angle_mae), 2)
            }
        except Exception:
            return {'area_std': 0.0, 'angle_mae': 0.0}

class BenchmarkRunner:
    def __init__(self, pdb_path, chain_a, chain_b, cutoff=5.0, output_root="benchmark_results"):
        self.pdb_path = pdb_path
        self.chain_a = chain_a
        self.chain_b = chain_b
        self.cutoff = cutoff
        self.output_root = output_root
        
        # 1. Extract PDB ID and Determine Category
        pdb_filename = os.path.basename(pdb_path)
        pdb_id = os.path.splitext(pdb_filename)[0].split('_')[0] # Handle suffixes if any
        # Default category if not provided externally (e.g. single file run)
        self.category = "Uncategorized"
        
        # 2. Create Output Directory (Base)
        # Category subfolder will be created in run() or set externally
        if not os.path.exists(self.output_root):
            os.makedirs(self.output_root)
            
        self.param = Parameterizer()
        
    def run(self):
        logger.info(f"--- Starting Benchmark for {os.path.basename(self.pdb_path)} ---")
        
        # Setup specific output path based on current category
        self.category_dir = os.path.join(self.output_root, self.category)
        if not os.path.exists(self.category_dir):
            os.makedirs(self.category_dir)
            
        pdb_id = os.path.splitext(os.path.basename(self.pdb_path))[0]
        self.output_csv = os.path.join(self.category_dir, f"{pdb_id}_benchmark.csv")

        # 1. Load & Process
        loader = PDBLoader(self.pdb_path)
        try:
            coords_A, _ = loader.get_chain_data(self.chain_a)
            coords_B, _ = loader.get_chain_data(self.chain_b)
        except Exception as e:
            logger.error(f"Failed to load PDB: {e}")
            return

        surf_gen = SurfaceGenerator(coords_A)
        mesh_A = surf_gen.generate_mesh(grid_resolution=1.0, sigma=1.5)
        if mesh_A is None: return

        topo = TopologyManager(mesh_A, coords_B)
        patches = topo.get_interface_patches(distance_cutoff=self.cutoff, min_patch_vertices=50)
        
        if not patches:
            logger.error("No valid interface patches found.")
            return

        # Prepare CSV Header
        fieldnames = [
            'PDB', 'Category', 'Patch_ID', 'Vertices', 
            'TopoPPI_Area_Std', 'TopoPPI_Angle_MAE', 'TopoPPI_Time',
            'Planar_Area_Std', 'Planar_Angle_MAE',
            'Spherical_Area_Std', 'Spherical_Angle_MAE',
            'Cylindrical_Area_Std', 'Cylindrical_Angle_MAE'
        ]
        
        patch_data_list = []
        results_rows = []

        # 2. Process Each Patch
        for i, patch in enumerate(patches):
            patch_name = f"{os.path.basename(self.pdb_path)}_P{i+1}"
            logger.info(f"Benchmarking Patch {i+1} ({len(patch.vertices)} verts)...")
            
            # --- Method A: TopoPPI (LSCM) ---
            t0 = time.time()
            uv_lscm = self.param.flatten_patch(patch)
            t_lscm = time.time() - t0
            res_lscm = MetricsCalculator.compute_metrics(patch, uv_lscm)
            
            # --- Method B: Planar ---
            uv_planar = BenchmarkAlgorithms.planar_projection(patch)
            res_planar = MetricsCalculator.compute_metrics(patch, uv_planar)

            # --- Method C: Spherical ---
            uv_sph = BenchmarkAlgorithms.spherical_projection(patch)
            res_sph = MetricsCalculator.compute_metrics(patch, uv_sph)
            
            # --- Method D: Cylindrical ---
            uv_cyl = BenchmarkAlgorithms.cylindrical_projection(patch)
            res_cyl = MetricsCalculator.compute_metrics(patch, uv_cyl)
            
            row = {
                'PDB': os.path.basename(self.pdb_path),
                'Category': self.category,
                'Patch_ID': str(i + 1),
                'Vertices': len(patch.vertices),
                
                'TopoPPI_Area_Std': res_lscm['area_std'],
                'TopoPPI_Angle_MAE': res_lscm['angle_mae'],
                'TopoPPI_Time': round(t_lscm, 4),
                
                'Planar_Area_Std': res_planar['area_std'],
                'Planar_Angle_MAE': res_planar['angle_mae'],

                'Spherical_Area_Std': res_sph['area_std'],
                'Spherical_Angle_MAE': res_sph['angle_mae'],
                
                'Cylindrical_Area_Std': res_cyl['area_std'],
                'Cylindrical_Angle_MAE': res_cyl['angle_mae']
            }
            results_rows.append(row)
            patch_data_list.append(row)
            
            print(f"  [Result P{i+1}] TopoPPI Angle Err: {row['TopoPPI_Angle_MAE']}° | Sphere: {row['Spherical_Angle_MAE']}°")

        # 3. Compute Aggregates
        if patch_data_list:
            total_verts = sum(p['Vertices'] for p in patch_data_list)
            
            # A. Weighted Average
            avg_row = {
                'PDB': os.path.basename(self.pdb_path), 
                'Category': self.category,
                'Patch_ID': 'weighted_avg', 
                'Vertices': total_verts
            }
            
            metrics = [
                'TopoPPI_Area_Std', 'TopoPPI_Angle_MAE', 'TopoPPI_Time',
                'Planar_Area_Std', 'Planar_Angle_MAE',
                'Spherical_Area_Std', 'Spherical_Angle_MAE',
                'Cylindrical_Area_Std', 'Cylindrical_Angle_MAE'
            ]
            
            for m in metrics:
                weighted_sum = sum(p[m] * p['Vertices'] for p in patch_data_list)
                avg_row[m] = round(weighted_sum / total_verts, 4)
            
            results_rows.append(avg_row)
            
            # B. Largest Patch
            largest_patch = max(patch_data_list, key=lambda x: x['Vertices'])
            largest_row = largest_patch.copy()
            largest_row['Patch_ID'] = 'largest'
            results_rows.append(largest_row)

            print(f"\n  [Summary] Weighted Avg TopoPPI Angle: {avg_row['TopoPPI_Angle_MAE']}°")

        # 4. Save to CSV
        try:
            with open(self.output_csv, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(results_rows)
            logger.info(f"Results saved to {self.output_csv}")
        except Exception as e:
            logger.error(f"Failed to save CSV: {e}")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="TopoPPI Benchmark Tool")
    parser.add_argument("path", help="Single PDB file path or Directory path")
    parser.add_argument("chain_a", nargs='?', help="Chain A ID (Required for Single File mode)")
    parser.add_argument("chain_b", nargs='?', help="Chain B ID (Required for Single File mode)")
    
    args = parser.parse_args()

    input_path = args.path
    
    if os.path.isdir(input_path):
        # --- FOLDER MODE (BATCH) ---
        print(f"Running Batch Benchmark in: {input_path}")
        
        # 1. Look for CSV in folder or fallback to root
        csv_path = os.path.join(input_path, "benchmark_targets.csv")
        if not os.path.exists(csv_path):
            csv_path = "benchmark_targets.csv" # Fallback to current directory
            
        if not os.path.exists(csv_path):
            print("Error: 'benchmark_targets.csv' not found in folder or current directory.")
            sys.exit(1)
            
        # 2. Read tasks
        tasks = []
        with open(csv_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                tasks.append(row)
        
        print(f"Loaded {len(tasks)} tasks from {csv_path}")
        
        # 3. Execute
        for i, task in enumerate(tasks):
            pdb_id = task['PDB']
            pdb_file = os.path.join(input_path, f"{pdb_id}.pdb")
            
            if not os.path.exists(pdb_file):
                print(f"[{i+1}/{len(tasks)}] Skipping {pdb_id}: File not found.")
                continue
                
            print(f"\n[{i+1}/{len(tasks)}] Processing {pdb_id}...")
            try:
                runner = BenchmarkRunner(pdb_file, task['Chain_A'], task['Chain_B'])
                # Override output root to be inside the input folder
                runner.output_root = os.path.join(input_path, "benchmark_results")
                # Set category explicitly from CSV
                runner.category = task['Category'] 
                
                runner.run()
            except Exception as e:
                print(f"Error processing {pdb_id}: {e}")

    elif os.path.isfile(input_path):
        # --- SINGLE FILE MODE ---
        if not args.chain_a or not args.chain_b:
            print("Error: Single file mode requires Chain A and Chain B arguments.")
            print("Usage: python benchmark.py <file.pdb> <ChainA> <ChainB>")
            sys.exit(1)
            
        runner = BenchmarkRunner(input_path, args.chain_a, args.chain_b)
        runner.category = "Single_Runs" # Default category for manual runs
        runner.run()
        
    else:
        print(f"Error: Path {input_path} not found.")
