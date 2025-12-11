import numpy as np
import trimesh
import trimesh.smoothing
from skimage import measure
from scipy.ndimage import gaussian_filter
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("Surface")

class SurfaceGenerator:
    """
    Generates a smooth Solvent-Excluded Surface (SES) approximation 
    from atomic coordinates using Gaussian Density Fields.
    """
    
    def __init__(self, coords: np.ndarray):
        """
        Args:
            coords: Atom coordinates (N, 3) numpy array.
        """
        self.coords = coords

    def generate_mesh(self, grid_resolution: float = 1.0, sigma: float = 1.5, level: float = 0.1) -> trimesh.Trimesh:
        """
        Run the pipeline: Voxelization -> Density -> Isosurface -> Mesh.
        Includes robust retry logic for sparse surfaces.
        """
        num_atoms = len(self.coords)
        logger.info(f"Generating surface from {num_atoms} atoms...")
        
        if num_atoms == 0:
            logger.error("No coordinates provided.")
            return None

        # 1. Define Grid Bounds with Padding
        padding = 10.0
        min_bound = self.coords.min(axis=0) - padding
        max_bound = self.coords.max(axis=0) + padding
        
        # Calculate grid shape
        shape = np.ceil((max_bound - min_bound) / grid_resolution).astype(int)
        logger.info(f"Grid shape: {shape}, Resolution: {grid_resolution}A")

        # 2. Fast Voxelization
        grid, edges = np.histogramdd(
            self.coords, 
            bins=shape, 
            range=[(min_bound[i], max_bound[i]) for i in range(3)]
        )
        
        # 3. Compute Density Field
        density_field = gaussian_filter(grid.astype(float), sigma=sigma)
        
        max_density = density_field.max()
        if max_density == 0:
            logger.error("Density field is empty (all zeros). Check coordinates.")
            return None

        # --- Iterative Level Adjustment (Smart Retry) ---
        # Fix for 1AHW issue: If max_density is driven by outliers (e.g. clashing atoms),
        # the default level (0.1) might be too high for the rest of the surface.
        
        current_level = level
        # Safety: Ensure we start below the max
        if current_level >= max_density:
            current_level = max_density * 0.5

        final_mesh = None
        
        # Heuristic: Expect at least 0.5 vertices per atom for a decent coarse surface
        min_expected_verts = min(500, num_atoms * 0.5) 

        for attempt in range(4): # Try up to 4 times
            try:
                verts, faces, normals, values = measure.marching_cubes(
                    density_field, 
                    level=current_level,
                    step_size=1
                )
                
                n_verts = len(verts)
                # Check if surface is suspiciously small
                if n_verts < min_expected_verts and attempt < 3:
                    logger.warning(f"Surface too small ({n_verts} verts) at level {current_level:.4f}. Reducing threshold...")
                    current_level *= 0.5 # Halve the threshold to capture more volume
                    continue 
                
                # If acceptable size or last attempt
                real_verts = min_bound + verts * grid_resolution
                final_mesh = trimesh.Trimesh(vertices=real_verts, faces=faces, vertex_normals=normals)
                break
                
            except ValueError as e:
                logger.warning(f"Marching Cubes failed at level {current_level}: {e}")
                current_level *= 0.5
            except Exception as e:
                logger.error(f"Unexpected error: {e}")
                break

        if final_mesh is None or len(final_mesh.vertices) == 0:
            logger.error("Failed to generate a valid mesh after retries.")
            return None

        # Optional: Basic smoothing
        try:
            trimesh.smoothing.filter_laplacian(final_mesh, iterations=3)
        except Exception as e:
            logger.warning(f"Mesh smoothing skipped: {e}")
        
        logger.info(f"Surface generated: {len(final_mesh.vertices)} verts, {len(final_mesh.faces)} faces.")
        return final_mesh

if __name__ == "__main__":
    pass
