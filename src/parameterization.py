import numpy as np
import igl
import trimesh
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("LSCM")

class Parameterizer:
    """
    Handles the flattening of 3D meshes into 2D UV coordinates using LSCM
    (Least Squares Conformal Maps), with a Harmonic fallback.
    """
    
    @staticmethod
    def flatten_patch(mesh: trimesh.Trimesh) -> np.ndarray:
        """
        Flatten a 3D mesh patch to 2D.
        
        Args:
            mesh: A single connected component (trimesh object).
            
        Returns:
            np.ndarray: UV coordinates of shape (N, 2), or None if failed.
        """
        # --- Step 0: Robust Mesh Sanitation ---
        # Critical Fix: Smoothing in TopologyManager can introduce degenerate (zero-area) faces.
        # These cause the LSCM linear solver to fail (singular matrix). 
        # We must clean the mesh immediately before parameterization to ensure stability.
        try:
            # 1. Merge Vertices (Topological Stitching)
            # Ensures that geometrically identical vertices are topologically merged.
            # This is crucial for LSCM which requires a single connected component.
            mesh.merge_vertices()

            # 2. Remove degenerate faces (area approx 0)
            valid_faces = mesh.area_faces > 1e-12
            if not np.all(valid_faces):
                mesh.update_faces(valid_faces)
            
            # 3. Safe remove duplicate faces (Fix for AttributeError on older trimesh versions)
            if hasattr(mesh, 'remove_duplicate_faces'):
                mesh.remove_duplicate_faces()
            else:
                # Fallback: Process usually handles duplicates if specific method is missing,
                # but we avoid full process() to keep vertex order if possible.
                # If method is missing, we trust merge_vertices did enough.
                pass

            # 4. Remove unused vertices (Critical so indices match 0..N-1 for IGL)
            mesh.remove_unreferenced_vertices()
            
            # 5. Component Check: Ensure we really have one connected component
            # Sometimes 'remove_degenerate_faces' can disconnect the mesh.
            # LSCM cannot handle multiple disconnected components in one call.
            components = mesh.split(only_watertight=False)
            if len(components) > 1:
                # Keep only the largest component by vertex count
                mesh = max(components, key=lambda m: len(m.vertices))
            
            # Check if mesh is still valid
            if len(mesh.vertices) < 3 or len(mesh.faces) == 0:
                logger.warning("Mesh became empty or degenerate after cleanup.")
                return None
                
        except Exception as e:
            logger.warning(f"Mesh sanitation in Parameterizer failed: {e}")
            return None

        # --- Step 1: Prepare IGL Data ---
        # IGL is strict about types. Ensure correct C++ compatible types.
        # Use np.ascontiguousarray to strip Trimesh wrappers and ensure C-order.
        v = np.ascontiguousarray(mesh.vertices, dtype=np.float64)
        f = np.ascontiguousarray(mesh.faces, dtype=np.int64)

        # 2. Find Boundary Loop (LSCM needs a boundary)
        # igl.boundary_loop returns the ordered vertex indices of the boundary
        try:
            bnd = igl.boundary_loop(f)
        except Exception as e:
            logger.error(f"Failed to detect boundary: {e}")
            return None
        
        # Handle case where multiple boundaries might be returned (list of lists)
        if len(bnd) > 0 and isinstance(bnd[0], (list, np.ndarray)):
            # Heuristic: Pick the longest boundary loop (the "outer" one)
            # Short loops are likely internal holes; LSCM can handle them as free boundaries.
            bnd = sorted(bnd, key=lambda x: len(x), reverse=True)[0]
        
        bnd = np.array(bnd, dtype=np.int64)

        if len(bnd) < 3:
            logger.error("Mesh has no valid boundary (closed or degenerate).")
            return None

        # 3. Fix Boundary Conditions for LSCM
        # Strategy: Pin the two most distant points on the boundary.
        b1_idx = bnd[0]
        
        # Find point on boundary furthest from A (Euclidean distance)
        boundary_coords = v[bnd]
        dists = np.linalg.norm(boundary_coords - v[b1_idx], axis=1)
        b2_idx = bnd[np.argmax(dists)]
        
        # Constraints inputs: b (indices), bc (target coords)
        # FIX: Ensure 'b' is also int64 to match 'f'
        b = np.array([b1_idx, b2_idx], dtype=np.int64)
        bc = np.array([[0.0, 0.0], [1.0, 0.0]], dtype=np.float64)

        # 4. Run LSCM
        uv_normalized = None
        try:
            # IGL LSCM signature: lscm(V, F, b, bc) -> (success, UV)
            ret = igl.lscm(v, f, b, bc)
            
            # Robust unpacking for different libigl versions
            if isinstance(ret, tuple) and len(ret) == 2:
                success, uv = ret
            else:
                success = True
                uv = ret

            # Handle numpy/bool ambiguity
            if isinstance(success, np.ndarray):
                is_success = success.all()
            else:
                is_success = bool(success)
            
            if is_success:
                uv_normalized = Parameterizer._normalize_uv(uv)
            else:
                logger.warning("IGL LSCM solver returned failure status (Matrix likely singular).")
                
        except Exception as e:
            logger.warning(f"LSCM Exception: {e}")

        # 5. Fallback: Harmonic Parameterization
        # Only triggered if LSCM absolutely fails
        if uv_normalized is None:
            logger.info("Attempting Harmonic Parameterization fallback...")
            uv_normalized = Parameterizer._flatten_harmonic(v, f, bnd)

        return uv_normalized

    @staticmethod
    def _flatten_harmonic(v, f, bnd):
        """
        Fallback method: Map boundary to circle and minimize Dirichlet energy.
        """
        try:
            # 1. Map boundary vertices to a circle
            # Ensure bnd is int64
            bnd = bnd.astype(np.int64)
            bnd_uv = igl.map_vertices_to_circle(v, bnd)
            
            # 2. Harmonic parameterization (power=1)
            # harmonic(V, F, b, bc, k)
            uv = igl.harmonic(v, f, bnd, bnd_uv, 1)
            
            return Parameterizer._normalize_uv(uv)
        except Exception as e:
            logger.error(f"Harmonic Parameterization failed: {e}")
            return None

    @staticmethod
    def _normalize_uv(uv):
        """Helper to normalize UV to [0,1] range."""
        if uv is None or len(uv) == 0:
            return None
        uv_min = uv.min(axis=0)
        uv_max = uv.max(axis=0)
        scale = uv_max - uv_min
        scale[scale < 1e-6] = 1.0 
        return (uv - uv_min) / scale

# --- Self-Contained Unit Test ---
if __name__ == "__main__":
    pass