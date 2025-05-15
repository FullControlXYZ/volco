import logging
import trimesh
import numpy as np
from skimage import measure

logger = logging.getLogger(__name__)

# Log NumPy version for debugging purposes
logger.info(f"[Mesh]: Using NumPy {np.__version__}")

def create_box_representation(voxel_space, voxel_size):
    """
    Create a box representation mesh from a voxel space.
    This is a custom implementation that doesn't rely on any trimesh functions
    that might use the `ptp` method, making it compatible with NumPy 2.x.
    
    This optimized version only creates triangles for voxel faces that are either:
    1. At the boundary of the voxel matrix
    2. Adjacent to an empty voxel
    
    Parameters:
    -----------
    voxel_space : numpy.ndarray
        The 3D voxel space array
    voxel_size : float
        The size of each voxel
        
    Returns:
    --------
    trimesh.Trimesh
        A mesh containing only visible faces of filled voxels
    """
    # Find the indices of filled voxels
    filled_voxels = np.where(voxel_space > 0)
    
    # If no filled voxels, return an empty scene
    if len(filled_voxels[0]) == 0:
        return trimesh.Scene()
    
    # Get the dimensions of the voxel space
    max_i, max_j, max_k = voxel_space.shape
    
    # Define a unit cube vertices (8 corners)
    unit_cube_vertices = np.array([
        [-0.5, -0.5, -0.5],  # 0: bottom, back, left
        [0.5, -0.5, -0.5],   # 1: bottom, back, right
        [0.5, 0.5, -0.5],    # 2: bottom, front, right
        [-0.5, 0.5, -0.5],   # 3: bottom, front, left
        [-0.5, -0.5, 0.5],   # 4: top, back, left
        [0.5, -0.5, 0.5],    # 5: top, back, right
        [0.5, 0.5, 0.5],     # 6: top, front, right
        [-0.5, 0.5, 0.5]     # 7: top, front, left
    ])
    
    # Define the faces for each side of the cube (2 triangles per face)
    # Each face is associated with a direction (negative or positive x, y, z)
    face_definitions = [
        # Face indices, Direction to check (di, dj, dk)
        ([[0, 2, 1], [0, 3, 2]], (0, 0, -1)),  # bottom face (negative z)
        ([[4, 5, 6], [4, 6, 7]], (0, 0, 1)),   # top face (positive z)
        ([[0, 1, 5], [0, 5, 4]], (0, -1, 0)),  # back face (negative y)
        ([[2, 3, 7], [2, 7, 6]], (0, 1, 0)),   # front face (positive y)
        ([[0, 4, 7], [0, 7, 3]], (-1, 0, 0)),  # left face (negative x)
        ([[1, 2, 6], [1, 6, 5]], (1, 0, 0))    # right face (positive x)
    ]
    
    # First pass: determine which voxels need to have vertices added
    # and which faces need to be rendered
    voxels_to_render = {}  # Maps voxel index to list of faces to render
    
    for idx, (i, j, k) in enumerate(zip(*filled_voxels)):
        faces_to_render = []
        
        # Check each face to see if it should be rendered
        for face_idx, (face_triangles, (di, dj, dk)) in enumerate(face_definitions):
            # Check if this face is at the boundary or adjacent to an empty voxel
            ni, nj, nk = i + di, j + dj, k + dk
            
            # If the neighbor is outside the voxel space or is empty, render this face
            if (ni < 0 or ni >= max_i or
                nj < 0 or nj >= max_j or
                nk < 0 or nk >= max_k or
                voxel_space[ni, nj, nk] == 0):
                
                faces_to_render.append(face_idx)
        
        # If this voxel has faces to render, add it to the dictionary
        if faces_to_render:
            voxels_to_render[(i, j, k)] = faces_to_render
    
    # Second pass: create vertices and faces
    all_vertices = []
    all_faces = []
    voxel_to_vertex_idx = {}  # Maps voxel coordinates to starting vertex index
    
    for idx, (i, j, k) in enumerate(voxels_to_render.keys()):
        # Calculate the center position of the voxel
        center = np.array([
            (i) * voxel_size,
            (j) * voxel_size,
            (k) * voxel_size
        ])
        
        # Scale and translate the unit cube vertices for this voxel
        voxel_vertices = unit_cube_vertices * voxel_size + center
        
        # Add vertices for this voxel
        vertex_start = len(all_vertices)
        all_vertices.extend(voxel_vertices)
        voxel_to_vertex_idx[(i, j, k)] = vertex_start
        
        # Add faces for this voxel
        for face_idx in voxels_to_render[(i, j, k)]:
            face_triangles = face_definitions[face_idx][0]
            for triangle in face_triangles:
                all_faces.append([t + vertex_start for t in triangle])
    
    # Convert lists to numpy arrays
    if all_vertices and all_faces:
        vertices_array = np.array(all_vertices)
        faces_array = np.array(all_faces)
        
        # Create a mesh from the vertices and faces
        mesh = trimesh.Trimesh(vertices=vertices_array, faces=faces_array)
        return mesh
    else:
        # If no faces were added, return an empty scene
        return trimesh.Scene()

def generate_mesh_from_voxels(voxel_space, voxel_size):
    """
    Generate a 3D mesh from the voxel data using optimized box representation.
    This implementation only creates triangles for voxel faces that are either:
    1. At the boundary of the voxel matrix
    2. Adjacent to an empty voxel
    
    Parameters:
    -----------
    voxel_space : numpy.ndarray
        The 3D voxel space array
    voxel_size : float
        The size of each voxel
        
    Returns:
    --------
    trimesh.Trimesh
        The generated mesh
    """
    if voxel_space is None:
        logger.warning("[Mesh]: No voxel space provided.")
        return None
    
    try:
        logger.info("[Mesh]: Using optimized box representation")
        # Use our optimized box representation function
        voxels = create_box_representation(voxel_space, voxel_size)
        return voxels
    except Exception as e:
        logger.error(f"[Mesh]: Failed to generate mesh with optimized box representation: {e}")
        
        # Fall back to using marching cubes
        try:
            logger.warning("[Mesh]: Falling back to marching cubes method")
            mesh = trimesh.voxel.ops.matrix_to_marching_cubes(voxel_space, pitch=voxel_size)
            return mesh
        except Exception as e:
            logger.warning(f"[Mesh]: Failed to generate mesh with trimesh: {e}")
            
            # Fall back to using skimage directly
            try:
                logger.warning("[Mesh]: Falling back to skimage marching cubes")
                verts, faces, normals, values = measure.marching_cubes(voxel_space, level=0.5)
                
                # Scale vertices by voxel size
                verts = verts * voxel_size
                
                # Create mesh from vertices and faces
                mesh = trimesh.Trimesh(vertices=verts, faces=faces)
                return mesh
            except Exception as e:
                logger.error(f"[Mesh]: Failed to generate mesh with skimage: {e}")
                return None

def export_mesh_to_stl(mesh, file_path, ascii_format=True):
    """
    Export the mesh to an STL file.
    
    Parameters:
    -----------
    mesh : trimesh.Trimesh or trimesh.Scene
        The mesh to export (can be a Trimesh object or a Scene object from as_boxes())
    file_path : str
        The path to save the STL file
    ascii_format : bool
        Whether to export in ASCII format (True) or binary format (False)
            
    Returns:
    --------
    str
        The path to the exported STL file
    """
    if mesh is None:
        logger.warning("[Mesh]: No mesh provided for export.")
        return None
        
    try:
        format_type = "ASCII" if ascii_format else "binary"
        logger.info(f"[Mesh]: Exporting STL in {format_type} format to path {file_path}...")
        
        # Set export options based on format
        export_options = {'file_type': 'stl_ascii' if ascii_format else 'stl'}
        
        # Use the generic export function which works for both Scene and Trimesh objects
        trimesh.exchange.export.export_mesh(mesh, file_path, **export_options)
            
        logger.info(f"[Mesh]: STL exported in {format_type} format!")
        return file_path
    except Exception as e:
        logger.error(f"[Mesh]: Failed to export STL: {e}")
        return None