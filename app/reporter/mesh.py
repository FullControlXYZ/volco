import logging
import trimesh
import numpy as np
from skimage import measure

logger = logging.getLogger(__name__)

# Check if we're using NumPy 2.x
USING_NUMPY2 = int(np.__version__.split('.')[0]) >= 2
if USING_NUMPY2:
    logger.info(f"[Mesh]: Detected NumPy {np.__version__}, using NumPy 2.x compatible methods")

def create_box_representation(voxel_space, voxel_size):
    """
    Create a box representation mesh from a voxel space.
    This is a custom implementation that doesn't rely on any trimesh functions
    that might use the `ptp` method, making it compatible with NumPy 2.x.
    
    Parameters:
    -----------
    voxel_space : numpy.ndarray
        The 3D voxel space array
    voxel_size : float
        The size of each voxel
        
    Returns:
    --------
    trimesh.Scene
        A scene containing a box for each filled voxel
    """
    # Find the indices of filled voxels
    filled_voxels = np.where(voxel_space > 0)
    
    # If no filled voxels, return an empty scene
    if len(filled_voxels[0]) == 0:
        return trimesh.Scene()
    
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
    
    # Define the 12 triangles (6 faces, 2 triangles per face)
    # Ensure normals point outward by ordering vertices counter-clockwise
    unit_cube_faces = np.array([
        [0, 2, 1], [0, 3, 2],  # bottom face
        [4, 5, 6], [4, 6, 7],  # top face
        [0, 1, 5], [0, 5, 4],  # back face
        [2, 3, 7], [2, 7, 6],  # front face
        [0, 4, 7], [0, 7, 3],  # left face
        [1, 2, 6], [1, 6, 5]   # right face
    ])
    
    # Calculate total number of vertices and faces
    num_voxels = len(filled_voxels[0])
    total_vertices = num_voxels * 8
    total_faces = num_voxels * 12
    
    # Pre-allocate arrays for all vertices and faces
    all_vertices = np.zeros((total_vertices, 3))
    all_faces = np.zeros((total_faces, 3), dtype=np.int64)
    
    # For each filled voxel, add a cube
    for idx, (i, j, k) in enumerate(zip(*filled_voxels)):
        # Calculate the center position of the voxel
        # This is set to match the trimesh output, it's not verified
        # against the original voxel space
        center = np.array([
            (i) * voxel_size,
            (j) * voxel_size,
            (k) * voxel_size
        ])
        
        # Scale and translate the unit cube vertices
        voxel_vertices = unit_cube_vertices * voxel_size + center
        
        # Add vertices to the array
        vertex_start = idx * 8
        vertex_end = vertex_start + 8
        all_vertices[vertex_start:vertex_end] = voxel_vertices
        
        # Add faces to the array (with offset for vertex indices)
        face_start = idx * 12
        face_end = face_start + 12
        all_faces[face_start:face_end] = unit_cube_faces + (idx * 8)
    
    # Create a single mesh from all vertices and faces
    mesh = trimesh.Trimesh(vertices=all_vertices, faces=all_faces)
    
    # Return the mesh directly (not wrapped in a Scene)
    return mesh

def generate_mesh_from_voxels(voxel_space, voxel_size):
    """
    Generate a 3D mesh from the voxel data using box representation.
    
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
    
    # For NumPy 2.x, use our custom box representation function
    if USING_NUMPY2:
        try:
            logger.info("[Mesh]: Using NumPy 2.x compatible box representation")
            # Use our custom box representation function
            voxels = create_box_representation(voxel_space, voxel_size)
            return voxels
        except Exception as e:
            logger.error(f"[Mesh]: Failed to generate mesh with custom box representation: {e}")
            return None
    else:
        # For NumPy 1.x, use the original approach
        try:
            # Create a VoxelGrid object
            mesh = trimesh.voxel.base.VoxelGrid(voxel_space)
            
            # Set origin and scale
            origin = (0.0, 0.0, 0.0)
            scale = (voxel_size, voxel_size, voxel_size)
            
            # Check if the mesh has origin and scale attributes before setting them
            if hasattr(mesh, 'origin'):
                mesh.origin[0], mesh.origin[1], mesh.origin[2] = origin
            if hasattr(mesh, 'scale'):
                mesh.scale.setflags(write=1)
                mesh.scale[0], mesh.scale[1], mesh.scale[2] = scale
            
            # Convert to boxes representation
            voxels = mesh.as_boxes(colors=None)
            return voxels
        except Exception as e:
            logger.warning(f"[Mesh]: Failed to generate mesh with box representation: {e}")
            # Fall back to using marching cubes
            try:
                mesh = trimesh.voxel.ops.matrix_to_marching_cubes(voxel_space, pitch=voxel_size)
                return mesh
            except Exception as e:
                logger.warning(f"[Mesh]: Failed to generate mesh with trimesh: {e}")
                # Fall back to using skimage directly
                try:
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