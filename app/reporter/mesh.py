import logging
import trimesh
import numpy as np
from skimage import measure

logger = logging.getLogger(__name__)

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
    
    # Create a mesh from the voxel data using box representation
    try:
        # Create a VoxelGrid object
        mesh = trimesh.voxel.base.VoxelGrid(voxel_space)
        
        # Set origin and scale
        origin = (0.0, 0.0, 0.0)
        scale = (voxel_size, voxel_size, voxel_size)
        
        mesh.origin[0], mesh.origin[1], mesh.origin[2] = origin
        mesh.scale.setflags(write=1)
        mesh.scale[0], mesh.scale[1], mesh.scale[2] = scale
        
        # Convert to boxes representation
        voxels = mesh.as_boxes(colors=None)
        return voxels
    except Exception as e:
        logger.error(f"[Mesh]: Failed to generate mesh with box representation: {e}")
        return None

    # NOTE: Previous implementation using marching cubes algorithm (commented for future reference)
    """
    # Create a mesh from the voxel data using trimesh
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
    """

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
        
        # Check if mesh is a Scene object (from as_boxes()) or a regular Trimesh
        if isinstance(mesh, trimesh.Scene):
            # For Scene objects, we need to handle differently
            # First convert to a Trimesh if possible
            export_options = {}
            if ascii_format:
                # Use the STLExporter's ASCII option
                export_options = {'file_type': 'stl_ascii'}
            else:
                export_options = {'file_type': 'stl'}
                
            _ = trimesh.exchange.export.export_mesh(mesh, file_path, **export_options)
        else:
            # For regular Trimesh objects
            if ascii_format:
                # ASCII STL export
                with open(file_path, 'w') as f:
                    f.write(trimesh.exchange.stl.export_stl_ascii(mesh))
            else:
                # Binary STL export (default)
                mesh.export(file_path, file_type='stl')
            
        logger.info(f"[Mesh]: STL exported in {format_type} format!")
        return file_path
    except Exception as e:
        logger.error(f"[Mesh]: Failed to export STL: {e}")
        return None