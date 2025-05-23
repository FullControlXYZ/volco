"""
Mesh generation module for FEA.

This module provides functionality to convert voxel matrices to hexahedral elements
for finite element analysis.
"""

import numpy as np
from typing import Tuple, List
from scipy.ndimage import label


def generate_mesh(voxel_matrix: np.ndarray, voxel_size: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert a voxel matrix to a hexahedral FE mesh.
    
    This function generates nodes and element connectivity matrices from a voxel matrix.
    Each voxel is converted to a hexahedral element with 8 nodes.
    
    Args:
        voxel_matrix: 3D numpy array representing the voxel model (1 = material, 0 = void)
        voxel_size: Size of each voxel in mm
        
    Returns:
        Tuple containing:
            - nodes: Array of node coordinates (N x 3)
            - elements: Array of element connectivity (E x 8)
    """
    # Find indices of non-zero voxels (material voxels)
    material_indices = np.argwhere(voxel_matrix > 0)
    
    # Initialize lists to store nodes and elements
    nodes_list = []
    elements_list = []
    
    # Dictionary to map node coordinates to node indices
    node_dict = {}
    node_counter = 0
    
    # For each material voxel, create a hexahedral element
    for idx in material_indices:
        i, j, k = idx
        
        # Define the 8 corners of the hexahedron (local node ordering)
        corners = [
            (i, j, k),           # Node 0: (x, y, z)
            (i+1, j, k),         # Node 1: (x+1, y, z)
            (i+1, j+1, k),       # Node 2: (x+1, y+1, z)
            (i, j+1, k),         # Node 3: (x, y+1, z)
            (i, j, k+1),         # Node 4: (x, y, z+1)
            (i+1, j, k+1),       # Node 5: (x+1, y, z+1)
            (i+1, j+1, k+1),     # Node 6: (x+1, y+1, z+1)
            (i, j+1, k+1)        # Node 7: (x, y+1, z+1)
        ]
        
        # Element connectivity (node indices)
        element = []
        
        # For each corner, get or create a node
        for corner in corners:
            if corner not in node_dict:
                # Convert voxel indices to physical coordinates
                x, y, z = corner
                node_coords = np.array([x, y, z]) * voxel_size
                nodes_list.append(node_coords)
                
                # Store the node index
                node_dict[corner] = node_counter
                node_counter += 1
            
            # Add the node index to the element connectivity
            element.append(node_dict[corner])
        
        # Add the element to the list
        elements_list.append(element)
    
    # Convert lists to numpy arrays
    nodes = np.array(nodes_list)
    elements = np.array(elements_list)
    
    return nodes, elements


def check_continuity(voxel_matrix: np.ndarray) -> bool:
    """
    Check if the voxel matrix represents a single continuous body.
    
    This function uses connected component labeling to determine if all material
    voxels form a single connected component.
    
    Args:
        voxel_matrix: 3D numpy array representing the voxel model (1 = material, 0 = void)
        
    Returns:
        True if the voxel matrix represents a single continuous body, False otherwise
    """
    # Create a binary matrix (1 = material, 0 = void)
    binary_matrix = (voxel_matrix > 0).astype(int)
    
    # Use connected component labeling to identify separate bodies
    # structure defines connectivity (26-connectivity for 3D)
    structure = np.ones((3, 3, 3))
    labeled_matrix, num_features = label(binary_matrix, structure=structure)
    
    # Check if there is only one feature (excluding background)
    return num_features == 1


def extract_surface_elements(elements: np.ndarray, voxel_matrix: np.ndarray) -> np.ndarray:
    """
    Extract elements that are on the surface of the model.
    
    This is useful for applying boundary conditions to specific surfaces.
    
    Args:
        elements: Array of element connectivity
        voxel_matrix: 3D numpy array representing the voxel model
        
    Returns:
        Array of indices of surface elements
    """
    # Get dimensions of the voxel matrix
    nx, ny, nz = voxel_matrix.shape
    
    # Find material voxels
    material_indices = np.argwhere(voxel_matrix > 0)
    
    # Initialize list to store surface element indices
    surface_elements = []
    
    # For each material voxel, check if it's on the surface
    for idx, (i, j, k) in enumerate(material_indices):
        # Check the six neighboring voxels
        neighbors = [
            (i-1, j, k), (i+1, j, k),
            (i, j-1, k), (i, j+1, k),
            (i, j, k-1), (i, j, k+1)
        ]
        
        # If any neighbor is outside the matrix or is void, this is a surface element
        for ni, nj, nk in neighbors:
            if (ni < 0 or ni >= nx or 
                nj < 0 or nj >= ny or 
                nk < 0 or nk >= nz or 
                voxel_matrix[ni, nj, nk] == 0):
                surface_elements.append(idx)
                break
    
    return np.array(surface_elements)


def get_top_surface_nodes(nodes: np.ndarray, elements: np.ndarray, voxel_matrix: np.ndarray) -> np.ndarray:
    """
    Get nodes on the top surface of the model.
    
    This is useful for applying displacement boundary conditions.
    
    Args:
        nodes: Array of node coordinates
        elements: Array of element connectivity
        voxel_matrix: 3D numpy array representing the voxel model
        
    Returns:
        Array of indices of nodes on the top surface
    """
    # Get the maximum z-coordinate (top surface)
    max_z = np.max(nodes[:, 2])
    
    # Find nodes with z-coordinate equal to max_z (with small tolerance)
    tolerance = 1e-6
    top_nodes = np.where(np.abs(nodes[:, 2] - max_z) < tolerance)[0]
    
    return top_nodes


def get_bottom_surface_nodes(nodes: np.ndarray, elements: np.ndarray, voxel_matrix: np.ndarray) -> np.ndarray:
    """
    Get nodes on the bottom surface of the model.
    
    This is useful for applying fixed boundary conditions.
    
    Args:
        nodes: Array of node coordinates
        elements: Array of element connectivity
        voxel_matrix: 3D numpy array representing the voxel model
        
    Returns:
        Array of indices of nodes on the bottom surface
    """
    # Get the minimum z-coordinate (bottom surface)
    min_z = np.min(nodes[:, 2])
    
    # Find nodes with z-coordinate equal to min_z (with small tolerance)
    tolerance = 1e-6
    bottom_nodes = np.where(np.abs(nodes[:, 2] - min_z) < tolerance)[0]
    
    return bottom_nodes