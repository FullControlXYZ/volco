"""
Boundary condition definitions for FEA.

This module provides functionality to define and apply boundary conditions
for finite element analysis.
"""

import numpy as np
from typing import Tuple, Dict, Any, Callable, Optional, List, Union
from .mesh import get_top_surface_nodes, get_bottom_surface_nodes


def apply_default_boundary_conditions(
    nodes: np.ndarray,
    voxel_matrix: np.ndarray,
    **kwargs
) -> Tuple[Dict[int, List[float]], np.ndarray]:
    """
    Apply default boundary conditions to the FE model.
    
    Default boundary conditions are:
    - Bottom surface: Fixed (all DOFs constrained)
    - Top surface: Prescribed displacement (1% downward)
    
    Args:
        nodes: Array of node coordinates (N x 3)
        voxel_matrix: 3D numpy array representing the voxel model
        **kwargs: Optional parameters to customize the default boundary conditions:
            - 'displacement_percentage': Percentage of model height for top displacement (default: 1.0)
            - 'custom_top_bc': Custom function for top surface boundary conditions
            - 'custom_bottom_bc': Custom function for bottom surface boundary conditions
            
    Returns:
        Tuple containing:
            - dofs: Dictionary mapping node indices to prescribed displacements
              {node_idx: [dx, dy, dz, rx, ry, rz]} where None means unconstrained
            - forces: Array of nodal forces (N x 6)
    """
    # Get model height (max z - min z)
    model_height = np.max(nodes[:, 2]) - np.min(nodes[:, 2])
    
    # Get displacement percentage (default: 1%)
    displacement_percentage = kwargs.get('displacement_percentage', 1.0)
    displacement_magnitude = model_height * displacement_percentage / 100.0
    
    # Get nodes on top and bottom surfaces
    top_nodes = get_top_surface_nodes(nodes, None, voxel_matrix)
    bottom_nodes = get_bottom_surface_nodes(nodes, None, voxel_matrix)
    
    # Initialize DOFs dictionary and forces array
    dofs = {}
    forces = np.zeros((nodes.shape[0], 6))  # [Fx, Fy, Fz, Mx, My, Mz]
    
    # Apply custom or default boundary conditions
    if 'custom_bottom_bc' in kwargs:
        custom_bottom_bc = kwargs['custom_bottom_bc']
        bottom_dofs = custom_bottom_bc(nodes, bottom_nodes)
        dofs.update(bottom_dofs)
    else:
        # Default: Fix bottom surface (constrain all DOFs)
        for node_idx in bottom_nodes:
            dofs[node_idx] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # [dx, dy, dz, rx, ry, rz]
    
    if 'custom_top_bc' in kwargs:
        custom_top_bc = kwargs['custom_top_bc']
        top_dofs = custom_top_bc(nodes, top_nodes, displacement_magnitude)
        dofs.update(top_dofs)
    else:
        # Default: Prescribe displacement on top surface (negative z direction)
        for node_idx in top_nodes:
            # Allow x and y displacement, constrain z to prescribed value
            # Allow rotations (None means unconstrained)
            dofs[node_idx] = [None, None, -displacement_magnitude, None, None, None]
    
    return dofs, forces


def apply_custom_boundary_conditions(
    nodes: np.ndarray,
    bc_function: Callable,
    **kwargs
) -> Tuple[Dict[int, List[float]], np.ndarray]:
    """
    Apply custom boundary conditions using a user-defined function.
    
    Args:
        nodes: Array of node coordinates (N x 3)
        bc_function: Function that takes nodes and kwargs and returns dofs and forces
        **kwargs: Additional parameters for the bc_function
        
    Returns:
        Tuple containing:
            - dofs: Dictionary mapping node indices to prescribed displacements
            - forces: Array of nodal forces
    """
    return bc_function(nodes, **kwargs)


def create_fixed_support(node_indices: List[int]) -> Dict[int, List[float]]:
    """
    Create fixed support boundary conditions for specified nodes.
    
    Args:
        node_indices: List of node indices to fix
        
    Returns:
        Dictionary mapping node indices to prescribed displacements (all zero)
    """
    dofs = {}
    for node_idx in node_indices:
        dofs[node_idx] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # [dx, dy, dz, rx, ry, rz]
    return dofs


def create_displacement_bc(
    node_indices: List[int],
    displacement: List[Optional[float]]
) -> Dict[int, List[Optional[float]]]:
    """
    Create displacement boundary conditions for specified nodes.
    
    Args:
        node_indices: List of node indices to apply displacement to
        displacement: List of 6 values [dx, dy, dz, rx, ry, rz] where None means unconstrained
        
    Returns:
        Dictionary mapping node indices to prescribed displacements
    """
    dofs = {}
    for node_idx in node_indices:
        dofs[node_idx] = displacement
    return dofs


def create_force_bc(
    node_indices: List[int],
    force: List[float],
    forces: np.ndarray
) -> np.ndarray:
    """
    Create force boundary conditions for specified nodes.
    
    Args:
        node_indices: List of node indices to apply force to
        force: List of 6 values [Fx, Fy, Fz, Mx, My, Mz]
        forces: Existing forces array to update
        
    Returns:
        Updated forces array
    """
    for node_idx in node_indices:
        forces[node_idx] = force
    return forces


def select_nodes_by_position(
    nodes: np.ndarray,
    x_range: Optional[Tuple[float, float]] = None,
    y_range: Optional[Tuple[float, float]] = None,
    z_range: Optional[Tuple[float, float]] = None
) -> List[int]:
    """
    Select nodes within specified coordinate ranges.
    
    Args:
        nodes: Array of node coordinates
        x_range: Tuple of (min_x, max_x) or None to ignore
        y_range: Tuple of (min_y, max_y) or None to ignore
        z_range: Tuple of (min_z, max_z) or None to ignore
        
    Returns:
        List of node indices that fall within the specified ranges
    """
    mask = np.ones(nodes.shape[0], dtype=bool)
    
    if x_range is not None:
        mask &= (nodes[:, 0] >= x_range[0]) & (nodes[:, 0] <= x_range[1])
    
    if y_range is not None:
        mask &= (nodes[:, 1] >= y_range[0]) & (nodes[:, 1] <= y_range[1])
    
    if z_range is not None:
        mask &= (nodes[:, 2] >= z_range[0]) & (nodes[:, 2] <= z_range[1])
    
    return np.where(mask)[0].tolist()