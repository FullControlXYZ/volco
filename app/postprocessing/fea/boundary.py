"""
Boundary condition definitions for FEA.

This module provides functionality to define and apply boundary conditions
for finite element analysis, with both Simple and Expert modes.
"""

import numpy as np
from enum import Enum, auto
from typing import Tuple, Dict, Any, Callable, Optional, List, Union, Set


class Surface(Enum):
    """Enumeration for standard surface identifiers"""
    PLUS_X = auto()  # +x surface
    MINUS_X = auto()  # -x surface
    PLUS_Y = auto()  # +y surface
    MINUS_Y = auto()  # -y surface
    PLUS_Z = auto()  # +z surface
    MINUS_Z = auto()  # -z surface


class BoundaryCondition:
    """Internal representation of boundary conditions"""
    def __init__(self, node_indices: List[int], dof_values: List[Optional[float]]):
        """
        Initialize a boundary condition.
        
        Args:
            node_indices: List of node indices to apply the boundary condition to
            dof_values: List of 6 values [ux, uy, uz, rx, ry, rz] where None means unconstrained
        """
        self.node_indices = node_indices
        self.dof_values = dof_values
        
    def apply(self, dofs: Dict[int, List[Optional[float]]]) -> Dict[int, List[Optional[float]]]:
        """
        Apply this boundary condition to the DOFs dictionary.
        
        Args:
            dofs: Dictionary mapping node indices to prescribed displacements
            
        Returns:
            Updated DOFs dictionary
        """
        for node_idx in self.node_indices:
            dofs[node_idx] = self.dof_values
        return dofs


def identify_surface_nodes(nodes: np.ndarray, surface: Surface) -> List[int]:
    """
    Identify nodes on a specified surface.
    
    Args:
        nodes: Array of node coordinates
        surface: Surface enumeration value
        
    Returns:
        List of node indices on the specified surface
    """
    # Get model dimensions
    x_min, y_min, z_min = np.min(nodes, axis=0)
    x_max, y_max, z_max = np.max(nodes, axis=0)
    
    # Small tolerance for floating-point comparisons
    tol = 1e-6
    
    if surface == Surface.PLUS_X:
        return [i for i, node in enumerate(nodes) if abs(node[0] - x_max) < tol]
    elif surface == Surface.MINUS_X:
        return [i for i, node in enumerate(nodes) if abs(node[0] - x_min) < tol]
    elif surface == Surface.PLUS_Y:
        return [i for i, node in enumerate(nodes) if abs(node[1] - y_max) < tol]
    elif surface == Surface.MINUS_Y:
        return [i for i, node in enumerate(nodes) if abs(node[1] - y_min) < tol]
    elif surface == Surface.PLUS_Z:
        return [i for i, node in enumerate(nodes) if abs(node[2] - z_max) < tol]
    elif surface == Surface.MINUS_Z:
        return [i for i, node in enumerate(nodes) if abs(node[2] - z_min) < tol]
    else:
        raise ValueError(f"Unknown surface: {surface}")


def select_nodes_by_predicate(nodes: np.ndarray, predicate: Callable) -> List[int]:
    """
    Select nodes that satisfy a given predicate function.
    
    Args:
        nodes: Array of node coordinates
        predicate: Function that takes a node coordinate and returns a boolean
        
    Returns:
        List of node indices that satisfy the predicate
    """
    return [i for i, node in enumerate(nodes) if predicate(node)]


def select_nodes_in_box(nodes: np.ndarray, min_coords: List[float], max_coords: List[float]) -> List[int]:
    """
    Select nodes within a bounding box.
    
    Args:
        nodes: Array of node coordinates
        min_coords: Minimum coordinates [x_min, y_min, z_min]
        max_coords: Maximum coordinates [x_max, y_max, z_max]
        
    Returns:
        List of node indices within the bounding box
    """
    return select_nodes_by_predicate(
        nodes,
        lambda node: all(min_coords[i] <= node[i] <= max_coords[i] for i in range(3))
    )


def select_nodes_on_plane(nodes: np.ndarray, point: List[float], normal: List[float], tolerance: float = 1e-6) -> List[int]:
    """
    Select nodes on a plane defined by a point and normal vector.
    
    Args:
        nodes: Array of node coordinates
        point: Point on the plane [x, y, z]
        normal: Normal vector to the plane [nx, ny, nz]
        tolerance: Distance tolerance for considering a node to be on the plane
        
    Returns:
        List of node indices on the plane
    """
    # Normalize the normal vector
    normal = np.array(normal)
    normal = normal / np.linalg.norm(normal)
    
    # Convert point to numpy array
    point = np.array(point)
    
    # Define predicate function for nodes on the plane
    def on_plane(node):
        # Calculate distance from node to plane
        vec = node - point
        distance = abs(np.dot(vec, normal))
        return distance < tolerance
    
    return select_nodes_by_predicate(nodes, on_plane)


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
        dofs[node_idx] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # [ux, uy, uz, rx, ry, rz]
    return dofs


def create_displacement_bc(
    node_indices: List[int],
    displacement: List[Optional[float]]
) -> Dict[int, List[Optional[float]]]:
    """
    Create displacement boundary conditions for specified nodes.
    
    Args:
        node_indices: List of node indices to apply displacement to
        displacement: List of 6 values [ux, uy, uz, rx, ry, rz] where None means unconstrained
        
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


def apply_boundary_conditions(
    nodes: np.ndarray,
    elements: np.ndarray,
    constraints: Dict[Union[Surface, str], Union[str, List[Optional[float]], Callable]],
    **kwargs
) -> Tuple[Dict[int, List[Optional[float]]], np.ndarray]:
    """
    Apply boundary conditions to the FE model using the enhanced system.
    
    Args:
        nodes: Array of node coordinates (N x 3)
        elements: Array of element connectivity (E x 8)
        constraints: Dictionary mapping surface identifiers or custom functions to constraints
        **kwargs: Additional parameters for customization
            
    Returns:
        Tuple containing:
            - dofs: Dictionary mapping node indices to prescribed displacements
            - forces: Array of nodal forces (N x 6)
    """
    # Initialize DOFs dictionary and forces array
    dofs = {}
    forces = np.zeros((nodes.shape[0], 6))  # [Fx, Fy, Fz, Mx, My, Mz]
    
    # Process each constraint
    for key, value in constraints.items():
        if isinstance(key, Surface):
            # Simple mode: Surface enumeration
            node_indices = identify_surface_nodes(nodes, key)
            
            if value == "fix":
                # Fix all DOFs
                bc = BoundaryCondition(node_indices, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
                bc.apply(dofs)
            elif isinstance(value, list):
                # Apply displacement vector
                bc = BoundaryCondition(node_indices, value)
                bc.apply(dofs)
            else:
                raise ValueError(f"Unsupported constraint value for surface {key}: {value}")
        
        elif key == "custom" and callable(value):
            # Expert mode: Custom function
            custom_dofs = value(nodes, elements)
            dofs.update(custom_dofs)
        
        else:
            raise ValueError(f"Unsupported constraint key: {key}")
    
    return dofs, forces


