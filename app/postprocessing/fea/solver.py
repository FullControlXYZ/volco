"""
Linear static solver for FEA.

This module provides functionality to solve linear static FEA problems
using NumPy and SciPy.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Any
from scipy.sparse import coo_matrix, lil_matrix
from scipy.sparse.linalg import spsolve


def solve_static_problem(
    nodes: np.ndarray,
    elements: np.ndarray,
    dofs: Dict[int, List[Optional[float]]],
    forces: np.ndarray,
    material_properties: Dict[str, float]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Solve a linear static FEA problem.
    
    Args:
        nodes: Array of node coordinates (N x 3)
        elements: Array of element connectivity (E x 8)
        dofs: Dictionary mapping node indices to prescribed displacements
        forces: Array of nodal forces (N x 6)
        material_properties: Dictionary containing material properties
            - 'young_modulus': Young's modulus in MPa
            - 'poisson_ratio': Poisson's ratio
            
    Returns:
        Tuple containing:
            - displacements: Nodal displacements (N x 3)
            - stresses: Element stresses (E x 6)
            - strains: Element strains (E x 6)
            - von_mises: von Mises stress for each element (E)
    """
    # Extract material properties
    E = material_properties['young_modulus']
    nu = material_properties['poisson_ratio']
    
    # Number of nodes and elements
    num_nodes = nodes.shape[0]
    num_elements = elements.shape[0]
    
    # Assemble global stiffness matrix and force vector
    K_global, f_global = assemble_system(nodes, elements, forces, E, nu)
    
    # Apply boundary conditions
    K_modified, f_modified = apply_boundary_conditions(K_global, f_global, dofs, num_nodes)
    
    # Solve the system of equations
    u_global = solve_system(K_modified, f_modified)
    
    # Reshape the solution to get nodal displacements
    displacements = u_global.reshape(num_nodes, 3)
    
    # Calculate strains, stresses, and von Mises stress
    strains, stresses, von_mises = calculate_stresses_and_strains(
        nodes, elements, displacements, E, nu
    )
    
    return displacements, stresses, strains, von_mises


def assemble_system(
    nodes: np.ndarray,
    elements: np.ndarray,
    forces: np.ndarray,
    E: float,
    nu: float
) -> Tuple[lil_matrix, np.ndarray]:
    """
    Assemble the global stiffness matrix and force vector.
    
    Args:
        nodes: Array of node coordinates
        elements: Array of element connectivity
        forces: Array of nodal forces
        E: Young's modulus
        nu: Poisson's ratio
        
    Returns:
        Tuple containing:
            - K_global: Global stiffness matrix
            - f_global: Global force vector
    """
    # Number of nodes and DOFs
    num_nodes = nodes.shape[0]
    num_dofs = num_nodes * 3  # 3 DOFs per node (x, y, z displacements)
    
    # Initialize global stiffness matrix and force vector
    K_global = lil_matrix((num_dofs, num_dofs))
    f_global = np.zeros(num_dofs)
    
    # Populate force vector from nodal forces
    for i in range(num_nodes):
        f_global[i*3:(i+1)*3] = forces[i, :3]  # Only translational forces
    
    # For each element, calculate element stiffness matrix and assemble into global matrix
    for el_idx, element in enumerate(elements):
        # Get element nodes
        el_nodes = nodes[element]
        
        # Calculate element stiffness matrix
        K_el = calculate_element_stiffness(el_nodes, E, nu)
        
        # Assemble element stiffness matrix into global stiffness matrix
        for i in range(8):  # 8 nodes per hexahedral element
            for j in range(8):
                # Get global indices for the nodes
                ni = element[i]
                nj = element[j]
                
                # For each DOF (x, y, z)
                for di in range(3):
                    for dj in range(3):
                        # Global DOF indices
                        gi = ni * 3 + di
                        gj = nj * 3 + dj
                        
                        # Add element stiffness to global stiffness
                        K_global[gi, gj] += K_el[i*3+di, j*3+dj]
    
    return K_global, f_global


def calculate_element_stiffness(
    el_nodes: np.ndarray,
    E: float,
    nu: float
) -> np.ndarray:
    """
    Calculate the stiffness matrix for a hexahedral element.
    
    This is a simplified implementation using the direct stiffness method
    for a linear hexahedral element.
    
    Args:
        el_nodes: Coordinates of the element nodes (8 x 3)
        E: Young's modulus
        nu: Poisson's ratio
        
    Returns:
        Element stiffness matrix (24 x 24)
    """
    # Calculate element dimensions
    dx = np.max(el_nodes[:, 0]) - np.min(el_nodes[:, 0])
    dy = np.max(el_nodes[:, 1]) - np.min(el_nodes[:, 1])
    dz = np.max(el_nodes[:, 2]) - np.min(el_nodes[:, 2])
    
    # Element volume
    V = dx * dy * dz
    
    # Material matrix (plane stress)
    D = np.zeros((6, 6))
    c = E / ((1 + nu) * (1 - 2*nu))
    
    # Fill material matrix
    D[0, 0] = D[1, 1] = D[2, 2] = c * (1 - nu)
    D[0, 1] = D[0, 2] = D[1, 0] = D[1, 2] = D[2, 0] = D[2, 1] = c * nu
    D[3, 3] = D[4, 4] = D[5, 5] = c * (1 - 2*nu) / 2
    
    # Initialize element stiffness matrix
    K_el = np.zeros((24, 24))  # 8 nodes x 3 DOFs
    
    # Integration points for 2x2x2 Gaussian quadrature
    gauss_points = np.array([
        [-1/np.sqrt(3), -1/np.sqrt(3), -1/np.sqrt(3)],
        [1/np.sqrt(3), -1/np.sqrt(3), -1/np.sqrt(3)],
        [-1/np.sqrt(3), 1/np.sqrt(3), -1/np.sqrt(3)],
        [1/np.sqrt(3), 1/np.sqrt(3), -1/np.sqrt(3)],
        [-1/np.sqrt(3), -1/np.sqrt(3), 1/np.sqrt(3)],
        [1/np.sqrt(3), -1/np.sqrt(3), 1/np.sqrt(3)],
        [-1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)],
        [1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)]
    ])
    
    # Weight for each integration point
    weight = 1.0
    
    # For each integration point
    for gp in gauss_points:
        xi, eta, zeta = gp
        
        # Calculate shape functions and derivatives
        N, dN = shape_functions(xi, eta, zeta)
        
        # Calculate Jacobian matrix
        J = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                for k in range(8):
                    J[i, j] += dN[k, i] * el_nodes[k, j]
        
        # Calculate determinant of Jacobian
        detJ = np.linalg.det(J)
        
        # Calculate inverse of Jacobian
        Jinv = np.linalg.inv(J)
        
        # Calculate B matrix (strain-displacement matrix)
        B = np.zeros((6, 24))
        for i in range(8):
            # Derivatives of shape functions with respect to x, y, z
            dNdx = Jinv[0, 0] * dN[i, 0] + Jinv[0, 1] * dN[i, 1] + Jinv[0, 2] * dN[i, 2]
            dNdy = Jinv[1, 0] * dN[i, 0] + Jinv[1, 1] * dN[i, 1] + Jinv[1, 2] * dN[i, 2]
            dNdz = Jinv[2, 0] * dN[i, 0] + Jinv[2, 1] * dN[i, 1] + Jinv[2, 2] * dN[i, 2]
            
            # Fill B matrix
            B[0, i*3] = dNdx      # du/dx
            B[1, i*3+1] = dNdy    # dv/dy
            B[2, i*3+2] = dNdz    # dw/dz
            B[3, i*3] = dNdy      # du/dy
            B[3, i*3+1] = dNdx    # dv/dx
            B[4, i*3+1] = dNdz    # dv/dz
            B[4, i*3+2] = dNdy    # dw/dy
            B[5, i*3] = dNdz      # du/dz
            B[5, i*3+2] = dNdx    # dw/dx
        
        # Calculate contribution to element stiffness matrix
        K_el += B.T @ D @ B * detJ * weight
    
    return K_el


def shape_functions(xi: float, eta: float, zeta: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate shape functions and their derivatives for a hexahedral element.
    
    Args:
        xi, eta, zeta: Natural coordinates (-1 to 1)
        
    Returns:
        Tuple containing:
            - N: Shape functions (8)
            - dN: Derivatives of shape functions with respect to xi, eta, zeta (8 x 3)
    """
    # Shape functions
    N = np.zeros(8)
    N[0] = (1 - xi) * (1 - eta) * (1 - zeta) / 8
    N[1] = (1 + xi) * (1 - eta) * (1 - zeta) / 8
    N[2] = (1 + xi) * (1 + eta) * (1 - zeta) / 8
    N[3] = (1 - xi) * (1 + eta) * (1 - zeta) / 8
    N[4] = (1 - xi) * (1 - eta) * (1 + zeta) / 8
    N[5] = (1 + xi) * (1 - eta) * (1 + zeta) / 8
    N[6] = (1 + xi) * (1 + eta) * (1 + zeta) / 8
    N[7] = (1 - xi) * (1 + eta) * (1 + zeta) / 8
    
    # Derivatives of shape functions
    dN = np.zeros((8, 3))
    
    # dN/dxi
    dN[0, 0] = -(1 - eta) * (1 - zeta) / 8
    dN[1, 0] = (1 - eta) * (1 - zeta) / 8
    dN[2, 0] = (1 + eta) * (1 - zeta) / 8
    dN[3, 0] = -(1 + eta) * (1 - zeta) / 8
    dN[4, 0] = -(1 - eta) * (1 + zeta) / 8
    dN[5, 0] = (1 - eta) * (1 + zeta) / 8
    dN[6, 0] = (1 + eta) * (1 + zeta) / 8
    dN[7, 0] = -(1 + eta) * (1 + zeta) / 8
    
    # dN/deta
    dN[0, 1] = -(1 - xi) * (1 - zeta) / 8
    dN[1, 1] = -(1 + xi) * (1 - zeta) / 8
    dN[2, 1] = (1 + xi) * (1 - zeta) / 8
    dN[3, 1] = (1 - xi) * (1 - zeta) / 8
    dN[4, 1] = -(1 - xi) * (1 + zeta) / 8
    dN[5, 1] = -(1 + xi) * (1 + zeta) / 8
    dN[6, 1] = (1 + xi) * (1 + zeta) / 8
    dN[7, 1] = (1 - xi) * (1 + zeta) / 8
    
    # dN/dzeta
    dN[0, 2] = -(1 - xi) * (1 - eta) / 8
    dN[1, 2] = -(1 + xi) * (1 - eta) / 8
    dN[2, 2] = -(1 + xi) * (1 + eta) / 8
    dN[3, 2] = -(1 - xi) * (1 + eta) / 8
    dN[4, 2] = (1 - xi) * (1 - eta) / 8
    dN[5, 2] = (1 + xi) * (1 - eta) / 8
    dN[6, 2] = (1 + xi) * (1 + eta) / 8
    dN[7, 2] = (1 - xi) * (1 + eta) / 8
    
    return N, dN


def apply_boundary_conditions(
    K_global: lil_matrix,
    f_global: np.ndarray,
    dofs: Dict[int, List[Optional[float]]],
    num_nodes: int
) -> Tuple[lil_matrix, np.ndarray]:
    """
    Apply boundary conditions to the global stiffness matrix and force vector.
    
    Args:
        K_global: Global stiffness matrix
        f_global: Global force vector
        dofs: Dictionary mapping node indices to prescribed displacements
        num_nodes: Number of nodes
        
    Returns:
        Tuple containing:
            - K_modified: Modified stiffness matrix
            - f_modified: Modified force vector
    """
    # Make a copy of the stiffness matrix and force vector
    K_modified = K_global.copy()
    f_modified = f_global.copy()
    
    # For each node with prescribed DOFs
    for node_idx, node_dofs in dofs.items():
        # For each DOF (x, y, z)
        for dof_idx in range(3):
            # Check if this DOF is prescribed
            if node_dofs[dof_idx] is not None:
                # Global DOF index
                global_dof = node_idx * 3 + dof_idx
                
                # Prescribed displacement value
                prescribed_value = node_dofs[dof_idx]
                
                # Modify force vector: f = f - K * u_prescribed
                f_modified -= K_modified[:, global_dof].toarray().flatten() * prescribed_value
                
                # Modify stiffness matrix: set row and column to zero, diagonal to 1
                K_modified[global_dof, :] = 0
                K_modified[:, global_dof] = 0
                K_modified[global_dof, global_dof] = 1
                
                # Set force vector entry to prescribed value
                f_modified[global_dof] = prescribed_value
    
    return K_modified, f_modified


def solve_system(K: lil_matrix, f: np.ndarray) -> np.ndarray:
    """
    Solve the system of equations K * u = f.
    
    Args:
        K: Stiffness matrix
        f: Force vector
        
    Returns:
        Solution vector u
    """
    # Convert to CSR format for efficient solving
    K_csr = K.tocsr()
    
    # Solve the system
    u = spsolve(K_csr, f)
    
    return u


def calculate_stresses_and_strains(
    nodes: np.ndarray,
    elements: np.ndarray,
    displacements: np.ndarray,
    E: float,
    nu: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate strains, stresses, and von Mises stress for each element.
    
    Args:
        nodes: Array of node coordinates
        elements: Array of element connectivity
        displacements: Nodal displacements
        E: Young's modulus
        nu: Poisson's ratio
        
    Returns:
        Tuple containing:
            - strains: Element strains (E x 6)
            - stresses: Element stresses (E x 6)
            - von_mises: von Mises stress for each element (E)
    """
    # Number of elements
    num_elements = elements.shape[0]
    
    # Initialize arrays for results
    strains = np.zeros((num_elements, 6))
    stresses = np.zeros((num_elements, 6))
    von_mises = np.zeros(num_elements)
    
    # Material matrix
    D = np.zeros((6, 6))
    c = E / ((1 + nu) * (1 - 2*nu))
    
    # Fill material matrix
    D[0, 0] = D[1, 1] = D[2, 2] = c * (1 - nu)
    D[0, 1] = D[0, 2] = D[1, 0] = D[1, 2] = D[2, 0] = D[2, 1] = c * nu
    D[3, 3] = D[4, 4] = D[5, 5] = c * (1 - 2*nu) / 2
    
    # For each element
    for el_idx, element in enumerate(elements):
        # Get element nodes
        el_nodes = nodes[element]
        
        # Get element displacements
        el_displacements = np.zeros(24)
        for i in range(8):
            node_idx = element[i]
            el_displacements[i*3:(i+1)*3] = displacements[node_idx]
        
        # Calculate strain at the center of the element
        xi = eta = zeta = 0.0
        N, dN = shape_functions(xi, eta, zeta)
        
        # Calculate Jacobian matrix
        J = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                for k in range(8):
                    J[i, j] += dN[k, i] * el_nodes[k, j]
        
        # Calculate inverse of Jacobian
        Jinv = np.linalg.inv(J)
        
        # Calculate B matrix (strain-displacement matrix)
        B = np.zeros((6, 24))
        for i in range(8):
            # Derivatives of shape functions with respect to x, y, z
            dNdx = Jinv[0, 0] * dN[i, 0] + Jinv[0, 1] * dN[i, 1] + Jinv[0, 2] * dN[i, 2]
            dNdy = Jinv[1, 0] * dN[i, 0] + Jinv[1, 1] * dN[i, 1] + Jinv[1, 2] * dN[i, 2]
            dNdz = Jinv[2, 0] * dN[i, 0] + Jinv[2, 1] * dN[i, 1] + Jinv[2, 2] * dN[i, 2]
            
            # Fill B matrix
            B[0, i*3] = dNdx      # du/dx
            B[1, i*3+1] = dNdy    # dv/dy
            B[2, i*3+2] = dNdz    # dw/dz
            B[3, i*3] = dNdy      # du/dy
            B[3, i*3+1] = dNdx    # dv/dx
            B[4, i*3+1] = dNdz    # dv/dz
            B[4, i*3+2] = dNdy    # dw/dy
            B[5, i*3] = dNdz      # du/dz
            B[5, i*3+2] = dNdx    # dw/dx
        
        # Calculate strain: ε = B * u
        strain = B @ el_displacements
        
        # Calculate stress: σ = D * ε
        stress = D @ strain
        
        # Store results
        strains[el_idx] = strain
        stresses[el_idx] = stress
        
        # Calculate von Mises stress
        s1, s2, s3 = stress[0], stress[1], stress[2]  # Normal stresses
        s4, s5, s6 = stress[3], stress[4], stress[5]  # Shear stresses
        
        von_mises[el_idx] = np.sqrt(0.5 * ((s1 - s2)**2 + (s2 - s3)**2 + (s3 - s1)**2 + 6*(s4**2 + s5**2 + s6**2)))
    
    return strains, stresses, von_mises