"""
Visualization module for FEA results.

This module provides functionality to visualize FEA results such as
displacements, strains, and stresses using Plotly. It features an optimized
mesh visualization approach (the "volco mesh visualization method") that only
renders visible faces for improved performance with large models.
"""

import numpy as np
from typing import Dict, Any, Optional, Tuple, List, Union
import plotly.graph_objects as go


def visualize_fea(
    nodes: np.ndarray,
    elements: np.ndarray,
    displacements: np.ndarray,
    von_mises: np.ndarray,
    result_type: str = 'von_mises',
    scale_factor: float = 1.0,
    show_undeformed: bool = True,
    original_opacity: float = 0.3
) -> Any:
    """
    Unified visualization function for FEA results using Plotly with optimized mesh rendering.
    
    This function combines the functionality of the previous visualization functions,
    allowing for visualization of either displacement or von Mises stress, with the
    option to show or hide the undeformed mesh alongside the deformed mesh.
    
    The function uses the "volco mesh visualization method" which:
    - Only renders visible faces of the mesh (faces that are exposed to the outside)
    - Creates a 3D grid representation of the voxel structure to determine visibility
    - Triangulates quadrilateral faces for better rendering
    - Significantly improves performance for large models by reducing the number of rendered faces
    
    Args:
        nodes: Array of node coordinates with shape (n_nodes, 3)
        elements: Array of element connectivity with shape (n_elements, 8) for hexahedral elements
        displacements: Nodal displacements with shape (n_nodes, 3)
        von_mises: von Mises stress for each element with shape (n_elements,)
        result_type: Type of result to visualize ('displacement', 'von_mises')
        scale_factor: Factor to scale displacements for visualization
        show_undeformed: Whether to show the original undeformed mesh alongside the deformed mesh
        original_opacity: Opacity of the original mesh when shown (0.0-1.0)
        
    Returns:
        Plotly figure with interactive 3D visualization of the FEA results
    """
    
    # Calculate deformed node positions
    deformed_nodes = nodes + displacements * scale_factor
    
    # Prepare data for visualization
    if result_type == 'displacement':
        # Calculate displacement magnitude
        disp_mag = np.sqrt(np.sum(displacements**2, axis=1))
        node_colors = disp_mag
        vmin, vmax = np.min(disp_mag), np.max(disp_mag)
        title = f'Displacement Magnitude (max = {vmax:.4f})'
        colorbar_title = 'Displacement'
    else:  # von_mises
        # Map element von Mises stress to nodes
        node_colors = np.zeros(nodes.shape[0])
        node_counts = np.zeros(nodes.shape[0])
        
        for i, element in enumerate(elements):
            for node_idx in element:
                node_colors[node_idx] += von_mises[i]
                node_counts[node_idx] += 1
        
        # Average the values
        node_colors = np.divide(node_colors, node_counts, where=node_counts>0)
        vmin, vmax = np.min(von_mises), np.max(von_mises)
        title = f'von Mises Stress (max = {vmax:.4f} MPa)'
        colorbar_title = 'Stress (MPa)'
    
    # Define faces for each hexahedral element and the direction to check for each face
    hex_faces = [
        ([0, 1, 2, 3], (0, 0, -1)),  # Bottom face (negative z)
        ([4, 5, 6, 7], (0, 0, 1)),   # Top face (positive z)
        ([0, 1, 5, 4], (0, -1, 0)),  # Back face (negative y)
        ([2, 3, 7, 6], (0, 1, 0)),   # Front face (positive y)
        ([0, 3, 7, 4], (-1, 0, 0)),  # Left face (negative x)
        ([1, 2, 6, 5], (1, 0, 0))    # Right face (positive x)
    ]
    
    # Create a 3D grid to represent the voxel-like structure
    # First, determine the grid dimensions by finding the min and max coordinates
    min_coords = np.min(nodes, axis=0)
    max_coords = np.max(nodes, axis=0)
    
    # Estimate the voxel size based on the average distance between adjacent nodes
    sample_element = elements[0]
    sample_nodes = nodes[sample_element]
    voxel_size = np.mean([
        np.linalg.norm(sample_nodes[1] - sample_nodes[0]),  # x-direction
        np.linalg.norm(sample_nodes[3] - sample_nodes[0]),  # y-direction
        np.linalg.norm(sample_nodes[4] - sample_nodes[0])   # z-direction
    ])
    
    # Calculate grid dimensions
    grid_dims = np.ceil((max_coords - min_coords) / voxel_size).astype(int) + 1
    
    # Create an empty grid
    grid = np.zeros(grid_dims, dtype=bool)
    
    # Map each element to a grid cell
    element_to_grid = {}
    grid_to_element = {}
    
    for el_idx, element in enumerate(elements):
        # Calculate the center of the element
        element_center = np.mean(nodes[element], axis=0)
        
        # Convert to grid coordinates
        grid_coords = np.floor((element_center - min_coords) / voxel_size).astype(int)
        
        # Ensure we're within bounds
        grid_coords = np.clip(grid_coords, 0, grid_dims - 1)
        
        # Mark this cell as filled
        grid_tuple = tuple(grid_coords)
        grid[grid_tuple] = True
        
        # Store the mapping
        element_to_grid[el_idx] = grid_tuple
        grid_to_element[grid_tuple] = el_idx
    
    # Create the figure
    fig = go.Figure()
    
    # If showing the undeformed mesh, add it first
    if show_undeformed:
        # Prepare data for original mesh (transparent grey)
        orig_i_faces = []
        orig_j_faces = []
        orig_k_faces = []
        
        # Create a mapping from original node indices to new indices for original mesh
        orig_node_map = {}
        orig_new_nodes = []
        orig_idx = 0
        
        # For each element
        for el_idx, element in enumerate(elements):
            grid_coords = element_to_grid[el_idx]
            
            # Check each face to see if it should be rendered
            for face_idx, (face_vertices, direction) in enumerate(hex_faces):
                # Calculate the neighboring cell coordinates
                neighbor_coords = (
                    grid_coords[0] + direction[0],
                    grid_coords[1] + direction[1],
                    grid_coords[2] + direction[2]
                )
                
                # Check if the neighboring cell is outside the grid or empty
                if (neighbor_coords[0] < 0 or neighbor_coords[0] >= grid_dims[0] or
                    neighbor_coords[1] < 0 or neighbor_coords[1] >= grid_dims[1] or
                    neighbor_coords[2] < 0 or neighbor_coords[2] >= grid_dims[2] or
                    not grid[neighbor_coords]):
                    
                    # This face is visible, render it for original mesh
                    face_nodes = [element[v] for v in face_vertices]
                    
                    # Add triangles (triangulate the quadrilateral face)
                    triangles = [[0, 1, 2], [0, 2, 3]]  # Two triangles per face
                    
                    for tri in triangles:
                        # Get the nodes for this triangle
                        tri_nodes = [face_nodes[t] for t in tri]
                        
                        # Add the nodes to the original mesh
                        for node_idx in tri_nodes:
                            if node_idx not in orig_node_map:
                                orig_node_map[node_idx] = orig_idx
                                orig_new_nodes.append(nodes[node_idx])
                                orig_idx += 1
                        
                        # Add the face indices for original mesh
                        orig_i_faces.append(orig_node_map[tri_nodes[0]])
                        orig_j_faces.append(orig_node_map[tri_nodes[1]])
                        orig_k_faces.append(orig_node_map[tri_nodes[2]])
        
        # Convert to numpy arrays
        orig_new_nodes = np.array(orig_new_nodes)
        
        # Add the original mesh (transparent grey)
        fig.add_trace(
            go.Mesh3d(
                x=orig_new_nodes[:, 0],
                y=orig_new_nodes[:, 1],
                z=orig_new_nodes[:, 2],
                i=orig_i_faces,
                j=orig_j_faces,
                k=orig_k_faces,
                color='grey',
                opacity=original_opacity,
                name='Original Mesh'
            )
        )
    
    # Prepare data for deformed mesh
    def_i_faces = []
    def_j_faces = []
    def_k_faces = []
    
    # Create a mapping from original node indices to new indices for deformed mesh
    def_node_map = {}
    def_new_nodes = []
    def_new_colors = []
    def_idx = 0
    
    # For each element
    for el_idx, element in enumerate(elements):
        grid_coords = element_to_grid[el_idx]
        
        # Check each face to see if it should be rendered
        for face_idx, (face_vertices, direction) in enumerate(hex_faces):
            # Calculate the neighboring cell coordinates
            neighbor_coords = (
                grid_coords[0] + direction[0],
                grid_coords[1] + direction[1],
                grid_coords[2] + direction[2]
            )
            
            # Check if the neighboring cell is outside the grid or empty
            if (neighbor_coords[0] < 0 or neighbor_coords[0] >= grid_dims[0] or
                neighbor_coords[1] < 0 or neighbor_coords[1] >= grid_dims[1] or
                neighbor_coords[2] < 0 or neighbor_coords[2] >= grid_dims[2] or
                not grid[neighbor_coords]):
                
                # This face is visible, render it for deformed mesh
                face_nodes = [element[v] for v in face_vertices]
                
                # Add triangles (triangulate the quadrilateral face)
                triangles = [[0, 1, 2], [0, 2, 3]]  # Two triangles per face
                
                for tri in triangles:
                    # Get the nodes for this triangle
                    tri_nodes = [face_nodes[t] for t in tri]
                    
                    # Add the nodes and colors to the deformed mesh
                    for node_idx in tri_nodes:
                        if node_idx not in def_node_map:
                            def_node_map[node_idx] = def_idx
                            def_new_nodes.append(deformed_nodes[node_idx])
                            def_new_colors.append(node_colors[node_idx])
                            def_idx += 1
                    
                    # Add the face indices for deformed mesh
                    def_i_faces.append(def_node_map[tri_nodes[0]])
                    def_j_faces.append(def_node_map[tri_nodes[1]])
                    def_k_faces.append(def_node_map[tri_nodes[2]])
    
    # Convert to numpy arrays
    def_new_nodes = np.array(def_new_nodes)
    def_new_colors = np.array(def_new_colors)
    
    # Add the deformed mesh (colored by result type)
    fig.add_trace(
        go.Mesh3d(
            x=def_new_nodes[:, 0],
            y=def_new_nodes[:, 1],
            z=def_new_nodes[:, 2],
            i=def_i_faces,
            j=def_j_faces,
            k=def_k_faces,
            intensity=def_new_colors,
            colorscale='Jet',
            colorbar=dict(title=colorbar_title),
            cmin=vmin,
            cmax=vmax,
            showscale=True,
            name='Deformed Mesh'
        )
    )
    
    # Update layout
    if show_undeformed:
        title = f'{title} (with original mesh)'
    
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z',
            aspectmode='data'
        ),
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01
        )
    )
    
    return fig
    


def export_visualization(
    fig: Any,
    filename: str
) -> None:
    """
    Export visualization to a file.
    
    Args:
        fig: Visualization figure (plotly figure)
        filename: Output filename
    """
    fig.write_html(filename)


def visualize_voxel_matrix(
    voxel_matrix: np.ndarray,
    voxel_size: float
) -> go.Scatter3d:
    """
    Visualize a voxel matrix using Plotly.
    
    Args:
        voxel_matrix: 3D numpy array representing the voxel model
        voxel_size: Size of each voxel
        
    Returns:
        Plotly Scatter3d trace
    """
    # Find indices of non-zero voxels
    voxels = np.argwhere(voxel_matrix > 0)
    
    # Convert to physical coordinates
    x = voxels[:, 0] * voxel_size
    y = voxels[:, 1] * voxel_size
    z = voxels[:, 2] * voxel_size
    
    # Create a colormap based on z-coordinate
    colors = z / np.max(z) if np.max(z) > 0 else z
    
    # Create a scatter3d trace for voxels
    trace = go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode='markers',
        marker=dict(
            size=voxel_size*10,
            color=colors,
            colorscale='Viridis',
            opacity=0.8
        )
    )
    
    return trace

