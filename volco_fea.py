"""
VOLCO FEA Module

This module provides a simplified interface to the Finite Element Analysis (FEA)
functionality in VOLCO. It allows users to perform structural analysis on voxel models,
visualize results, and save/load analysis data.

Basic usage:
    from volco_fea import analyze_voxel_matrix, Surface
    
    # Define boundary conditions
    boundary_conditions = {
        'constraints': {
            Surface.MINUS_Z: "fix",  # Fix bottom surface
            Surface.PLUS_Z: [None, None, -0.1, None, None, None]  # Apply displacement on top
        }
    }
    
    # Run analysis
    results = analyze_voxel_matrix(
        voxel_matrix=voxel_matrix,
        voxel_size=voxel_size,
        boundary_conditions=boundary_conditions
    )
"""

# Re-export core FEA functionality
from app.postprocessing.fea.core import analyze_voxel_matrix, load_fea_results
from app.postprocessing.fea.viz import visualize_fea, export_visualization, visualize_voxel_matrix
from app.postprocessing.fea.boundary import (
    Surface,
    select_nodes_by_predicate,
    select_nodes_in_box,
    select_nodes_on_plane,
    select_nodes_by_position
)
from app.postprocessing.fea.io import save_results, load_results

# Define what gets imported with "from volco_fea import *"
__all__ = [
    # Core functionality
    'analyze_voxel_matrix',
    'load_fea_results',
    
    # Visualization
    'visualize_fea',
    'export_visualization',
    'visualize_voxel_matrix',
    
    # Boundary conditions
    'Surface',
    'select_nodes_by_predicate',
    'select_nodes_in_box',
    'select_nodes_on_plane',
    'select_nodes_by_position',
    
    # I/O operations
    'save_results',
    'load_results'
]