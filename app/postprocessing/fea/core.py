"""
Core FEA implementation for VOLCO.

This module provides the main entry point for FEA analysis of voxel models.
"""

import numpy as np
import os
from typing import Dict, Any, Optional, Tuple, Union

from .mesh import generate_mesh, check_continuity
from .solver import solve_static_problem
from .boundary import apply_boundary_conditions, Surface
from .viz import visualize_fea
from .io import save_results


def analyze_voxel_matrix(
    voxel_matrix: np.ndarray,
    voxel_size: float,
    material_properties: Optional[Dict[str, float]] = None,
    boundary_conditions: Optional[Dict[str, Any]] = None,
    **options
) -> Dict[str, Any]:
    """
    Perform finite element analysis on a voxel matrix.
    
    This function serves as the main entry point for FEA analysis. It takes a voxel matrix
    (typically from a VOLCO simulation), converts it to a FE mesh, applies boundary conditions,
    solves the static problem, and returns the results.
    
    Args:
        voxel_matrix: 3D numpy array representing the voxel model (1 = material, 0 = void)
        voxel_size: Size of each voxel in mm
        material_properties: Dictionary containing material properties:
            - 'young_modulus': Young's modulus in MPa (default: 2000 for typical PLA)
            - 'poisson_ratio': Poisson's ratio (default: 0.3 for typical PLA)
        boundary_conditions: Dictionary specifying boundary conditions (if None, defaults are used)
        **options: Additional options for the analysis:
            - 'visualization': Whether to generate visualization (default: True)
            - 'result_type': Type of result to visualize ('displacement', 'strain', 'stress', 'von_mises')
            - 'scale_factor': Factor to scale displacements for visualization (default: 1.0)
            - 'save_results': Whether to save results to a file (default: False)
            - 'save_path': Path to save results (default: 'Results_volco/fea/results')
            - 'save_format': Format to save results ('pickle', 'json', or 'hdf5') (default: 'pickle')
            - 'include_visualization': Whether to include visualization in saved file (default: False)
            
    Returns:
        Dictionary containing analysis results:
            - 'nodes': Node coordinates
            - 'elements': Element connectivity
            - 'displacements': Nodal displacements
            - 'stresses': Element stresses
            - 'strains': Element strains
            - 'von_mises': von Mises stress for each element
            - 'max_displacement': Maximum displacement magnitude
            - 'max_von_mises': Maximum von Mises stress
            - 'visualization': Visualization object (if visualization=True)
            - 'saved_file': Path to saved results file (if save_results=True)
    """
    # Set default material properties if not provided
    if material_properties is None:
        material_properties = {
            'young_modulus': 2000.0,  # MPa (typical for PLA)
            'poisson_ratio': 0.3      # Typical for PLA
        }
    
    # Check if the voxel matrix represents a single continuous body
    if not check_continuity(voxel_matrix):
        raise ValueError("The voxel matrix does not represent a single continuous body")
    
    # Generate FE mesh from voxel matrix
    nodes, elements = generate_mesh(voxel_matrix, voxel_size)
    
    # Apply boundary conditions (must be supplied)
    if boundary_conditions is None:
        raise ValueError("Boundary conditions must be supplied. "
                        "Please use the format: {'constraints': {Surface.MINUS_Z: 'fix', ...}}")
    
    # Check if constraints are specified in the new format
    if 'constraints' in boundary_conditions:
        constraints = boundary_conditions.get('constraints', {})
        dofs, forces = apply_boundary_conditions(nodes, elements, constraints)
    else:
        raise ValueError("Boundary conditions must be specified using the 'constraints' key. "
                        "Please use the format: {'constraints': {Surface.MINUS_Z: 'fix', ...}}")
    
    # Solve the static problem
    displacements, stresses, strains, von_mises = solve_static_problem(
        nodes, elements, dofs, forces, material_properties
    )
    
    # Prepare results dictionary
    results = {
        'nodes': nodes,
        'elements': elements,
        'displacements': displacements,
        'stresses': stresses,
        'strains': strains,
        'von_mises': von_mises,
        'max_displacement': np.max(np.sqrt(np.sum(displacements**2, axis=1))),
        'max_von_mises': np.max(von_mises)
    }
    
    # Generate visualization if requested
    visualization_enabled = options.get('visualization', True)
    if visualization_enabled:
        result_type = options.get('result_type', 'von_mises')
        scale_factor = options.get('scale_factor', 1.0)
        show_undeformed = options.get('show_undeformed', False)
        original_opacity = options.get('original_opacity', 0.3)
        viz = visualize_fea(
            nodes, elements, displacements, von_mises,
            result_type=result_type,
            scale_factor=scale_factor,
            show_undeformed=show_undeformed,
            original_opacity=original_opacity
        )
        results['visualization'] = viz
    
    # Save results if requested
    if options.get('save_results', False):
        save_path = options.get('save_path', 'Results_volco/fea/results')
        save_format = options.get('save_format', 'pickle')
        include_visualization = options.get('include_visualization', False)
        
        # Ensure the save path has the correct extension
        if not save_path.endswith(f'.{save_format}'):
            if save_format == 'pickle':
                save_path += '.pkl'
            else:
                save_path += f'.{save_format}'
        
        # Create metadata with material properties and boundary conditions
        metadata = {
            'voxel_size': voxel_size,
            'material_properties': material_properties,
            'boundary_conditions': boundary_conditions,
            'options': {k: v for k, v in options.items() if k not in ['save_results', 'save_path', 'save_format']}
        }
        
        # Save the results
        saved_file = save_results(
            results=results,
            filename=save_path,
            format=save_format,
            include_visualization=include_visualization,
            metadata=metadata
        )
        
        results['saved_file'] = saved_file
    
    return results


def load_fea_results(file_path: str) -> Dict[str, Any]:
    """
    Load FEA results from a file.
    
    Args:
        file_path: Path to the file containing saved FEA results
        
    Returns:
        Dictionary containing the loaded FEA results
    """
    from .io import load_results
    
    # Load the results
    data = load_results(file_path)
    
    # Return just the results part
    return data['results']