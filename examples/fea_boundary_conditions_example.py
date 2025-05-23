"""
Example script demonstrating the enhanced boundary condition system for VOLCO FEA.

This script shows how to:
1. Use Simple Mode with Surface enumerations
2. Use Expert Mode with custom functions
3. Use the functional utilities for node selection
"""

import os
import numpy as np
from volco import run_simulation
from volco_fea import analyze_voxel_matrix, Surface, select_nodes_by_predicate, select_nodes_in_box, export_visualization


def run_simple_mode_example(voxel_matrix, voxel_size):
    """Run FEA analysis using Simple Mode boundary conditions."""
    print("\n=== Simple Mode Example ===")
    
    # Create output directory if it doesn't exist
    os.makedirs("Results_volco/fea/boundary_examples", exist_ok=True)
    
    # Define material properties for PLA
    material_properties = {
        'young_modulus': 2000.0,  # MPa (typical for PLA)
        'poisson_ratio': 0.3      # Typical for PLA
    }
    
    # Define boundary conditions using Simple Mode
    # Get model dimensions - use the z-dimension of the voxel matrix
    model_height = voxel_matrix.shape[2] * voxel_size
    
    # Calculate displacement as 1% of model height (similar to default)
    displacement_magnitude = model_height * 0.01
    
    boundary_conditions = {
        'constraints': {
            Surface.MINUS_Z: "fix",  # Fix bottom surface
            Surface.PLUS_Z: [None, None, -displacement_magnitude, None, None, None]  # Apply 1% compression on top
        }
    }
    
    # Run the analysis
    results = analyze_voxel_matrix(
        voxel_matrix=voxel_matrix,
        voxel_size=voxel_size,
        material_properties=material_properties,
        boundary_conditions=boundary_conditions,
        visualization=True,
        result_type='von_mises',
        scale_factor=10.0  # Exaggerate deformation for visualization
    )
    
    print("FEA analysis complete.")
    print(f"Maximum displacement: {results['max_displacement']:.6f} mm")
    print(f"Maximum von Mises stress: {results['max_von_mises']:.2f} MPa")
    
    # Save the visualization
    if 'visualization' in results:
        export_visualization(
            results['visualization'],
            "Results_volco/fea/boundary_examples/simple_mode.html"
        )
    
    print("Simple Mode example complete. Results saved to Results_volco/fea/boundary_examples/simple_mode.html")


def run_simple_mode_tension_example(voxel_matrix, voxel_size):
    """Run FEA analysis using Simple Mode for tension in X direction."""
    print("\n=== Simple Mode Tension Example ===")
    
    # Get model dimensions - use the x-dimension of the voxel matrix
    model_width = voxel_matrix.shape[0] * voxel_size
    
    # Calculate displacement as 1% of model width
    displacement_magnitude = model_width * 0.01
    
    # Define boundary conditions for tension in X direction
    boundary_conditions = {
        'constraints': {
            Surface.MINUS_X: "fix",  # Fix -X surface
            Surface.PLUS_X: [displacement_magnitude, None, None, None, None, None]  # Apply 1% tension on +X
        }
    }
    
    # Run the analysis
    results = analyze_voxel_matrix(
        voxel_matrix=voxel_matrix,
        voxel_size=voxel_size,
        boundary_conditions=boundary_conditions,
        visualization=True,
        result_type='von_mises',
        scale_factor=5.0
    )
    
    print("FEA analysis complete.")
    print(f"Maximum displacement: {results['max_displacement']:.6f} mm")
    print(f"Maximum von Mises stress: {results['max_von_mises']:.2f} MPa")
    
    # Save the visualization
    if 'visualization' in results:
        export_visualization(
            results['visualization'],
            "Results_volco/fea/boundary_examples/simple_mode_tension.html"
        )
    
    print("Simple Mode Tension example complete. Results saved to Results_volco/fea/boundary_examples/simple_mode_tension.html")


def run_expert_mode_example(voxel_matrix, voxel_size):
    """Run FEA analysis using Expert Mode with custom function."""
    print("\n=== Expert Mode Example ===")
    
    # Define custom boundary conditions using expert mode with list comprehensions
    def custom_constraint_function(nodes, elements):
        # Get model dimensions directly
        x_coords = nodes[:, 0]
        y_coords = nodes[:, 1]
        z_coords = nodes[:, 2]
        
        # Create constraints dictionary using list comprehension
        # Calculate model height
        model_height = np.max(z_coords) - np.min(z_coords)
        
        # Calculate displacement as 1% of model height (to match original example)
        displacement_magnitude = model_height * 0.01
        
        # Get min and max coordinates
        x_min, y_min, z_min = np.min(nodes, axis=0)
        x_max, y_max, z_max = np.max(nodes, axis=0)
        
        # Create constraints dictionary with bottom nodes using list comprehension
        constraints = {
            # Fix nodes in bottom 1% using list comprehension
            i: [0, 0, 0, 0, 0, 0] for i in range(len(nodes))
            if nodes[i, 2] < np.min(z_coords) + 0.01 * model_height
        }
        
        # ALTERNATIVE APPROACH: This demonstrates how to use select_nodes_in_box instead of list comprehension
        # Apply displacement to nodes in the top 1% of the model
        top_nodes = select_nodes_in_box(
            nodes,
            [x_min, y_min, z_max - 0.01 * model_height],  # Min coordinates for top box
            [x_max, y_max, z_max]                         # Max coordinates for top box
        )
        
        # Apply displacement to top nodes
        for node_idx in top_nodes:
            constraints[node_idx] = [None, None, -displacement_magnitude, None, None, None]
        
        return constraints
    
    # Define boundary conditions
    boundary_conditions = {
        'constraints': {
            "custom": custom_constraint_function
        }
    }
    
    # Run the analysis
    results = analyze_voxel_matrix(
        voxel_matrix=voxel_matrix,
        voxel_size=voxel_size,
        boundary_conditions=boundary_conditions,
        visualization=True,
        result_type='von_mises',
        scale_factor=10.0
    )
    
    print("FEA analysis complete.")
    print(f"Maximum displacement: {results['max_displacement']:.6f} mm")
    print(f"Maximum von Mises stress: {results['max_von_mises']:.2f} MPa")
    
    # Save the visualization
    if 'visualization' in results:
        export_visualization(
            results['visualization'],
            "Results_volco/fea/boundary_examples/expert_mode.html"
        )
    
    print("Expert Mode example complete. Results saved to Results_volco/fea/boundary_examples/expert_mode.html")


def run_expert_mode_with_utilities_example(voxel_matrix, voxel_size):
    """Run FEA analysis using Expert Mode with functional utilities."""
    print("\n=== Expert Mode with Utilities Example ===")
    
    # Define custom boundary conditions using functional utilities
    def custom_constraint_function(nodes, elements):
        # Get model dimensions
        x_min, y_min, z_min = np.min(nodes, axis=0)
        x_max, y_max, z_max = np.max(nodes, axis=0)

        # Select nodes on two diagonal planes:
        # 1. One at 25% to one side of the main diagonal (x ≈ y + 0.25*(x_max - x_min))
        # 2. Another at 25% to the other side of the main diagonal (x ≈ y - 0.25*(x_max - x_min))
        diagonal_offset = 0.25 * (x_max - x_min)
        tolerance = 0.05 * (x_max - x_min)  # Tolerance for selecting nodes near the diagonal
        
        diagonal1_nodes = select_nodes_by_predicate(
            nodes,
            lambda node: abs((node[0] - node[1]) - diagonal_offset) < tolerance
        )
        
        diagonal2_nodes = select_nodes_by_predicate(
            nodes,
            lambda node: abs((node[0] - node[1]) + diagonal_offset) < tolerance
        )
        
        # Create constraints dictionary
        constraints = {}
        
        # Calculate displacement as 1% of model height
        displacement_magnitude = (z_max - z_min) * 0.01
        
        # Apply displacement to first diagonal nodes that are in the top half (move up)
        for node_idx in diagonal1_nodes:
            if nodes[node_idx, 2] > (z_min + z_max) / 2:
                constraints[node_idx] = [displacement_magnitude,
                                         displacement_magnitude, displacement_magnitude, 0, 0, 0]
        
        # Apply displacement to second diagonal nodes that are in the top half (move down)
        for node_idx in diagonal2_nodes:
            if nodes[node_idx, 2] > (z_min + z_max) / 2:
                constraints[node_idx] = [-displacement_magnitude, -displacement_magnitude, -displacement_magnitude, 0, 0, 0]
        
        return constraints
    
    # Define boundary conditions
    boundary_conditions = {
        'constraints': {
            "custom": custom_constraint_function
        }
    }
    
    # Run the analysis
    results = analyze_voxel_matrix(
        voxel_matrix=voxel_matrix,
        voxel_size=voxel_size,
        boundary_conditions=boundary_conditions,
        visualization=True,
        result_type='displacement',
        scale_factor=10.0
    )
    
    print("FEA analysis complete.")
    print(f"Maximum displacement: {results['max_displacement']:.6f} mm")
    print(f"Maximum von Mises stress: {results['max_von_mises']:.2f} MPa")
    
    # Save the visualization
    if 'visualization' in results:
        export_visualization(
            results['visualization'],
            "Results_volco/fea/boundary_examples/expert_mode_utilities.html"
        )
    
    print("Expert Mode with Utilities example complete. Results saved to Results_volco/fea/boundary_examples/expert_mode_utilities.html")


def main():
    """Run all boundary condition examples."""
    print("Running Enhanced Boundary Condition System Examples")
    
    # Run VOLCO simulation once
    print("Running VOLCO simulation...")
    output = run_simulation(
        gcode_path='examples/gcode_example.gcode',
        printer_config_path='examples/printer_settings.json',
        sim_config_path='examples/simulation_settings.json'
    )
    
    # Get the cropped voxel matrix from the simulation output
    voxel_matrix = output.cropped_voxel_space
    voxel_size = output._simulation.voxel_size
    
    print(f"Simulation complete. Cropped voxel matrix shape: {voxel_matrix.shape}")
    print(f"Voxel size: {voxel_size} mm")
    
    # Run Simple Mode example
    run_simple_mode_example(voxel_matrix, voxel_size)
    
    # Run Simple Mode Tension example
    run_simple_mode_tension_example(voxel_matrix, voxel_size)
    
    # Run Expert Mode example
    run_expert_mode_example(voxel_matrix, voxel_size)
    
    # Run Expert Mode with Utilities example
    run_expert_mode_with_utilities_example(voxel_matrix, voxel_size)
    
    print("\nAll examples complete. Results saved to Results_volco/fea/boundary_examples/")
    print("Files generated:")
    print("  - Results_volco/fea/boundary_examples/simple_mode.html")
    print("  - Results_volco/fea/boundary_examples/simple_mode_tension.html")
    print("  - Results_volco/fea/boundary_examples/expert_mode.html")
    print("  - Results_volco/fea/boundary_examples/expert_mode_utilities.html")


if __name__ == "__main__":
    main()