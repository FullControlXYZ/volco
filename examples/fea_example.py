"""
Example script demonstrating the use of the FEA module with VOLCO.

This script shows how to:
1. Run a VOLCO simulation
2. Pass the resulting voxel matrix to the FEA module
3. Visualize the stress results
"""

import os
import numpy as np
from volco import run_simulation
from app.postprocessing.fea import analyze_voxel_matrix
from app.postprocessing.fea.viz import export_visualization, visualize_fea


def main():
    """Run the FEA example."""
    print("Running VOLCO simulation...")
    
    # Run VOLCO simulation
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
    
    # Create output directory if it doesn't exist
    os.makedirs("Results_volco/fea", exist_ok=True)
    
    # Note: Original voxel matrix visualization is now handled by the Plotly-based functions
    
    # Run FEA analysis
    print("Running FEA analysis...")
    
    # Define material properties for PLA
    material_properties = {
        'young_modulus': 2000.0,  # MPa (typical for PLA)
        'poisson_ratio': 0.3      # Typical for PLA
    }
    
    # Define custom boundary conditions (optional)
    boundary_conditions = {
        'displacement_percentage': 1.0  # 1% displacement
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
    
    # Create von Mises stress visualization from the results
    print("Creating default von Mises visualization...")
    if 'visualization' in results:
        export_visualization(
            results['visualization'],
            "Results_volco/fea/von_mises.html"
        )
    
    # Create visualizations using the new unified visualization function
    
    # Create von Mises stress visualization with undeformed mesh
    print("Creating von Mises stress visualization with undeformed mesh...")
    von_mises_fig = visualize_fea(
        nodes=results['nodes'],
        elements=results['elements'],
        displacements=results['displacements'],
        von_mises=results['von_mises'],
        result_type='von_mises',
        scale_factor=10.0,  # Exaggerate deformation for visualization
        show_undeformed=True,
        original_opacity=0.3  # Set transparency of original mesh
    )
    
    # Save the von Mises visualization
    export_visualization(
        von_mises_fig,
        "Results_volco/fea/von_mises_with_undeformed.html"
    )
    
    # Create displacement visualization with undeformed mesh
    print("Creating displacement visualization with undeformed mesh...")
    displacement_fig = visualize_fea(
        nodes=results['nodes'],
        elements=results['elements'],
        displacements=results['displacements'],
        von_mises=results['von_mises'],
        result_type='displacement',
        scale_factor=10.0,  # Exaggerate deformation for visualization
        show_undeformed=True,
        original_opacity=0.3  # Set transparency of original mesh
    )
    
    # Save the displacement visualization
    export_visualization(
        displacement_fig,
        "Results_volco/fea/displacement_with_undeformed.html"
    )
    
    # Create displacement visualization without undeformed mesh
    print("Creating displacement visualization without undeformed mesh...")
    displacement_only_fig = visualize_fea(
        nodes=results['nodes'],
        elements=results['elements'],
        displacements=results['displacements'],
        von_mises=results['von_mises'],
        result_type='displacement',
        scale_factor=10.0,  # Exaggerate deformation for visualization
        show_undeformed=False
    )
    
    # Save the displacement-only visualization
    export_visualization(
        displacement_only_fig,
        "Results_volco/fea/displacement_only.html"
    )
    
    print("Example complete. Results saved to Results_volco/fea/")
    print("Files generated:")
    print("  - Results_volco/fea/von_mises.html")
    print("  - Results_volco/fea/von_mises_with_undeformed.html")
    print("  - Results_volco/fea/displacement_with_undeformed.html")
    print("  - Results_volco/fea/displacement_only.html")


if __name__ == "__main__":
    main()