"""
Example script demonstrating how to save FEA results to a file.

This script shows how to:
1. Run a VOLCO simulation
2. Perform FEA analysis on the resulting voxel matrix
3. Save the FEA results to a file in different formats
"""

import os
import numpy as np
from volco import run_simulation
from app.postprocessing.fea import analyze_voxel_matrix
from app.postprocessing.fea.viz import export_visualization


def main():
    """Run the FEA example with result saving."""
    print("Running VOLCO simulation...")
    
    # Run VOLCO simulation
    output = run_simulation(
        gcode_path='examples/gcode_example.gcode',
        printer_config_path='examples/printer_settings.json',
        sim_config_path='examples/simulation_settings.json'
    )
    
    # Get the voxel matrix from the simulation output
    voxel_matrix = output.voxel_space.space
    voxel_size = output._simulation.voxel_size
    
    print(f"Simulation complete. Voxel matrix shape: {voxel_matrix.shape}")
    print(f"Voxel size: {voxel_size} mm")
    
    # Create output directory if it doesn't exist
    os.makedirs("Results_volco/fea", exist_ok=True)
    
    # Define material properties for PLA
    material_properties = {
        'young_modulus': 2000.0,  # MPa (typical for PLA)
        'poisson_ratio': 0.3      # Typical for PLA
    }
    
    # Define custom boundary conditions
    boundary_conditions = {
        'displacement_percentage': 1.0  # 1% displacement
    }
    
    # Save results in pickle format (default)
    print("\nRunning FEA analysis and saving results in pickle format...")
    results_pickle = analyze_voxel_matrix(
        voxel_matrix=voxel_matrix,
        voxel_size=voxel_size,
        material_properties=material_properties,
        boundary_conditions=boundary_conditions,
        visualization=True,
        result_type='von_mises',
        scale_factor=10.0,
        save_results=True,
        save_path='Results_volco/fea/results_pickle',
        save_format='pickle',
        include_visualization=False
    )
    
    print(f"Results saved to: {results_pickle['saved_file']}")
    print(f"Maximum displacement: {results_pickle['max_displacement']:.6f} mm")
    print(f"Maximum von Mises stress: {results_pickle['max_von_mises']:.2f} MPa")
    
    # Save results in JSON format
    print("\nRunning FEA analysis and saving results in JSON format...")
    results_json = analyze_voxel_matrix(
        voxel_matrix=voxel_matrix,
        voxel_size=voxel_size,
        material_properties=material_properties,
        boundary_conditions=boundary_conditions,
        visualization=True,
        result_type='von_mises',
        scale_factor=10.0,
        save_results=True,
        save_path='Results_volco/fea/results_json',
        save_format='json',
        include_visualization=False
    )
    
    print(f"Results saved to: {results_json['saved_file']}")
    
    # Try to save results in HDF5 format (if h5py is available)
    try:
        import h5py
        print("\nRunning FEA analysis and saving results in HDF5 format...")
        results_hdf5 = analyze_voxel_matrix(
            voxel_matrix=voxel_matrix,
            voxel_size=voxel_size,
            material_properties=material_properties,
            boundary_conditions=boundary_conditions,
            visualization=True,
            result_type='von_mises',
            scale_factor=10.0,
            save_results=True,
            save_path='Results_volco/fea/results_hdf5',
            save_format='hdf5',
            include_visualization=False
        )
        
        print(f"Results saved to: {results_hdf5['saved_file']}")
    except ImportError:
        print("\nSkipping HDF5 format (h5py not available)")
    
    # Save visualization
    if 'visualization' in results_pickle:
        print("\nSaving FEA visualization...")
        export_visualization(
            results_pickle['visualization'],
            "Results_volco/fea/fea_results_saved.html"
        )
    
    print("\nExample complete. Results saved to Results_volco/fea/")
    print("Files generated:")
    print("  - Results_volco/fea/results_pickle.pkl")
    print("  - Results_volco/fea/results_json.json")
    print("  - Results_volco/fea/fea_results_saved.html")
    if 'h5py' in globals():
        print("  - Results_volco/fea/results_hdf5.h5")


if __name__ == "__main__":
    main()