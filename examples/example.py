"""
Updated Example for Volco Simulation

This script demonstrates how to use the updated Volco simulation with the new
SimulationOutput class and reorganized code structure.
"""

import json
import os
import sys
import numpy as np

# Add the parent directory to the path so we can import volco
sys.path.append(os.path.abspath('..'))

from volco import run_simulation
from app.reporter.mesh import generate_mesh_from_voxels, export_mesh_to_stl
from app.reporter.visualization import color_mesh, visualize_with_trimesh, visualize_with_plotly


# Define the custom simulation config at module level for use in multiple functions
custom_sim_config = {
    "voxel_size": 0.05,
    "step_size": 0.1,
    "x_offset": 2.0,
    "y_offset": 2.0,
    "z_offset": 0.5,
    "simulation_name": "Test_volco",
    "results_folder": "Results_volco",
    "radius_increment": 0.1,
    "solver_tolerance": 0.0001,
    "x_crop": ["all", "all"],
    "y_crop": ["all", "all"],
    "z_crop": [0.0, "all"],
    "consider_acceleration": True
}


def example_1_using_file_paths():
    """Example 1: Using File Paths"""
    print("\n=== Example 1: Using File Paths ===")
    
    # Define file paths
    gcode_path = 'examples/gcode_example.gcode'
    printer_config_path = 'examples/printer_settings.json'
    sim_config_path = 'examples/simulation_settings.json'
    
    # Run the simulation using file paths
    voxel_space, output = run_simulation(
        gcode_path=gcode_path,
        printer_config_path=printer_config_path,
        sim_config_path=sim_config_path
    )
    
    print(f"Voxel space dimensions: {voxel_space.dimensions}")
    return voxel_space, output


def example_2_using_python_variables():
    """Example 2: Using Python Variables"""
    print("\n=== Example 2: Using Python Variables ===")
    
    # Load G-code content from file
    with open('examples/gcode_example.gcode', 'r') as f:
        gcode_content = f.read()
    
    # Add new lines to the G-code
    gcode_content += '\nG0 F7200 X12.0 Y10.0 Z1.2'
    gcode_content += '\nG1 F1000 X14.0 E0.133041'
    
    # Load printer configuration from file
    with open('examples/printer_settings.json', 'r') as f:
        printer_config = json.load(f)
    
    # Load simulation configuration from file
    with open('examples/simulation_settings.json', 'r') as f:
        sim_config = json.load(f)
    
    # Run the simulation using Python variables
    voxel_space, output = run_simulation(
        gcode=gcode_content,
        printer_config=printer_config,
        sim_config=sim_config
    )
    
    print(f"Voxel space dimensions: {voxel_space.dimensions}")
    return voxel_space, output


def example_3_custom_configurations():
    """Example 3: Creating Configuration Dictionaries Directly"""
    print("\n=== Example 3: Creating Configuration Dictionaries Directly ===")
    
    # Create printer configuration dictionary
    custom_printer_config = {
        "nozzle_jerk_speed": 1.0,
        "extruder_jerk_speed": 5.0,
        "nozzle_acceleration": 120.0,
        "extruder_acceleration": 1200.0,
        "feedstock_filament_diameter": 1.75,
        "nozzle_diameter": 0.4
    }
    
    # Use the global custom_sim_config
    global custom_sim_config
    
    gcode_content = """M83
G0 X10 Y10 Z0.3
G1 F7200 X10.0 Y10.0 Z0.3
G1 F1000 X14.0 E0.133041
G0 F7200 X12.0 Y8.0 Z0.5
G1 F1000 Y12.0 E0.133041
G1 F7200 X10.0 Y10.0 Z0.7
G1 F1000 X12.0 E0.067
G1 F1000 X14 Y8.0 E0.1
G1 Z4 E0.2"""
    
    # Run the simulation with custom configurations
    voxel_space, output = run_simulation(
        gcode=gcode_content,
        printer_config=custom_printer_config,
        sim_config=custom_sim_config
    )
    
    print(f"Voxel space dimensions: {voxel_space.dimensions}")
    return voxel_space, output


def working_with_simulation_output(output, sim_config):
    """Working with the SimulationOutput object"""
    print("\n=== Working with SimulationOutput ===")
    
    # Generate the mesh using the mesh module functionality
    print("Generating mesh...")
    mesh = output.generate_mesh()
    
    # Export to STL if needed
    print("Exporting mesh to STL...")
    stl_path = output.export_mesh_to_stl(mesh)
    print(f"STL exported to: {stl_path}")
    
    print("\nDirect access to mesh and visualization modules:")
    
    # Example: Generate a mesh directly from voxel space
    print("- Generating mesh directly from voxel space")
    custom_mesh = generate_mesh_from_voxels(output.cropped_voxel_space, voxel_size=sim_config["voxel_size"])
    
    # Example: Apply custom coloring
    print("- Applying custom coloring")
    colored_mesh = color_mesh(custom_mesh, color_scheme='viridis')
    
    # Example: Create a custom visualization (not shown, but code is ready)
    print("- Creating custom visualization (not displayed)")
    # To display, you would need a GUI environment or Jupyter notebook
    
    return mesh


if __name__ == "__main__":
    # Run the examples
    voxel_space1, output1 = example_1_using_file_paths()
    voxel_space2, output2 = example_2_using_python_variables()
    voxel_space3, output3 = example_3_custom_configurations()
    
    # Work with the output from the last example
    mesh = working_with_simulation_output(output3, custom_sim_config)
    
    print("\nExamples completed successfully!")
    print("Note: To visualize the results, run this code in a Jupyter notebook or GUI environment.")