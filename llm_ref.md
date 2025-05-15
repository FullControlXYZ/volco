# VOLCO - LLM Reference Document

This document provides a comprehensive overview of the Volco repository structure and functionality to minimize context needed for LLM code editing. It serves as a reference for large language models to understand the codebase without requiring all modules to be uploaded as context.

> **Note**: This reference document is based on examination of key files including volco.py, app/geometry/voxel_space.py, app/instructions/gcode.py, app/reporter/report.py, app/configs/simulation.py, app/configs/printer.py, app/geometry/sphere.py, and app/physics/acceleration modules. It may not cover aspects contained in files that were not examined. As the repository evolves, this document should be updated accordingly.

## 1. Repository Overview
VOLCO (VOLume COnserving) is a 3D printing simulation model that predicts the final shape of printed materials by simulating the printing process in a voxelized space. It takes G-code as input and produces a voxel 3D-matrix or STL file as output, enabling accurate prediction of material distribution during extrusion-based 3D printing.

## 2. Directory Structure
```
volco/
├── app/                    # Core application modules
│   ├── configs/            # Configuration handling
│   ├── geometry/           # Spatial calculations and voxel operations
│   ├── instructions/       # G-code parsing and instruction handling
│   ├── physics/            # Physics simulation (acceleration, volume)
│   │   └── acceleration/   # Speed profiles and acceleration models
│   ├── reporter/           # Output generation and visualization
│   └── solvers/            # Mathematical solvers
├── examples/               # Example files and usage demonstrations
├── tests/                  # Test suite mirroring app structure
└── volco.py                # Main entry point
```

## 3. Core Components

### Main Entry Point (volco.py)
- Provides `run_simulation()` function that orchestrates the entire simulation process
- Accepts G-code, printer configuration, and simulation configuration as inputs
- Returns a `SimulationOutput` object containing the voxel space and mesh data

### Configuration (app/configs/)
- `Printer` class: Manages printer hardware settings (nozzle diameter, acceleration, etc.)
- `Simulation` class: Manages simulation parameters (voxel size, step size, etc.)
- `Arguments` class: Handles command-line argument parsing

### G-code Processing (app/instructions/)
- `Instruction` class: Abstract base class for all instruction types
- `Gcode` class: Parses G-code files and extracts movement coordinates and extrusion information
- Key methods: `read()`, `_define_movement()`, `_max_min_extru_coordinates()`

### Geometry (app/geometry/)
- `VoxelSpace` class: Manages the 3D voxel grid where the simulation occurs
- `Sphere` class: Models material deposition as spheres in the voxel space
- `GeometryMath` class: Provides utility functions for geometric calculations
- Key methods: `VoxelSpace.initialize_space()`, `VoxelSpace.print()`, `Sphere.deposit_sphere()`

### Physics (app/physics/)
- `Volume` class: Handles volume calculations for material deposition
- `AccelerationManager` class: Calculates speed profiles for nozzle and extruder
- `NozzleSpeed` class: Models the nozzle movement speed
- `ExtruderSpeed` class: Models the extruder feed rate
- `SpeedProfile` class: Abstract base class for speed profiles
- `TrapezoidalSpeedProfile` class: Implements a trapezoidal speed profile
- `FlatSpeedProfile` class: Implements a constant speed profile

### Reporting (app/reporter/)
- `SimulationOutput` class: Processes and visualizes simulation results
- `mesh` module: Generates and exports 3D meshes
- `visualization` module: Creates visualizations of simulation results
- Key methods: `crop_voxel_space()`, `generate_mesh()`, `export_mesh_to_stl()`, `visualize_mesh()`

### Solvers (app/solvers/)
- `BisectionMethod` class: Implements the bisection method for finding roots
- Key method: `execute()`: Executes the bisection method algorithm

## 4. Key Data Flows

### G-code to Voxel Space
1. `Gcode.read()` → Parses G-code into movement coordinates
2. `VoxelSpace.initialize_space()` → Creates empty voxel grid
3. `VoxelSpace.print()` → For each filament segment:
   - Calculate direction vector and filament length
   - Determine volume distribution with `Volume.get_volumes_for_filament()`
   - Deposit material with `_deposit_filament()`
4. `Sphere.deposit_sphere()` → Models material deposition at each step

### Voxel Space to Mesh Output
1. `SimulationOutput.crop_voxel_space()` → Crops voxel space to region of interest
2. `SimulationOutput.generate_mesh()` → Converts voxels to 3D mesh using marching cubes
3. `SimulationOutput.export_mesh_to_stl()` → Exports mesh to STL file
4. `SimulationOutput.visualize_mesh()` → Creates visualization with Trimesh or Plotly

### Acceleration-Based Volume Distribution
1. `AccelerationManager.calculate_speed_profiles()` → Calculates nozzle and extruder speed profiles
2. `NozzleSpeed.calculate_displacements()` → Computes nozzle position over time
3. `ExtruderSpeed.calculate_displacements()` → Computes extruder position over time
4. `Volume.get_volumes_for_filament()` → Distributes volume based on speed profiles
5. Material deposition adjusted according to calculated volumes

## 5. Configuration Options

See README.md for default values

### Simulation Configuration
- `voxel_size`: Size of each voxel in mm (smaller = higher accuracy but more memory usage)
- `step_size`: Distance between simulation steps in mm (smaller = higher accuracy but slower)
- `x_offset`, `y_offset`, `z_offset`: Offsets for the voxel space boundaries (default: 5 * nozzle_diameter)
- `sphere_z_offset`: Distance to offset the sphere center below the nozzle (default: 0.5 * nozzle_diameter)
- `x_crop`, `y_crop`, `z_crop`: Cropping boundaries for the final output (default: ["all", "all"] for each)
- `radius_increment`: Increment for sphere radius in bisection method
- `solver_tolerance`: Tolerance for volume conservation in bisection method (default: 0.0001)
- `consider_acceleration`: Whether to consider acceleration in volume distribution (default: false)
- `stl_ascii`: Whether to export STL in ASCII format (default: false)

### Printer Configuration
- `nozzle_diameter`: Diameter of the printer nozzle in mm
- `feedstock_filament_diameter`: Diameter of the input filament in mm
- `nozzle_jerk_speed`: Maximum instantaneous speed change for nozzle in mm/s
- `extruder_jerk_speed`: Maximum instantaneous speed change for extruder in mm/s
- `nozzle_acceleration`: Acceleration of the nozzle in mm/s²
- `extruder_acceleration`: Acceleration of the extruder in mm/s²

## 6. Usage Patterns

### Command Line Usage
```bash
python volco.py --gcode=examples/gcode_example.gcode --sim=examples/simulation_settings.json --printer=examples/printer_settings.json
```

### Python API Usage
```python
from volco import run_simulation

# Using file paths
output = run_simulation(
    gcode_path='examples/gcode_example.gcode',
    printer_config_path='examples/printer_settings.json',
    sim_config_path='examples/simulation_settings.json'
)

# Using Python variables
output = run_simulation(
    gcode=gcode_content,  # G-code as a string
    printer_config=printer_config_dict,  # Printer config as a dictionary
    sim_config=sim_config_dict  # Simulation config as a dictionary
)

# Working with output
output.export_mesh_to_stl()
fig = output.visualize_mesh(visualizer='plotly')
```

## 7. Key Algorithms

- **G-code Parsing**: Converts G-code commands to 3D coordinates and extrusion values, handling different positioning modes and unit conversions
- **Voxel Space Initialization**: Creates a 3D grid based on printing dimensions and voxel size
- **Material Deposition**: Uses sphere-based model to deposit material in voxel space, with radius determined by volume conservation
- **Volume Conservation**: Bisection method to determine sphere radius for accurate volume deposition
- **Acceleration Profiles**: Trapezoidal speed profiles for realistic nozzle and extruder movement
- **Mesh Generation**: Marching cubes algorithm to convert voxel space to 3D mesh

## 8. Module Dependencies

- `volco.py` → `app/configs`, `app/instructions`, `app/geometry`, `app/reporter`
- `app/geometry/voxel_space.py` → `app/geometry/sphere.py`, `app/physics/volume.py`
- `app/physics/volume.py` → `app/physics/acceleration/acceleration_profiles.py`
- `app/reporter/report.py` → `app/reporter/mesh.py`, `app/reporter/visualization.py`
- `app/instructions/gcode.py` → `app/instructions/instruction.py`
- `app/geometry/sphere.py` → `app/solvers/bisection_method.py`

## 9. Common Modification Patterns

### Adding New Speed Profiles
1. Create a new class in `app/physics/acceleration/` that extends `SpeedProfile`
2. Implement the `calculate_displacements()` method
3. Update `Volume.get_volumes_for_filament()` to use the new profile

### Supporting Additional G-code Commands
1. Modify `app/instructions/gcode.py` to recognize new commands
2. Update the `read()` method to handle the new commands
3. Implement any necessary processing in `_define_movement()`

### Adding New Visualization Options
1. Add new visualization function in `app/reporter/visualization.py`
2. Update `SimulationOutput.visualize_mesh()` to support the new option
3. Add any necessary parameters to the visualization method

### Modifying Material Deposition Model
1. Update `app/geometry/sphere.py` or create a new deposition model class
2. Modify `VoxelSpace._deposit_filament()` to use the new model
3. Update volume calculation methods if necessary
4. Adjust the `sphere_z_offset` centroid of spheres below nozzle

## 10. Performance Considerations

- **Voxel Size**: Smaller voxel sizes increase accuracy but exponentially increase memory usage and computation time
- **Step Size**: Smaller step sizes increase simulation accuracy but increase computation time
- **Acceleration Consideration**: Enabling acceleration profiles increases accuracy but adds computational overhead
- **Mesh Generation**: Marching cubes algorithm can be memory-intensive for large voxel spaces
- **Cropping**: Proper cropping of the voxel space before mesh generation can significantly improve performance
- **Bisection Method Tolerance**: Lower tolerance values increase volume accuracy but require more iterations
- **Sphere Radius Increment**: Smaller increments increase accuracy but require more iterations in the bisection method