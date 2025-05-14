# Volco LLM Reference Guide

This document provides a comprehensive overview of the Volco repository structure and functionality to minimize context needed for LLM code editing. It serves as a reference for large language models to understand the codebase without requiring all modules to be uploaded as context.

> **Note**: This reference document is based on examination of 18 key files in the repository that are listed at the end of this document. It may not cover aspects contained in files that were not examined. As the repository evolves, this document should be updated accordingly.

## Repository Overview

Volco (VOLume COnserving) is a 3D printing simulation model that predicts the final shape of 3D printed materials. It simulates the printing process in a voxelized space, taking G-code as input and producing a voxel 3D-matrix or STL file as output.

## Directory Structure

```
volco/
├── app/                      # Main application code
│   ├── configs/              # Configuration handling
│   ├── geometry/             # Geometric operations and voxel space
│   ├── instructions/         # G-code parsing and instruction handling
│   ├── physics/              # Physics simulation components
│   │   └── acceleration/     # Speed and acceleration profiles
│   ├── reporter/             # Output generation and visualization
│   └── solvers/              # Numerical solvers
├── examples/                 # Example files and usage demonstrations
├── tests/                    # Test suite
└── volco.py                  # Main entry point
```

## Core Components

### Main Entry Point (volco.py)

The main entry point that orchestrates the simulation process. It provides the `run_simulation()` function which:
- Takes G-code, printer configuration, and simulation configuration as inputs
- Initializes the simulation components
- Processes the G-code instructions
- Generates the voxel space
- Returns a `SimulationOutput` object with the results

**Key Functions:**
- `run_simulation()`: Orchestrates the entire simulation process

### Configuration (app/configs/printer.py, app/configs/simulation.py, examples/printer_settings.json, examples/simulation_settings.json)

Handles loading and managing configuration settings for the simulation.

**Key Classes:**
- `Printer`: Manages printer hardware settings (nozzle diameter, acceleration, etc.)
- `Simulation`: Manages simulation parameters (voxel size, step size, etc.)
- `Arguments`: Handles command-line argument parsing

### G-code Processing (app/instructions/gcode.py, app/instructions/instruction.py)

Parses and processes G-code instructions for the simulation.

**Key Classes:**
- `Instruction`: Abstract base class for all instruction types
- `Gcode`: Parses G-code files and extracts movement coordinates and extrusion information

**Key Methods:**
- `Gcode.read()`: Parses G-code content and extracts movement data
- `Gcode._define_movement()`: Calculates movement coordinates based on G-code commands
- `Gcode._max_min_extru_coordinates()`: Determines printing coordinate limits and filament segments

### Geometry (app/geometry/voxel_space.py, app/geometry/sphere.py, app/geometry/geometry_math.py)

Handles geometric operations and the voxel space representation.

**Key Classes:**
- `VoxelSpace`: Manages the 3D voxel grid where the simulation occurs
- `Sphere`: Represents material deposition as spheres in the voxel space
- `GeometryMath`: Provides utility functions for geometric calculations

**Key Methods:**
- `VoxelSpace.initialize_space()`: Creates the voxel grid
- `VoxelSpace.print()`: Simulates the printing process in the voxel space
- `Sphere.deposit_sphere()`: Deposits material as a sphere in the voxel space
- `GeometryMath.calculate_filled_volume()`: Calculates the volume of filled voxels

### Physics (app/physics/volume.py, app/physics/acceleration/acceleration_profiles.py, app/physics/acceleration/speed_profile.py, app/physics/acceleration/nozzle_speed.py, app/physics/acceleration/trapezoidal_speed_profile.py, app/physics/acceleration/volume.py)

Models the physical aspects of the 3D printing process.

**Key Classes:**
- `Volume`: Handles volume calculations for material deposition
- `AccelerationVolume`: Calculates volume distribution considering acceleration
- `AccelerationManager`: Manages speed profiles for nozzle and extruder
- `NozzleSpeed`: Models the nozzle movement speed
- `ExtruderSpeed`: Models the extruder feed rate
- `SpeedProfile`: Abstract base class for speed profiles
- `TrapezoidalSpeedProfile`: Implements a trapezoidal speed profile
- `FlatSpeedProfile`: Implements a constant speed profile

**Key Methods:**
- `Volume.get_volumes_for_filament()`: Calculates volumes for filament segments
- `AccelerationVolume.calculate_volumes_with_acceleration()`: Distributes volume considering acceleration
- `NozzleSpeed.calculate_displacements()`: Calculates nozzle position over time
- `ExtruderSpeed.calculate_displacements()`: Calculates extruder position over time

### Reporting (app/reporter/report.py, app/reporter/mesh.py, app/reporter/visualization.py)

Processes simulation results and generates output.

**Key Classes:**
- `SimulationOutput`: Processes and visualizes simulation results
- `mesh`: Module for generating and exporting 3D meshes
- `visualization`: Module for visualizing simulation results

**Key Methods:**
- `SimulationOutput.crop_voxel_space()`: Crops the voxel space to the region of interest
- `SimulationOutput.generate_mesh()`: Generates a 3D mesh from the voxel data
- `SimulationOutput.export_mesh_to_stl()`: Exports the mesh to an STL file
- `SimulationOutput.visualize_mesh()`: Creates visualizations of the mesh
- `generate_mesh_from_voxels()`: Converts voxel data to a 3D mesh
- `export_mesh_to_stl()`: Exports a mesh to an STL file
- `visualize_with_trimesh()`: Creates a visualization using trimesh
- `visualize_with_plotly()`: Creates a visualization using Plotly

### Solvers (app/solvers/bisection_method.py)

Provides numerical solvers for the simulation.

**Key Classes:**
- `BisectionMethod`: Implements the bisection method for finding roots

**Key Methods:**
- `BisectionMethod.execute()`: Executes the bisection method algorithm

## Key Data Flows

1. **Input Processing**:
   - G-code file → `Gcode` class → Movement coordinates and extrusion data
   - Configuration files → `Printer` and `Simulation` classes → Simulation parameters

2. **Simulation Process**:
   - Movement coordinates → `VoxelSpace.print()` → Voxel space with deposited material
   - For each filament segment:
     - Calculate volumes using `Volume.get_volumes_for_filament()`
     - Deposit material using `Sphere.deposit_sphere()`

3. **Output Generation**:
   - Voxel space → `SimulationOutput.generate_mesh()` → 3D mesh
   - 3D mesh → `SimulationOutput.export_mesh_to_stl()` → STL file
   - 3D mesh → `SimulationOutput.visualize_mesh()` → Visual representation

## Configuration Options

### Printer Settings (`printer_settings.json`)
Parameters: nozzle_jerk_speed, extruder_jerk_speed, nozzle_acceleration, extruder_acceleration, feedstock_filament_diameter, nozzle_diameter

### Simulation Settings (`simulation_settings.json`)
Parameters: voxel_size, step_size, x_offset, y_offset, z_offset, simulation_name, results_folder, radius_increment, solver_tolerance, x_crop, y_crop, z_crop, consider_acceleration

## Usage Patterns

### Command Line Usage
`python volco.py --gcode=examples/gcode_example.gcode --sim=examples/simulation_settings.json --printer=examples/printer_settings.json`

### Python API Usage
- Import: `from volco import run_simulation`
- Run with files: `output = run_simulation(gcode_path=path, printer_config_path=path, sim_config_path=path)`
- Run with variables: `output = run_simulation(gcode=content, printer_config=dict, sim_config=dict)`
- Export: `output.export_mesh_to_stl()`
- Visualize: `output.visualize_mesh(visualizer='plotly', color_scheme='cyan_blue')`

## Key Algorithms

- **G-code Parsing**: Extract movement coordinates, convert extrusion length (E) to volume
- **Voxel Space Simulation**: Deposit material as spheres, adjust radius to match volume
- **Speed Profile Calculation**: Model nozzle and extruder movement with acceleration
- **Mesh Generation**: Convert voxels to mesh, apply colors, export to STL

## Module Dependencies

- `volco.py` → `app/configs/`, `app/instructions/`, `app/geometry/`, `app/reporter/`
- `app/geometry/voxel_space.py` → `app/geometry/sphere.py`, `app/physics/volume.py`
- `app/physics/volume.py` → `app/physics/acceleration/volume.py` (when acceleration is considered)
- `app/reporter/report.py` → `app/reporter/mesh.py`, `app/reporter/visualization.py`
- `app/instructions/gcode.py` → `app/instructions/instruction.py`
- `app/physics/acceleration/` modules have interdependencies for speed profile calculations

## Common Modification Patterns

1. **Adding a new speed profile**:
   - Create a new class that extends `SpeedProfile`
   - Implement the required abstract methods
   - Update `NozzleSpeed` or `ExtruderSpeed` to use the new profile

2. **Modifying material deposition**:
   - Update the `Sphere` class or create a new deposition shape
   - Modify `VoxelSpace._deposit_filament()` to use the new deposition method

3. **Adding visualization options**:
   - Add new visualization methods to `app/reporter/visualization.py`
   - Update `SimulationOutput.visualize_mesh()` to support the new options

4. **Adding G-code commands**:
   - Update `Gcode.read()` to handle new G-code commands
   - Add corresponding logic to process the commands

## Performance Considerations

- Voxel size significantly impacts memory usage and processing time
- The `consider_acceleration` option adds computational complexity
- Large G-code files may require significant memory
- Mesh generation can be memory-intensive for high-resolution simulations

> **Note**: This reference also examined README.md for project overview. Other files in the repository may contain additional functionality not covered here.