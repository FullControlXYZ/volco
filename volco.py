import logging
import json

from app.configs.printer import Printer
from app.configs.simulation import Simulation
from app.configs.arguments import Arguments
from app.geometry.voxel_space import VoxelSpace
from app.instructions.gcode import Gcode
from app.reporter.report import SimulationOutput


logging.basicConfig(format="%(levelname)s %(asctime)s %(message)s", level=logging.INFO)


def run_simulation(
    gcode=None,
    printer_config=None,
    sim_config=None,
    gcode_path=None,
    printer_config_path=None,
    sim_config_path=None,
):
    """
    Run the Volco simulation from a Python environment (e.g., Jupyter notebook).
    
    Parameters:
    -----------
    gcode : str
        G-code content as a string
    printer_config : dict
        Printer configuration as a dictionary
    sim_config : dict
        Simulation configuration as a dictionary
    gcode_path : str
        Path to G-code file (alternative to gcode parameter)
    printer_config_path : str
        Path to printer configuration file (alternative to printer_config parameter)
    sim_config_path : str
        Path to simulation configuration file (alternative to sim_config parameter)
        
    Returns:
    --------
    tuple
        (voxel_space, output) - The VoxelSpace and SimulationOutput objects for further use
        
    Notes:
    ------
    After getting the simulation output object, you can:
    1. Generate a mesh: mesh = output.generate_mesh()
    2. Export to STL: output.export_mesh_to_stl(mesh)
    3. Visualize with trimesh: scene = output.visualize_mesh(mesh)
    4. Visualize with Plotly: fig = output.visualize_mesh(mesh, visualizer='plotly')
    """
    # Initialize printer configuration
    if printer_config:
        printer = Printer(config_dict=printer_config)
    elif printer_config_path:
        printer = Printer(config_path=printer_config_path)
    else:
        raise ValueError("Either printer_config or printer_config_path must be provided")
    
    # Initialize G-code instruction
    default_nozzle_speed = 50.0  # Default nozzle speed in mm/s
    if gcode:
        instruction = Gcode(gcode_content=gcode, default_nozzle_speed=default_nozzle_speed, printer=printer)
    elif gcode_path:
        instruction = Gcode(gcode_path=gcode_path, default_nozzle_speed=default_nozzle_speed, printer=printer)
    else:
        raise ValueError("Either gcode or gcode_path must be provided")
    
    instruction.read()
    print(f"Number of printed filaments: {instruction.number_printed_filaments}")
    
    # Initialize simulation configuration
    if sim_config:
        simulation = Simulation(config_dict=sim_config)
    elif sim_config_path:
        simulation = Simulation(config_path=sim_config_path)
    else:
        raise ValueError("Either sim_config or sim_config_path must be provided")
    
    # Initialize voxel space
    voxel_space = VoxelSpace(
        instruction=instruction,
        simulation_config=simulation,
        printer=printer,
    )
    
    voxel_space.initialize_space()
    voxel_space.print()
    
    # Process simulation output
    output = SimulationOutput(voxel_space=voxel_space, simulation=simulation)
    
    # Crop the voxel space
    output.crop_voxel_space()
    
    # For CLI usage, generate and export STL automatically
    if __name__ == "__main__":
        mesh = output.generate_mesh()
        output.export_mesh_to_stl(mesh)
    
    return voxel_space, output


if __name__ == "__main__":
    options = Arguments.get_options()

    run_simulation(
        gcode_path=options.gcode,
        printer_config_path=options.printer,
        sim_config_path=options.sim
    )
