import logging

from app.configs.printer import Printer
from app.configs.simulation import Simulation
from app.configs.arguments import Arguments
from app.geometry.voxel_space import VoxelSpace
from app.instructions.gcode import Gcode
from app.reporter.report import Report


logging.basicConfig(format="%(levelname)s %(asctime)s %(message)s", level=logging.INFO)

if __name__ == "__main__":
    options = Arguments.get_options()

    printer = Printer(config_path=options.printer)
    
    instruction = Gcode(gcode_path=options.gcode, default_nozzle_speed=40.0, printer=printer)
    instruction.read()
    print(instruction.number_printed_filaments)

    simulation = Simulation(config_path=options.sim)

    voxel_space = VoxelSpace(
        instruction=instruction,
        simulation_config=simulation,
        printer=printer,
    )

    voxel_space.initialize_space()
    voxel_space.print()

    report = Report(voxel_space=voxel_space, simulation=simulation)

    report.crop_voxel_space()
    report.generate_stl()
