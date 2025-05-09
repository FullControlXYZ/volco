# this text needs updating and additional tests created since the report module was split into multiple modules

# import pytest
# import numpy as np
# from app.geometry.voxel_space import VoxelSpace
# from app.instructions.gcode import Gcode

# from app.reporter.report import SimulationOutput


# class FakePrinter:
#     def __init__(self):
#         self.nozzle_jerk_speed = 40.0
#         self.extruder_jerk_speed = 5.0
#         self.nozzle_acceleration = 1200.0
#         self.extruder_acceleration = 1200.0
#         self.feedstock_filament_diameter = 1.75
#         self.nozzle_diameter = 0.4


# class FakeSimulationPrint:
#     def __init__(self):
#         self.x_offset = 2.0
#         self.y_offset = 2.0
#         self.z_offset = 0.5
#         self.voxel_size = 0.05
#         self.step_size = 0.2
#         self.radius_increment = 0.1
#         self.solver_tolerance = 0.001
#         self.x_crop = [11.0, 13.0]
#         self.y_crop = ["all", "all"]
#         self.z_crop = [0.0, "all"]


# class TestCropVoxelSpace:
#     def test_should_crop_voxel_space(self):
#         test_instruction = Gcode(
#             gcode_path="tests/fixtures/gcode_example.gcode", default_nozzle_speed=40.0
#         )
#         test_instruction.read()

#         simulation_config = FakeSimulationPrint()
#         printer = FakePrinter()

#         voxel_space = VoxelSpace(
#             instruction=test_instruction,
#             simulation_config=simulation_config,
#             printer=printer,
#         )

#         voxel_space.initialize_space()

#         report = Report(voxel_space=voxel_space, simulation=simulation_config)

#         report.crop_voxel_space()

#         assert report.cropped_voxel_space.shape == (41, 160, 23)

