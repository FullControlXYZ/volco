import pytest
from app.geometry.geometry_math import GeometryMath
from app.instructions.gcode import Gcode

from app.instructions.instruction import Instruction
from app.geometry.voxel_space import VoxelSpace


class FakeInstruction(Instruction):
    def __init__(self):
        self._movements = list()
        self._coordinate_limits = {"x": [2.0, 10.0], "y": [4.0, 16.0], "z": [0.0, 5.0]}
        self._number_printed_filaments = 0
        self._filaments_coordinates = list()
        self._default_nozzle_speed = 10.0

    def read(self):
        pass

    @property
    def movements(self):
        return self._movements

    @property
    def coordinate_limits(self):
        return self._coordinate_limits

    @property
    def number_printed_filaments(self):
        return self._number_printed_filaments

    @property
    def filaments_coordinates(self):
        return self._filaments_coordinates

    @property
    def default_nozzle_speed(self):
        return self._default_nozzle_speed


class FakeSimulation:
    def __init__(self):
        self.x_offset = 2.0
        self.y_offset = 3.0
        self.z_offset = 1.0
        self.voxel_size = 0.1
        self.step_size = 0.2
        self.radius_increment = 0.1
        self.solver_tolerance = 0.001


class FakePrinter:
    def __init__(self):
        self.nozzle_jerk_speed = 40.0
        self.extruder_jerk_speed = 5.0
        self.nozzle_acceleration = 1200.0
        self.extruder_acceleration = 1200.0
        self.bulk_filament_diameter = 1.75
        self.nozzle_diameter = 0.4


class TestVoxelSpace:
    def test_should_calculate_dimensions_and_translations(self):
        test_instruction = FakeInstruction()
        simulation_config = FakeSimulation()
        printer = FakePrinter()

        voxel_space = VoxelSpace(
            instruction=test_instruction,
            simulation_config=simulation_config,
            printer=printer,
        )

        assert voxel_space.dimensions == {"x": 120, "y": 180, "z": 60}

        assert voxel_space.filament_translations == {"x": 0.0, "y": -1.0}

    def test_should_initialize_voxel_space(self):
        test_instruction = FakeInstruction()
        simulation_config = FakeSimulation()
        printer = FakePrinter()

        voxel_space = VoxelSpace(
            instruction=test_instruction,
            simulation_config=simulation_config,
            printer=printer,
        )

        voxel_space.initialize_space()

        assert voxel_space.space.shape == (120, 180, 60)

    def test_should_print(self):
        test_instruction = FakeInstruction()
        simulation_config = FakeSimulation()
        printer = FakePrinter()

        voxel_space = VoxelSpace(
            instruction=test_instruction,
            simulation_config=simulation_config,
            printer=printer,
        )

        voxel_space.initialize_space()
        voxel_space.print()


class FakeSimulationPrint:
    def __init__(self):
        self.x_offset = 2.0
        self.y_offset = 2.0
        self.z_offset = 0.5
        self.voxel_size = 0.05
        self.step_size = 0.2
        self.radius_increment = 0.1
        self.solver_tolerance = 0.001


class TestPrint:
    def test_should_print(self):
        test_instruction = Gcode(
            gcode_path="tests/fixtures/gcode_example.gcode", default_nozzle_speed=40.0
        )
        test_instruction.read()

        simulation_config = FakeSimulationPrint()
        printer = FakePrinter()

        voxel_space = VoxelSpace(
            instruction=test_instruction,
            simulation_config=simulation_config,
            printer=printer,
        )

        voxel_space.initialize_space()
        voxel_space.print()

        total_deposited_volume = GeometryMath.calculate_filled_volume(
            voxel_space.space, simulation_config.voxel_size
        )
        assert (total_deposited_volume - 0.9595) < 1e-6
