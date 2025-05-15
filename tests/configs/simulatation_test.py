import pytest

from app.configs.simulation import Simulation


class TestSimulation:
    def test_should_load_simulation_settings(self):
        simulation = Simulation(config_path="tests/fixtures/simulation_settings.json")

        assert simulation.voxel_size == 0.05
        assert simulation.step_size == 0.2
        assert simulation.x_offset == 2.0
        assert simulation.y_offset == 2.0
        assert simulation.z_offset == 0
        assert simulation.sphere_z_offset == 0.2
        assert simulation.simulation_name == "Test_volco"
        assert simulation.results_folder == "Results_volco"
        assert simulation.radius_increment == 0.1
        assert simulation.solver_tolerance == 0.0001
        assert simulation.x_crop == [11, 13]
        assert simulation.y_crop == ["all", "all"]
        assert simulation.z_crop == [0.0, "all"]
        assert simulation.consider_acceleration == False
        assert simulation.stl_ascii == True
        
    def test_should_use_default_values_when_not_provided(self):
        # Create a minimal config dictionary with only required fields
        minimal_config = {
            "voxel_size": 0.1,
            "step_size": 0.2,
            "simulation_name": "Test_defaults",
            "results_folder": "Results_defaults",
            "radius_increment": 0.1
        }
        
        # Create a mock printer with a nozzle diameter
        class MockPrinter:
            def __init__(self):
                self.nozzle_diameter = 0.4
                
        mock_printer = MockPrinter()
        
        # Initialize simulation with minimal config and mock printer
        simulation = Simulation(config_dict=minimal_config, printer=mock_printer)
        
        # Check that defaults are set correctly
        assert simulation.x_offset == mock_printer.nozzle_diameter * 5.0
        assert simulation.y_offset == mock_printer.nozzle_diameter * 5.0
        assert simulation.z_offset == 0
        assert simulation.sphere_z_offset == mock_printer.nozzle_diameter * 0.5
        assert simulation.solver_tolerance == 0.0001
        assert simulation.x_crop == ["all", "all"]
        assert simulation.y_crop == ["all", "all"]
        assert simulation.z_crop == ["all", "all"]
        assert simulation.consider_acceleration == False
        assert simulation.stl_ascii == False
