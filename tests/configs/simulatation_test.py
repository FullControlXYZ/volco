import pytest

from app.configs.simulation import Simulation


class TestSimulation:
    def test_should_load_simulation_settings(self):
        simulation = Simulation(config_path="tests/fixtures/simulation_settings.json")

        assert simulation.voxel_size == 0.05
        assert simulation.step_size == 0.2
        assert simulation.x_offset == 2.0
        assert simulation.y_offset == 2.0
        assert simulation.z_offset == 0.5
        assert simulation.simulation_name == "Test_volco"
        assert simulation.results_folder == "Results_volco"
        assert simulation.radius_increment == 0.1
        assert simulation.solver_tolerance == 0.0001
        assert simulation.x_crop == [11, 13]
        assert simulation.y_crop == ["all", "all"]
        assert simulation.z_crop == [0.0, "all"]
        assert simulation.consider_acceleration == False
        assert simulation.stl_ascii == True
