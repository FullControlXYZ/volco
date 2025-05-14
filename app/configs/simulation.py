import json


class Simulation:
    def __init__(self, config_path=None, config_dict=None):
        if config_path:
            self._read_config_file(config_path)
        elif config_dict:
            self._load_config_from_dict(config_dict)
        else:
            raise ValueError("Either config_path or config_dict must be provided")

    def _read_config_file(self, config_path):
        with open(config_path) as f:
            config = json.load(f)
            self._load_config_from_dict(config)

    def _load_config_from_dict(self, config):
        self.voxel_size = config["voxel_size"]
        self.step_size = config["step_size"]
        self.x_offset = config["x_offset"]
        self.y_offset = config["y_offset"]
        self.z_offset = config["z_offset"]
        self.simulation_name = config["simulation_name"]
        self.results_folder = config["results_folder"]
        self.radius_increment = config["radius_increment"]
        self.solver_tolerance = config["solver_tolerance"]
        self.x_crop = config["x_crop"]
        self.y_crop = config["y_crop"]
        self.z_crop = config["z_crop"]
        self.consider_acceleration = config["consider_acceleration"]
        self.stl_ascii = config["stl_ascii"]
