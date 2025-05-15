import json


class Simulation:
    def __init__(self, config_path=None, config_dict=None, printer=None):
        self.printer = printer
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
        
        # Set default offsets based on nozzle diameter if printer is available
        nozzle_diameter_factor = 5.0
        default_offset = 2.0  # Default if printer not available
        
        if self.printer:
            default_offset = self.printer.nozzle_diameter * nozzle_diameter_factor
        
        self.x_offset = config.get("x_offset", default_offset)
        self.y_offset = config.get("y_offset", default_offset)
        self.z_offset = config.get("z_offset", 0)
        
        # New parameter for sphere z-offset with default of half nozzle diameter
        default_sphere_z_offset = 0.5
        if self.printer:
            default_sphere_z_offset = self.printer.nozzle_diameter * 0.5
        
        self.sphere_z_offset = config.get("sphere_z_offset", default_sphere_z_offset)
        
        self.simulation_name = config["simulation_name"]
        self.results_folder = config["results_folder"]
        self.radius_increment = config["radius_increment"]
        
        # Default solver tolerance
        self.solver_tolerance = config.get("solver_tolerance", 0.0001)
        
        # Default crop values to "all"
        self.x_crop = config.get("x_crop", ["all", "all"])
        self.y_crop = config.get("y_crop", ["all", "all"])
        self.z_crop = config.get("z_crop", ["all", "all"])
        
        # Default acceleration and STL settings
        self.consider_acceleration = config.get("consider_acceleration", False)
        self.stl_ascii = config.get("stl_ascii", False)
