import json


class Printer:
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
        self.nozzle_jerk_speed = config["nozzle_jerk_speed"]
        self.extruder_jerk_speed = config["extruder_jerk_speed"]
        self.nozzle_acceleration = config["nozzle_acceleration"]
        self.extruder_acceleration = config["extruder_acceleration"]
        self.feedstock_filament_diameter = config["feedstock_filament_diameter"]
        self.nozzle_diameter = config["nozzle_diameter"]
