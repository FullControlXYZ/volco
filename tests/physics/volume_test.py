import pytest
from app.physics.extruder_speed import ExtruderSpeed
from app.physics.nozzle_speed import NozzleSpeed
from app.physics.volume import Volume

import math


class TestVolume:
    def test_should_calculate_the_speed_profile(self):
        target_speed = 10
        threshold_speed = 40
        acceleration = 1200.0

        filament_length = 100.0

        extrusion_length = 10.0

        nozzle_speed = NozzleSpeed(
            filament_length=filament_length,
            target_speed=target_speed,
            threshold_speed=threshold_speed,
            acceleration=acceleration,
        )
        nozzle_speed.calculate_displacements()

        extruder_speed = ExtruderSpeed(
            extrusion_length=extrusion_length,
            threshold_speed=threshold_speed,
            acceleration=acceleration,
            total_time=nozzle_speed.total_time,
        )
        extruder_speed.calculate_displacements()

        volumes = Volume.calculate_volumes_per_step(
            number_simulation_steps=5,
            step_size=20,
            filament_diameter=1.75,
            nozzle_profile=nozzle_speed,
            extruder_profile=extruder_speed,
        )

        expected_volume = extrusion_length * math.pi * 1.75**2 / 4

        assert volumes[-1] == expected_volume
