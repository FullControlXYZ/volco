import pytest
import numpy as np

from app.physics.extruder_speed import ExtruderSpeed


class TestExtruderSpeed:
    def test_should_calculate_the_speed_profile(self):
        threshold_speed = 50
        acceleration = 100.0

        extrusion_length = 100

        extruder_speed = ExtruderSpeed(
            extrusion_length=extrusion_length,
            threshold_speed=threshold_speed,
            acceleration=acceleration,
            total_time=1.25,
        )
        extruder_speed.calculate_displacements()

        assert extruder_speed.target_speed == 100.0

        assert (extruder_speed.speed_profile.displacements[-1] - 100.0) < 1e-3
