import pytest

from app.physics.nozzle_speed import NozzleSpeed


class TestNozzleSpeed:
    def test_should_calculate_the_speed_profile(self):
        target_speed = 10
        threshold_speed = 40
        acceleration = 1200.0

        filament_length = 100

        nozzle_speed = NozzleSpeed(
            filament_length=filament_length,
            target_speed=target_speed,
            threshold_speed=threshold_speed,
            acceleration=acceleration,
        )
        nozzle_speed.calculate_displacements()

        assert (nozzle_speed.speed_profile.displacements[-1] - 100.0) < 1e-3
