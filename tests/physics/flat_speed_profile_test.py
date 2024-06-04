import pytest
import numpy as np

from app.physics.flat_speed_profile import FlatSpeedProfile


class TestFlatSpeedProfile:
    def test_should_calculate_the_speed_profile(self):
        target_speed = 10
        threshold_speed = 40
        acceleration = 1200.0

        filament_length = 100

        final_time = filament_length / target_speed

        length_time = int(5)
        discrete_time = np.linspace(0, final_time, num=length_time)

        flat_speed_profile = FlatSpeedProfile()

        flat_speed_profile.calculate_displacements_in_time(
            discrete_time=discrete_time,
            final_time=final_time,
            target_speed=target_speed,
            threshold_speed=threshold_speed,
            acceleration=acceleration,
        )

        assert flat_speed_profile.displacements[-1] == 100
        assert flat_speed_profile.simulation_length == 100
