import pytest
import numpy as np

from app.physics.trapezoidal_speed_profile import TrapezoidalSpeedProfile


class TestTrapezoidalSpeedProfile:
    def test_should_calculate_the_speed_profile(self):
        target_speed = 100
        threshold_speed = 50
        acceleration = 100.0

        filament_length = 100

        final_time = (
            target_speed / acceleration
            + filament_length / target_speed
            + threshold_speed**2 / target_speed / acceleration
            - 2 * threshold_speed / acceleration
        )

        length_time = int(6)
        discrete_time = np.linspace(0, final_time, num=length_time)

        trapezoidal_speed_profile = TrapezoidalSpeedProfile()

        trapezoidal_speed_profile.calculate_displacements_in_time(
            discrete_time=discrete_time,
            final_time=final_time,
            target_speed=target_speed,
            threshold_speed=threshold_speed,
            acceleration=acceleration,
        )

        assert trapezoidal_speed_profile.displacements[0] == 0.0
        assert trapezoidal_speed_profile.displacements[2] == (100.0 + 50.0) * 0.5 / 2.0
        assert (
            trapezoidal_speed_profile.displacements[3]
            == trapezoidal_speed_profile.displacements[2] + 100.0 * 0.25
        )

        assert trapezoidal_speed_profile.displacements[-1] == 100
        assert trapezoidal_speed_profile.simulation_length == 100
