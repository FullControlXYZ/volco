from app.physics.flat_speed_profile import FlatSpeedProfile
from app.physics.speed import Speed
from app.physics.trapezoidal_speed_profile import TrapezoidalSpeedProfile

import numpy as np


class NozzleSpeed(Speed):
    def __init__(self, filament_length, target_speed, threshold_speed, acceleration):
        self._travel_length = filament_length
        self._target_speed = target_speed
        self._threshold_speed = threshold_speed
        self._acceleration = acceleration
        self._total_time = None
        self._speed_profile = None

    @property
    def travel_length(self):
        return self._travel_length

    @property
    def target_speed(self):
        return self._target_speed

    @property
    def threshold_speed(self):
        return self._threshold_speed

    @property
    def acceleration(self):
        return self._acceleration

    @property
    def total_time(self):
        return self._total_time

    @property
    def speed_profile(self):
        return self._speed_profile

    def calculate_displacements(self):
        self._find_total_time_and_speed_profile()

        length_time = int(1e5)
        time_vec = np.linspace(0, self.total_time, num=length_time)

        self.speed_profile.calculate_displacements_in_time(
            discrete_time=time_vec,
            final_time=self.total_time,
            target_speed=self.target_speed,
            threshold_speed=self.threshold_speed,
            acceleration=self.acceleration,
        )

    def _find_total_time_and_speed_profile(self):
        if self.target_speed < self.threshold_speed:
            self._total_time = self.travel_length / self.target_speed
            self._speed_profile = FlatSpeedProfile()
            return

        self._total_time = (
            self.target_speed / self.acceleration
            + self.travel_length / self.target_speed
            + self.threshold_speed**2 / self.target_speed / self.acceleration
            - 2 * self.threshold_speed / self.acceleration
        )

        self._speed_profile = TrapezoidalSpeedProfile()
        return
