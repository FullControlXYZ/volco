from app.physics.flat_speed_profile import FlatSpeedProfile
from app.physics.speed import Speed
from app.physics.trapezoidal_speed_profile import TrapezoidalSpeedProfile

import numpy as np
import math


class ExtruderSpeed(Speed):
    def __init__(self, volume, threshold_speed, acceleration, total_time, printer):
        # Convert volume to extrusion length
        filament_area = math.pi * (printer.bulk_filament_diameter/2)**2
        self._travel_length = volume / filament_area
        
        self._target_speed = None
        self._threshold_speed = threshold_speed
        self._acceleration = acceleration
        self._total_time = total_time
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
        self._find_target_speed_and_speed_profile()

        length_time = int(1e5)
        time_vec = np.linspace(0, self.total_time, num=length_time)

        self.speed_profile.calculate_displacements_in_time(
            discrete_time=time_vec,
            final_time=self.total_time,
            target_speed=self.target_speed,
            threshold_speed=self.threshold_speed,
            acceleration=self.acceleration,
        )

    def _find_target_speed_and_speed_profile(self):
        if self.travel_length / self.total_time < self.threshold_speed:
            self._speed_profile = FlatSpeedProfile()
            self._target_speed = self.travel_length / self.total_time
            return

        coeff = [
            1.0 / self.acceleration,
            -self.total_time - 2.0 * self.threshold_speed / self.acceleration,
            self.threshold_speed**2 / self.acceleration + self.travel_length,
        ]
        possible_target_speeds = np.roots(coeff)

        t1 = (possible_target_speeds - self.threshold_speed) / self.acceleration
        t2 = (
            self.total_time
            - 2.0 * possible_target_speeds / self.acceleration
            + 2.0 * self.threshold_speed / self.acceleration
            + t1
        )

        for (t1_i, t2_i, Ve_i) in zip(t1, t2, possible_target_speeds):
            if t2_i > t1_i and t1_i > 0 and t2_i > 0:
                self._target_speed = Ve_i

        self._speed_profile = TrapezoidalSpeedProfile()
        return
