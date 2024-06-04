from app.physics.speed_profile import SpeedProfile

import numpy as np


class FlatSpeedProfile(SpeedProfile):
    def calculate_displacements_in_time(
        self, discrete_time, final_time, target_speed, threshold_speed, acceleration
    ):
        speed = np.zeros(discrete_time.size)
        speed = np.add(speed, target_speed)

        displacements_nozzle, total_length = self.calculate_integral_trapezoidal(
            discrete_time, speed
        )

        self._displacements = displacements_nozzle
        self._simulation_length = total_length

    @property
    def displacements(self):
        return self._displacements

    @property
    def simulation_length(self):
        return self._simulation_length
