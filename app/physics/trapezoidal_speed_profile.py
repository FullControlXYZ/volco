from app.physics.speed_profile import SpeedProfile

import numpy as np


class TrapezoidalSpeedProfile(SpeedProfile):
    def calculate_displacements_in_time(
        self, discrete_time, final_time, target_speed, threshold_speed, acceleration
    ):
        ascending_time = (target_speed - threshold_speed) / acceleration

        cruising_time = final_time - ascending_time

        speed = np.zeros(discrete_time.size)

        cont = -1

        for time_i in discrete_time:
            cont += 1

            if time_i <= ascending_time:
                speed[cont] = acceleration * time_i + threshold_speed
            elif time_i > ascending_time and time_i <= cruising_time:
                speed[cont] = target_speed
            elif time_i > cruising_time:
                speed[cont] = (
                    -acceleration * time_i + threshold_speed + acceleration * final_time
                )

        displacements, total_length = self.calculate_integral_trapezoidal(
            discrete_time, speed
        )

        self._displacements = displacements
        self._simulation_length = total_length

    @property
    def displacements(self):
        return self._displacements

    @property
    def simulation_length(self):
        return self._simulation_length
