from abc import ABC, abstractmethod, abstractproperty

import numpy as np


class SpeedProfile(ABC):
    @abstractmethod
    def calculate_displacements_in_time(
        self, discrete_time, final_time, target_speed, threshold_speed, acceleration
    ):
        pass

    def calculate_integral_trapezoidal(self, x, y):
        len_x = len(x)

        intg = 0.0

        intg_vec = np.zeros(len_x)

        for i in range(1, len_x):
            intg = (y[i] + y[i - 1]) * (x[i] - x[i - 1]) * 0.5 + intg

            intg_vec[i] = intg

        return intg_vec, intg

    @abstractproperty
    def displacements(self):
        pass

    @abstractproperty
    def simulation_length(self):
        pass
