from app.physics.extruder_speed import ExtruderSpeed
from app.physics.nozzle_speed import NozzleSpeed

import numpy as np
import math


class Volume:
    @staticmethod
    def calculate_volumes_per_step(
        number_simulation_steps,
        step_size,
        filament_diameter,
        nozzle_profile: NozzleSpeed,
        extruder_profile: ExtruderSpeed,
    ):
        nozzle_pos = 0.0
        vol_deposited = list()

        for _ in range(0, number_simulation_steps):

            nozzle_pos += step_size

            disp_e = np.interp(
                nozzle_pos,
                nozzle_profile.speed_profile.displacements,
                extruder_profile.speed_profile.displacements,
            )

            vol_deposited.append(disp_e * math.pi * filament_diameter**2 / 4)

        return vol_deposited
