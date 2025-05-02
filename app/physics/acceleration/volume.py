import math
import numpy as np

from app.physics.acceleration.acceleration_profiles import AccelerationManager
from app.physics.acceleration.extruder_speed import ExtruderSpeed
from app.physics.acceleration.nozzle_speed import NozzleSpeed


class AccelerationVolume:
    """
    Handles volume calculations that consider acceleration.
    This class is separate from the main Volume class to maintain a clean separation
    between acceleration and non-acceleration code.
    """
    
    @staticmethod
    def calculate_volumes_with_acceleration(
        number_simulation_steps,
        step_size,
        feedstock_filament_diameter,
        nozzle_profile,
        extruder_profile,
    ):
        """
        Calculate volumes considering acceleration profiles.
        
        Args:
            number_simulation_steps: Number of simulation steps
            step_size: Size of each simulation step
            feedstock_filament_diameter: Diameter of the filament
            nozzle_profile: Nozzle speed profile
            extruder_profile: Extruder speed profile
            
        Returns:
            List of volumes for each step
        """
        nozzle_pos = 0.0
        vol_deposited = list()

        for _ in range(0, number_simulation_steps):
            nozzle_pos += step_size

            disp_e = np.interp(
                nozzle_pos,
                nozzle_profile.speed_profile.displacements,
                extruder_profile.speed_profile.displacements,
            )

            vol_deposited.append(disp_e * math.pi * feedstock_filament_diameter**2 / 4)

        return vol_deposited
    
    @staticmethod
    def get_volumes_for_filament(
        number_simulation_steps,
        total_volume,
        filament_length,
        printing_speed,
        printer,
    ):
        """
        Calculate volumes for a filament segment considering acceleration.
        
        Args:
            number_simulation_steps: Number of simulation steps
            total_volume: Total volume to distribute
            filament_length: Length of the filament
            printing_speed: Printing speed
            printer: Printer configuration
            
        Returns:
            List of volumes for each step
        """
        if None in (filament_length, printing_speed, printer):
            raise ValueError("filament_length, printing_speed, and printer must be provided for acceleration-based volume calculation")
        
        # Calculate speed profiles
        nozzle_speed, extruder_speed = AccelerationManager.calculate_speed_profiles(
            filament_length=filament_length,
            volume=total_volume,
            printing_speed=printing_speed,
            printer=printer,
        )
        
        # Calculate volumes with acceleration
        return AccelerationVolume.calculate_volumes_with_acceleration(
            number_simulation_steps=number_simulation_steps,
            step_size=filament_length / number_simulation_steps,
            feedstock_filament_diameter=printer.feedstock_filament_diameter,
            nozzle_profile=nozzle_speed,
            extruder_profile=extruder_speed,
        )