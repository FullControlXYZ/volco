import math
import numpy as np

from app.physics.acceleration.nozzle_speed import NozzleSpeed
from app.physics.acceleration.extruder_speed import ExtruderSpeed


class AccelerationManager:
    """
    Manages the calculation of speed profiles considering acceleration.
    """

    @staticmethod
    def calculate_speed_profiles(
        filament_length,
        volume,
        printing_speed,
        printer,
    ):
        """
        Calculate nozzle and extruder speed profiles considering acceleration.
        
        Args:
            filament_length: Length of the filament segment
            volume: Volume to be deposited
            printing_speed: Target printing speed
            printer: Printer configuration
            
        Returns:
            Tuple of (nozzle_speed, extruder_speed) profiles
        """
        # Calculate nozzle speed profile
        nozzle_speed = NozzleSpeed(
            filament_length=filament_length,
            target_speed=printing_speed,
            threshold_speed=printer.nozzle_jerk_speed,
            acceleration=printer.nozzle_acceleration,
        )
        nozzle_speed.calculate_displacements()

        # Calculate extruder speed profile
        extruder_speed = ExtruderSpeed(
            volume=volume,
            threshold_speed=printer.extruder_jerk_speed,
            acceleration=printer.extruder_acceleration,
            total_time=nozzle_speed.total_time,
            printer=printer,
        )
        extruder_speed.calculate_displacements()

        return nozzle_speed, extruder_speed