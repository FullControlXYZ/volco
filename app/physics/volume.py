import math


class Volume:
    @staticmethod
    def calculate_volumes_incremental(number_simulation_steps, total_volume):
        """Calculate incremental volumes (each step gets progressively larger)"""
        volume_per_step = total_volume / number_simulation_steps
        return [volume_per_step * (i + 1) for i in range(number_simulation_steps)]
    
    @staticmethod
    def get_volumes_for_filament(
        number_simulation_steps, total_volume, consider_acceleration=False,
        filament_length=None, printing_speed=None, printer=None,
    ):
        """Get volumes for a filament segment (handles both acceleration modes)"""
        if not consider_acceleration:
            # Simple case: incremental volume distribution
            return Volume.calculate_volumes_incremental(
                number_simulation_steps, total_volume
            )
        
        # Complex case: use acceleration-specific implementation
        from app.physics.acceleration.volume import AccelerationVolume
        return AccelerationVolume.get_volumes_for_filament(
            number_simulation_steps, total_volume,
            filament_length, printing_speed, printer
        )
