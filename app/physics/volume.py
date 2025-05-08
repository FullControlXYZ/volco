import math


class Volume:
    @staticmethod
    def get_volumes_for_filament(
        number_simulation_steps, total_volume, consider_acceleration=False,
        filament_length=None, printing_speed=None, printer=None,
    ):
        """Get volumes for a filament segment (handles both acceleration modes)
        
        There are two approaches to calculating filament volume distribution:
        
        1. Simple approach (consider_acceleration=False):
           Almost all users can simply use this method, which divides the total volume of
           filament evenly across simulation steps. This assumes uniform deposition along
           the length of the filament. The more complex second approach below is one example
           of a special case and there may be loads of different variants of it for completely
           different factors than acceleration.
           
        2. Acceleration-based approach (consider_acceleration=True):
           For users studying how volumetric deposition varies along the filament
           due to mismatches between XYZ movement acceleration and extruder motor
           feeding material into the nozzle. This requires additional parameters
           like filament_length, printing_speed, and printer specifications.
        """
        if not consider_acceleration:
            # Simple case: evenly distribute the volume across steps
            volume_per_step = total_volume / number_simulation_steps
            return [volume_per_step * (i + 1) for i in range(number_simulation_steps)]
        
        # Complex case: use acceleration-specific implementation
        from app.physics.acceleration.volume import AccelerationVolume
        return AccelerationVolume.get_volumes_for_filament(
            number_simulation_steps, total_volume,
            filament_length, printing_speed, printer
        )
