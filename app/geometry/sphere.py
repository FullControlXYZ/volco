import math
import numpy as np

from app.geometry.geometry_math import GeometryMath
from app.solvers.bisection_method import BisectionMethod


class Sphere:
    def __init__(self, centre_coordinates, voxel_size):
        self.centre_coordinates = centre_coordinates
        self.voxel_size = voxel_size

    def deposit_sphere(
        self,
        voxel_space,
        nozzle_height,
        sphere_volume,
        voxel_space_target_volume,
        solver_tolerance,
        radius_increment,
    ):
        initial_radius = self.estimate_initial_radius(sphere_volume)

        _, voxel_space = BisectionMethod().execute(
            self._deposit_sphere,
            initial_point=initial_radius,
            tolerance=solver_tolerance,
            increment=radius_increment,
            fun_increase_tolerance=self._increase_solver_tolerance,
            args=(
                voxel_space,
                nozzle_height,
                voxel_space_target_volume,
            ),
        )

        return voxel_space

    def fill_voxels(self, voxel_space, radius, lower_indexes, upper_indexes):
        empty_voxels = GeometryMath.find_empty_voxels_in_space(
            voxel_space, lower_indexes, upper_indexes
        )

        for voxel in empty_voxels:
            voxel_coordinate = GeometryMath.find_coordinates(voxel, self.voxel_size)

            distance_to_centre = GeometryMath.distance(
                voxel_coordinate, self.centre_coordinates
            )

            if distance_to_centre <= radius + self.voxel_size * 1e-8:
                voxel_space.itemset(tuple(voxel), 1)

        return voxel_space

    def estimate_initial_radius(self, volume):
        return (3.0 * volume / (4.0 * math.pi)) ** (1.0 / 3.0)

    def find_sphere_limits(self, radius, nozzle_height):
        # subtract 1 because the voxel space starts at position [0,0]
        min_indexes = [
            self._find_index(centre_coordinate - radius)
            for centre_coordinate in self.centre_coordinates
        ]

        if min_indexes[2] < 0:
            min_indexes[2] = 0

        max_indexes = [
            self._find_index(centre_coordinate + radius)
            for centre_coordinate in self.centre_coordinates
        ]

        _, _, z0 = self.centre_coordinates
        if z0 + radius > nozzle_height:
            max_indexes[2] = self._find_index(nozzle_height)
        else:
            max_indexes[2] = self._find_index(z0 + radius)

        return min_indexes, max_indexes

    """
    Checks if the sphere violates the boundaries of the voxel space. If it does,
    the voxel space is expanded to accommodate the sphere.
    """

    def deform_voxel_space_for_big_spheres(self, voxel_space, radius):
        max_indexes = [
            self._find_index(coord + radius) for coord in self.centre_coordinates
        ]

        for axis_number in range(0, 3):
            voxel_space = self._maybe_expand_voxel_space(
                voxel_space, max_indexes[axis_number], axis_number
            )

        return voxel_space

    def _maybe_expand_voxel_space(self, voxel_space, max_index, axis_number):
        size = np.shape(voxel_space)

        index_size = size[axis_number]

        if max_index < index_size:
            return voxel_space

        number_to_be_added = max_index - index_size + 1

        new_size = list(size)
        new_size[axis_number] = number_to_be_added

        mat_add = np.zeros(new_size, dtype=np.int8)

        voxel_space = np.concatenate((voxel_space, mat_add), axis=axis_number)

        return voxel_space

    def _deposit_sphere(self, radius, voxel_space, nozzle_height, target_volume):
        copy_voxel_space = voxel_space.copy()

        copy_voxel_space = self.deform_voxel_space_for_big_spheres(
            copy_voxel_space, radius
        )

        lower_indexes, upper_indexes = self.find_sphere_limits(radius, nozzle_height)

        copy_voxel_space = self.fill_voxels(
            copy_voxel_space, radius, lower_indexes, upper_indexes
        )

        current_volume = GeometryMath.calculate_filled_volume(
            copy_voxel_space, self.voxel_size
        )

        volume_overshoot = current_volume / target_volume - 1.0

        return volume_overshoot, copy_voxel_space

    def _increase_solver_tolerance(self, radius_a, radius_b):
        return radius_b - radius_a < self.voxel_size * 0.5

    def _find_index(self, coordinate):
        return GeometryMath.find_index(coordinate, self.voxel_size)
