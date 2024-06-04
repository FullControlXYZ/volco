import pytest
import numpy as np

from app.geometry.geometry_math import GeometryMath


class TestGeometryMath:
    def test_should_calculate_distance(self):
        point_a = [1.0, -1.0, 3.0]
        point_b = [4.0, 5.0, 2.0]

        assert GeometryMath.distance(point_a, point_b) == (9 + 36 + 1) ** 0.5

    def test_should_calculate_the_direction_vector(self):
        point_a = [4.0, 2.0, 5.0]
        point_b = [3.0, 2.0, 5.0]

        assert GeometryMath.direction_vector(point_a, point_b) == [-1.0, 0.0, 0.0]


class TestFindEmptyVoxels:
    def test_should_return_an_empty_list(self):
        voxel_space = np.ones((10, 10, 10), dtype=np.int8)

        empty_voxels = GeometryMath.find_empty_voxels_in_space(
            voxel_space, [0, 0, 0], [9, 9, 9]
        )

        assert empty_voxels == []

    def test_should_return_a_list_of_empty_voxels(self):
        voxel_space = np.ones((10, 10, 10), dtype=np.int8)

        voxel_space[4][4][4] = 0

        empty_voxels = GeometryMath.find_empty_voxels_in_space(
            voxel_space, [0, 0, 0], [9, 9, 9]
        )

        assert empty_voxels == [[4, 4, 4]]


class TestFindCoordinates:
    def test_should_find_coordinates_based_on_indexes(self):
        indexes = [1, 2, 3]
        voxel_size = 2.0

        assert GeometryMath.find_coordinates(indexes, voxel_size) == [3.0, 5.0, 7.0]
