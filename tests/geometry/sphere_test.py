import pytest
import numpy as np

from app.geometry.sphere import Sphere


class TestSphere:
    def test_should_calculate_the_limits(self):
        sphere = Sphere(centre_coordinates=[2.0, 4.0, 0.099], voxel_size=0.05)

        min_coordinates, max_coordinates = sphere.find_sphere_limits(0.15, 0.3)

        assert min_coordinates == [36, 76, 0]
        assert max_coordinates == [42, 82, 4]

    def test_should_deform_voxel_space(self):
        voxel_space = np.zeros((10, 15, 5), dtype=np.int8)

        sphere = Sphere(centre_coordinates=[5.0, 10.0, 3.0], voxel_size=1.0)

        new_voxel_space = sphere.deform_voxel_space_for_big_spheres(voxel_space, 10.0)

        length_i, length_j, length_k = new_voxel_space.shape

        assert length_i == 15
        assert length_j == 20
        assert length_k == 13


class TestFillVoxes:
    def test_fill_voxels_inside_sphere(self):
        voxel_space = np.zeros((10, 10, 10), dtype=np.int8)

        sphere = Sphere(centre_coordinates=[4.0, 4.0, 4.0], voxel_size=1.0)

        voxel_space = sphere.fill_voxels(voxel_space, 2.0, [0, 0, 0], [9, 9, 9])

        assert voxel_space[4][4][4] == 1
