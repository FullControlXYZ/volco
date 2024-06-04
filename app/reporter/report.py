import logging
import os
import trimesh
import numpy as np

from app.configs.simulation import Simulation
from app.geometry.geometry_math import GeometryMath
from app.geometry.voxel_space import VoxelSpace


logger = logging.getLogger(__name__)


class Report:
    def __init__(self, voxel_space: VoxelSpace, simulation: Simulation):
        self._voxel_space = voxel_space
        self._simulation = simulation
        self.cropped_voxel_space = None

    def crop_voxel_space(self):
        crop_coordinates_for_axes = [
            self._simulation.x_crop,
            self._simulation.y_crop,
            self._simulation.z_crop,
        ]

        filament_translations = [
            self._voxel_space.filament_translations["x"],
            self._voxel_space.filament_translations["y"],
            0.0,
        ]

        axes_length = list(np.shape(self._voxel_space.space))

        indexes_to_crop = list()

        for (crop_coordinates, filament_translation, axis_length) in zip(
            crop_coordinates_for_axes, filament_translations, axes_length
        ):
            indexes_for_axis = self._crop_axis(
                crop_coordinates, filament_translation, axis_length
            )
            indexes_to_crop.append(indexes_for_axis)

        cropped_voxel_matrix = self._voxel_space.space.copy()
        self.cropped_voxel_space = cropped_voxel_matrix[
            indexes_to_crop[0][0] : indexes_to_crop[0][1] + 1,
            indexes_to_crop[1][0] : indexes_to_crop[1][1] + 1,
            indexes_to_crop[2][0] : indexes_to_crop[2][1] + 1,
        ]

    def generate_stl(self):
        result_path = self._get_result_folder_path()

        stl_file_name = self._simulation.simulation_name + ".stl"

        stl_path = os.path.join(result_path, stl_file_name)

        logger.info(f"[Reporter]: exporting stl to path {stl_path}...")

        origin = (0.0, 0.0, 0.0)
        scale = (
            self._simulation.voxel_size,
            self._simulation.voxel_size,
            self._simulation.voxel_size,
        )

        mesh = trimesh.voxel.base.VoxelGrid(self.cropped_voxel_space)
        mesh.origin[0], mesh.origin[1], mesh.origin[2] = origin
        mesh.scale.setflags(write=1)
        mesh.scale[0], mesh.scale[1], mesh.scale[2] = scale

        voxels = mesh.as_boxes(colors=None)
        _ = trimesh.exchange.export.export_mesh(voxels, stl_path, file_type="stl")
        logger.info(f"[Reporter]: stl exported!")

    def _crop_axis(self, crop_coordinates, filament_translation, axis_length):
        indexes_to_crop = [0, 0]

        for index in range(0, 2):
            crop_coordinate = crop_coordinates[index]

            if crop_coordinate == "all":
                indexes_to_crop[index] = index * (axis_length - 1)
            else:
                coordinate = crop_coordinate + filament_translation

                if coordinate < 0.0 and index == 0:
                    indexes_to_crop[index] = 0
                elif coordinate < 0.0 and index == 1:
                    indexes_to_crop[index] = axis_length - 1
                else:
                    candidate_index = GeometryMath.find_index(
                        coordinate, self._simulation.voxel_size
                    )
                    if candidate_index == -1:
                        candidate_index = 0

                    if candidate_index > axis_length - 1:
                        indexes_to_crop[1] = axis_length - 1
                    else:
                        indexes_to_crop[index] = candidate_index

        return indexes_to_crop

    def _get_result_folder_path(self):
        folder_name = self._simulation.results_folder

        cwd = os.getcwd()
        path = os.path.join(cwd, folder_name)

        if not os.path.exists(path):
            os.mkdir(path)
        return path
