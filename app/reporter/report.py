import logging
import os
import numpy as np
import trimesh

from app.configs.simulation import Simulation
from app.geometry.geometry_math import GeometryMath
from app.geometry.voxel_space import VoxelSpace
from app.reporter.visualization import color_mesh, visualize_with_trimesh, visualize_with_plotly
from app.reporter.mesh import generate_mesh_from_voxels, export_mesh_to_stl


logger = logging.getLogger(__name__)


class SimulationOutput:
    """
    Processes and visualizes the output of a 3D printing simulation.
    
    This class handles post-processing of simulation data, including:
    - Cropping the voxel space
    - Generating meshes from voxel data
    - Exporting meshes to STL files
    - Visualizing the results
    """
    def __init__(self, voxel_space: VoxelSpace, simulation: Simulation):
        self._voxel_space = voxel_space
        self._simulation = simulation
        self.cropped_voxel_space = None
        self.mesh = None

    def crop_voxel_space(self):
        """
        Crop the voxel space based on the simulation configuration.
        """
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

    def generate_mesh(self):
        """
        Generate a 3D mesh from the voxel data using marching cubes algorithm.
        Returns the mesh object for further use.
        """
        if self.cropped_voxel_space is None:
            logger.warning("[SimulationOutput]: No cropped voxel space available. Run crop_voxel_space() first.")
            return None

        # Get voxel size
        voxel_size = self._simulation.voxel_size

        # Create a mesh from the voxel data using the mesh module
        mesh = generate_mesh_from_voxels(self.cropped_voxel_space, voxel_size)

        # Store the mesh for later use
        self.mesh = mesh
        return mesh

    def color_mesh(self, mesh=None, color_scheme='cyan_blue'):
        """
        Apply colors to the mesh based on height (z-value).
        Returns the colored mesh.
        
        This method uses the color_mesh function from the visualization submodule.
        """
        if mesh is None:
            if self.mesh is None:
                logger.warning("[SimulationOutput]: No mesh available. Generate mesh first.")
                return None
            mesh = self.mesh

        # Check if mesh is a Scene object (from box representation)
        if isinstance(mesh, trimesh.Scene):
            # For box representation, we can't easily color individual vertices
            # Return the original scene without coloring
            return mesh
            
        # For Trimesh objects (from marching cubes)
        if not hasattr(mesh, 'vertices') or len(mesh.vertices) == 0:
            logger.warning("[SimulationOutput]: Mesh has no vertices to color.")
            return mesh
            
        return color_mesh(mesh, color_scheme)

    def export_mesh_to_stl(self, mesh=None):
        """
        Export the mesh to an STL file.
        If mesh is not provided, uses the stored mesh.
        If file_path is not provided, uses the default path based on simulation settings.
        """
        if mesh is None:
            if self.mesh is None:
                logger.warning("[SimulationOutput]: No mesh available. Generate mesh first.")
                return
            mesh = self.mesh
            
        
        result_path = self._get_result_folder_path()
        stl_file_name = self._simulation.simulation_name + ".stl"
        file_path = os.path.join(result_path, stl_file_name)
        
        return export_mesh_to_stl(mesh, file_path)

    def visualize_mesh(self, mesh=None, visualizer='trimesh', color_scheme='cyan_blue'):
        """
        Create a 3D visualization of the mesh with coloring applied.
        If mesh is not provided, uses the stored mesh.
        
        Parameters:
        -----------
        mesh : trimesh.Trimesh
            The mesh to visualize
        visualizer : str
            The visualizer to use ('trimesh' or 'plotly')
        color_scheme : str
            The color scheme to use ('cyan_blue', 'viridis', or None for no coloring)
            
        Returns:
        --------
        trimesh.Scene or plotly.graph_objects.Figure
            A visualization object that can be displayed
        """
        if mesh is None:
            if self.mesh is None:
                logger.warning("[SimulationOutput]: No mesh available. Generate mesh first.")
                return None
            mesh = self.mesh
        
        # Apply coloring if requested
        if color_scheme is not None:
            mesh = self.color_mesh(mesh, color_scheme=color_scheme)
            
        if visualizer.lower() == 'plotly':
            return visualize_with_plotly(mesh)
        else:  # default to trimesh
            return visualize_with_trimesh(mesh)

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
