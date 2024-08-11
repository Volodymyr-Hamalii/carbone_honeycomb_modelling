import numpy as np
from numpy import ndarray

import matplotlib.pyplot as plt
from matplotlib.pyplot import Figure, Axes  # type: ignore

from .lines_builder import LinesBuilder
from .visualization_parameters import VisualizationParameters, StructureVisualParameters
from ..utils import Logger


logger = Logger(__name__)


class StructureVisualizer:
    @classmethod
    def show_structure(
            cls,
            coordinates: ndarray,
            to_build_bonds: bool = True,
            visual_parameters: StructureVisualParameters = VisualizationParameters.carbone,
            num_of_min_distances: int = 3,
            skip_first_distances: int = 0
    ) -> None:

        # Prepare to visualize
        fig: Figure = plt.figure()
        ax: Axes = fig.add_subplot(111, projection='3d')

        x: ndarray = coordinates[:, 0]
        y: ndarray = coordinates[:, 1]
        z: ndarray = coordinates[:, 2]

        # Plot the atoms
        ax.scatter(x, y, z, s=visual_parameters.size, c=visual_parameters.color_atoms)  # type: ignore
        cls._set_equal_scaling(ax, x, y, z)

        if to_build_bonds:
            LinesBuilder.add_lines(
                coordinates=coordinates, ax=ax,
                color_bonds=visual_parameters.color_bonds,
                num_of_min_distances=num_of_min_distances,
                skip_first_distances=skip_first_distances)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')  # type: ignore

        plt.show()

    @classmethod
    def show_two_structures(
        cls,
        coordinates_first: ndarray,
        coordinates_second: ndarray,
        to_build_bonds: bool = False,
        visual_parameters_first: StructureVisualParameters = VisualizationParameters.carbone,
        visual_parameters_second: StructureVisualParameters = VisualizationParameters.al,
    ) -> None:
        """ Show 3D plot with 2 structures (by default there are carbone and aluminium) """

        # Prepare to visualize
        fig: Figure = plt.figure()
        ax: Axes = fig.add_subplot(111, projection='3d')

        # Plot first structure atoms (carbone by default)
        x_first: ndarray = coordinates_first[:, 0]
        y_first: ndarray = coordinates_first[:, 1]
        z_first: ndarray = coordinates_first[:, 2]
        cls._set_equal_scaling(ax, x_first, y_first, z_first)

        ax.scatter(
            x_first, y_first, z_first,
            c=visual_parameters_first.color_atoms,
            label=visual_parameters_first.label,
            s=visual_parameters_first.size,  # type: ignore
            alpha=visual_parameters_first.transparency)

        # Plot second structure atoms (al by default)
        x_second: ndarray = coordinates_second[:, 0]
        y_second: ndarray = coordinates_second[:, 1]
        z_second: ndarray = coordinates_second[:, 2]

        ax.scatter(
            x_second, y_second, z_second,
            c=visual_parameters_second.color_atoms,
            label=visual_parameters_second.label,
            s=visual_parameters_second.size,  # type: ignore
            alpha=visual_parameters_second.transparency)

        if to_build_bonds:
            # Carbone
            LinesBuilder.add_lines(
                coordinates=coordinates_first, ax=ax,
                color_bonds=visual_parameters_first.color_bonds,
                num_of_min_distances=3)

            # Aluminium
            LinesBuilder.add_lines(
                coordinates=coordinates_second, ax=ax,
                color_bonds=visual_parameters_second.color_bonds,
                num_of_min_distances=1,
                skip_first_distances=2)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')  # type: ignore
        ax.legend()

        plt.show()

    @staticmethod
    def _set_equal_scaling(ax: Axes, x_coor: ndarray, y_coor: ndarray, z_coor: ndarray) -> None:

        max_range = np.array([
            x_coor.max() - x_coor.min(),
            y_coor.max() - y_coor.min(),
            z_coor.max() - z_coor.min(),
        ]).max()

        mid_x = (x_coor.max() + x_coor.min()) * 0.5
        mid_y = (y_coor.max() + y_coor.min()) * 0.5
        mid_z = (z_coor.max() + z_coor.min()) * 0.5

        ax.set_xlim(mid_x - max_range/2, mid_x + max_range/2)
        ax.set_ylim(mid_y - max_range/2, mid_y + max_range/2)
        ax.set_zlim(mid_z - max_range/2, mid_z + max_range/2)  # type: ignore
