from numpy import ndarray

import matplotlib.pyplot as plt
from matplotlib.pyplot import Figure, Axes  # type: ignore

from .lines_builder import LinesBuilder
from .visualization_parameters import VisualizationParameters
from ..utils import Logger


logger = Logger(__name__)


class StructureVisualizer:
    @classmethod
    def show_structure(
            cls,
            coordinates: ndarray,
            to_build_bonds: bool = True,
            color_atoms: str = VisualizationParameters.carbone.color_atoms,
            color_bonds: str = VisualizationParameters.carbone.color_bonds,
    ) -> None:

        # Prepare to visualize
        fig: Figure = plt.figure()
        ax: Axes = fig.add_subplot(111, projection='3d')

        x: ndarray = coordinates[:, 0]
        y: ndarray = coordinates[:, 1]
        z: ndarray = coordinates[:, 2]

        # Plot the atoms
        ax.scatter(x, y, z, s=VisualizationParameters.carbone.size, c=color_atoms)  # type: ignore

        if to_build_bonds:
            LinesBuilder.add_lines(
                coordinates=coordinates, ax=ax,
                color_bonds=color_bonds)

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

        # color_atoms_first: str = VisualizationParameters.C_ATOMS,
        # color_bonds_first: str = VisualizationParameters.C_BONDS,
        # color_first_transparency: float = 0.75,

        # color_atoms_second: str = VisualizationParameters.AL_ATOMS,
        # color_bonds_second: str = VisualizationParameters.AL_BONDS,
        # color_second_transparency: float = 0.75,

    ) -> None:
        # Prepare to visualize
        fig: Figure = plt.figure()
        ax: Axes = fig.add_subplot(111, projection='3d')

        # Plot first structure atoms
        x_first: ndarray = coordinates_first[:, 0]
        y_first: ndarray = coordinates_first[:, 1]
        z_first: ndarray = coordinates_first[:, 2]

        ax.scatter(
            x_first, y_first, z_first, c=VisualizationParameters.carbone.color_atoms,
            label='Carbon', s=VisualizationParameters.carbone.size,  # type: ignore
            alpha=VisualizationParameters.carbone.transparency)

        # Plot second structure atoms
        x_second: ndarray = coordinates_second[:, 0]
        y_second: ndarray = coordinates_second[:, 1]
        z_second: ndarray = coordinates_second[:, 2]

        ax.scatter(
            x_second, y_second, z_second, c=VisualizationParameters.al.color_atoms,
            label='Aluminum', s=VisualizationParameters.al.size,  # type: ignore
            alpha=VisualizationParameters.al.transparency)

        if to_build_bonds:
            LinesBuilder.add_lines(
                coordinates=coordinates_first, ax=ax,
                color_bonds=VisualizationParameters.carbone.color_bonds)

            LinesBuilder.add_lines(
                coordinates=coordinates_second, ax=ax,
                color_bonds=VisualizationParameters.al.color_bonds,
                num_of_min_distances=1)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')  # type: ignore
        ax.legend()

        plt.show()

    @staticmethod
    def _filter_one_plane_coordinates(coordinates: ndarray, x: float = 0, y: float = 0) -> ndarray:
        return coordinates[(coordinates[:, 0] == x) | (coordinates[:, 1] == y)]
