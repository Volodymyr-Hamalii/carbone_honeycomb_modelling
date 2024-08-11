from numpy import ndarray

import matplotlib.pyplot as plt
from matplotlib.pyplot import Figure, Axes  # type: ignore
from mpl_toolkits.mplot3d.art3d import Line3DCollection

from .lines_builder import LinesBuilder
from ..utils import Logger


logger = Logger(__name__)


class StructureVisualizer:
    @classmethod
    def show_structure(
            cls,
            coordinates: ndarray,
            to_build_bonds: bool = True,
            coordinates_color: str = "blue",
            ) -> None:

        # Prepare to visualize
        fig: Figure = plt.figure()
        ax: Axes = fig.add_subplot(111, projection='3d')

        x: ndarray = coordinates[:, 0]
        y: ndarray = coordinates[:, 1]
        z: ndarray = coordinates[:, 2]

        # Plot the atoms
        ax.scatter(x, y, z, s=100, c=coordinates_color)  # type: ignore

        if to_build_bonds:
            lines = LinesBuilder.build_lines(coordinates)
            lc = Line3DCollection(lines, colors='black', linewidths=1)
            ax.add_collection3d(lc)  # type: ignore

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')  # type: ignore

        plt.show()

    @classmethod
    def show_two_structures(
        cls,
        coordinates_first: ndarray,
        coordinates_second: ndarray,
        coordinates_first_color: str = "blue",
        coordinates_second_color: str = "red",
        coordinates_first_transparency: float = 0.75,
        coordinates_second_transparency: float = 0.75,
        to_build_bonds: bool = False
    ) -> None:
        # Prepare to visualize
        fig: Figure = plt.figure()
        ax: Axes = fig.add_subplot(111, projection='3d')

        # Plot first structure atoms
        x_first: ndarray = coordinates_first[:, 0]
        y_first: ndarray = coordinates_first[:, 1]
        z_first: ndarray = coordinates_first[:, 2]

        ax.scatter(
            x_first, y_first, z_first, c=coordinates_first_color,
            label='Carbon', s=50, alpha=coordinates_first_transparency)  # type: ignore

        # Plot second structure atoms
        x_second: ndarray = coordinates_second[:, 0]
        y_second: ndarray = coordinates_second[:, 1]
        z_second: ndarray = coordinates_second[:, 2]

        ax.scatter(
            x_second, y_second, z_second, c=coordinates_second_color,
            label='Aluminum', s=25, alpha=coordinates_second_transparency)  # type: ignore

        if to_build_bonds:
            lines_first = LinesBuilder.build_lines(coordinates_first)
            lines_second = LinesBuilder.build_lines(coordinates_second)
            lc = Line3DCollection(lines_first + lines_second, colors='black', linewidths=1)
            ax.add_collection3d(lc)  # type: ignore

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')  # type: ignore
        ax.legend()

        plt.show()

    @staticmethod
    def _filter_one_plane_coordinates(coordinates: ndarray, x: float = 0, y: float = 0) -> ndarray:
        return coordinates[(coordinates[:, 0] == x) | (coordinates[:, 1] == y)]
