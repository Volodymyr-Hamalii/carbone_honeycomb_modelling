import numpy as np
from numpy import ndarray

import matplotlib.pyplot as plt
from matplotlib.pyplot import Figure, Axes  # type: ignore

from .lines_builder import LinesBuilder
from .visualization_parameters import VisualizationParameters, StructureVisualParameters
from ..utils import Logger


logger = Logger("StructureVisualizer")


class StructureVisualizer:
    @classmethod
    def show_structure(
            cls,
            coordinates: ndarray,
            to_build_bonds: bool = True,
            visual_parameters: StructureVisualParameters = VisualizationParameters.carbon,
            num_of_min_distances: int = 3,
            skip_first_distances: int = 0,
            set_equal_scale: bool | None = None,
            title: str | None = None,
    ) -> None:
        """ Show 3D plot with 1 structure. """

        # Prepare to visualize
        fig: Figure = plt.figure()
        ax: Axes = fig.add_subplot(111, projection='3d')

        cls._plot_atoms_3d(
            ax=ax,
            coordinates=coordinates,
            visual_parameters=visual_parameters,
            to_build_bonds=to_build_bonds,
            num_of_min_distances=num_of_min_distances,
            skip_first_distances=skip_first_distances,
            set_equal_scale=set_equal_scale,
        )

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')  # type: ignore

        if title is not None:
            ax.set_title(title)

        plt.show()

    @classmethod
    def show_two_structures(
        cls,
        coordinates_first: ndarray,
        coordinates_second: ndarray,
        to_build_bonds: bool = False,
        visual_parameters_first: StructureVisualParameters = VisualizationParameters.carbon,
        visual_parameters_second: StructureVisualParameters = VisualizationParameters.al,
        title: str | None = None,
    ) -> None:
        """ Show 3D plot with 2 structures (by default there are carbon and aluminium) """

        # Prepare to visualize
        fig: Figure = plt.figure()
        ax: Axes = fig.add_subplot(111, projection='3d')

        # Plot first structure atoms (carbon by default)
        cls._plot_atoms_3d(
            ax=ax,
            coordinates=coordinates_first,
            visual_parameters=visual_parameters_first,
            to_build_bonds=to_build_bonds)

        # Plot second structure atoms (al by default)
        cls._plot_atoms_3d(
            ax=ax,
            coordinates=coordinates_second,
            visual_parameters=visual_parameters_second,
            to_build_bonds=to_build_bonds,
            num_of_min_distances=1,
            skip_first_distances=0)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')  # type: ignore
        ax.legend()

        if title is not None:
            ax.set_title(title)

        plt.show()


    @classmethod
    def show_2d_graph(
        cls,
        coordinates: np.ndarray,
        title: str | None = None,
        visual_parameters: StructureVisualParameters = VisualizationParameters.carbon,
        show_coordinates: bool = False,
    ) -> None:
        # Prepare to visualize in 2D
        fig: Figure = plt.figure()
        ax: Axes = fig.add_subplot(111)  # No 3D projection here, just 2D

        # Plot points
        x: np.ndarray = coordinates[:, 0]
        y: np.ndarray = coordinates[:, 1]
        ax.scatter(x, y, color=visual_parameters.color_atoms, label='Points')

        # If requested, show coordinates near each point
        if show_coordinates:
            for xx, yy in zip(x, y):
                ax.annotate(
                    f"({xx:.2f}, {yy:.2f})",
                    (xx, yy),
                    textcoords="offset points",
                    xytext=(5, 5),  # Offset from the point
                    fontsize=6,
                )

        ax.set_xlabel('X')
        ax.set_ylabel('Y')

        if title is not None:
            ax.set_title(title)

        ax.legend()
        plt.grid(True)
        plt.show()

    @classmethod
    def _plot_atoms_3d(
            cls,
            ax: Axes,
            coordinates: ndarray,
            visual_parameters: StructureVisualParameters,
            set_equal_scale: bool | None = None,
            to_build_bonds: bool = True,
            num_of_min_distances: int = 3,
            skip_first_distances: int = 0,
    ) -> None:

        if coordinates.size == 0:
            logger.warning(f"No points to plot for {visual_parameters.label}.")
            return

        x_first: ndarray = coordinates[:, 0]
        y_first: ndarray = coordinates[:, 1]
        z_first: ndarray = coordinates[:, 2]

        ax.scatter(
            x_first, y_first, z_first,
            c=visual_parameters.color_atoms,
            label=visual_parameters.label,
            s=visual_parameters.size,  # type: ignore
            alpha=visual_parameters.transparency)

        if set_equal_scale is None:
            set_equal_scale = visual_parameters.set_equal_scale

        if set_equal_scale:
            cls._set_equal_scale(ax, x_first, y_first, z_first)

        if to_build_bonds:
            # Carbon
            LinesBuilder.add_lines(
                coordinates=coordinates, ax=ax,
                color_bonds=visual_parameters.color_bonds,
                num_of_min_distances=num_of_min_distances,
                skip_first_distances=skip_first_distances)

    @staticmethod
    def _set_equal_scale(ax: Axes, x_coor: ndarray, y_coor: ndarray, z_coor: ndarray) -> None:
        """Set equal scaling for all axis."""

        x_min, x_max = x_coor.min(), x_coor.max()
        y_min, y_max = y_coor.min(), y_coor.max()
        z_min, z_max = z_coor.min(), z_coor.max()

        min_lim = np.min([x_min, y_min, z_min])
        max_lim = np.max([x_max, y_max, z_max])

        delta = (max_lim + min_lim) / 2
        delta_plus = delta * 1.25  # to stretch specific axis

        x_mid = (x_max + x_min) / 2
        y_mid = (y_max + y_min) / 2
        z_mid = (z_max + z_min) / 2

        ax.set_xlim(x_mid - delta_plus, x_mid + delta_plus)
        ax.set_ylim(y_mid - delta_plus, y_mid + delta_plus)
        ax.set_zlim(z_mid - delta, z_mid + delta)  # type: ignore
