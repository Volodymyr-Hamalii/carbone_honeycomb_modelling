import numpy as np
from numpy import ndarray

import matplotlib.pyplot as plt
from matplotlib.pyplot import Figure, Axes  # type: ignore
from matplotlib.collections import PathCollection

# from src.coordinate_operations import DistanceMeasure
from .lines_builder import LinesBuilder
from .visualization_params import VisualizationParams, StructureVisualParams
from ..utils import Logger


logger = Logger("StructureVisualizer")


class StructureVisualizer:
    @classmethod
    def show_structure(
            cls,
            coordinates: ndarray,
            to_build_bonds: bool = True,
            visual_params: StructureVisualParams = VisualizationParams.carbon,
            num_of_min_distances: int = 2,
            skip_first_distances: int = 0,
            set_equal_scale: bool | None = None,
            title: str | None = None,
            show_coordinates=False,
            interactive_mode: bool = False,
    ) -> None:
        """ Show 3D plot with 1 structure. """

        # Prepare to visualize
        fig: Figure = plt.figure()
        ax: Axes = fig.add_subplot(111, projection='3d')

        cls._plot_atoms_3d(
            fig=fig,
            ax=ax,
            coordinates=coordinates,
            visual_params=visual_params,
            to_build_bonds=to_build_bonds,
            num_of_min_distances=num_of_min_distances,
            skip_first_distances=skip_first_distances,
            set_equal_scale=set_equal_scale,
            show_coordinates=show_coordinates,
            interactive_mode=interactive_mode,
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
            visual_params_first: StructureVisualParams = VisualizationParams.carbon,
            visual_params_second: StructureVisualParams = VisualizationParams.al,
            title: str | None = None,
            show_coordinates: bool | None = None,
            show_indexes: bool | None = None,
            num_of_min_distances: int = 2,
            skip_first_distances: int = 0,
            interactive_mode: bool = False,
    ) -> None:
        """ Show 3D plot with 2 structures (by default there are carbon and aluminium) """

        # Prepare to visualize
        fig: Figure = plt.figure()
        ax: Axes = fig.add_subplot(111, projection='3d')

        # Plot first structure atoms (not interactive)
        cls._plot_atoms_3d(
            fig=fig,
            ax=ax,
            coordinates=coordinates_first,
            visual_params=visual_params_first,
            to_build_bonds=to_build_bonds,
            num_of_min_distances=num_of_min_distances,
            skip_first_distances=skip_first_distances,
            show_coordinates=show_coordinates,
            show_indexes=show_indexes,
            interactive_mode=False,
        )

        # Plot second structure atoms (interactive if enabled)
        cls._plot_atoms_3d(
            fig=fig,
            ax=ax,
            coordinates=coordinates_second,
            visual_params=visual_params_second,
            to_build_bonds=False,
            num_of_min_distances=1,
            skip_first_distances=0,
            show_coordinates=show_coordinates,
            show_indexes=show_indexes,
            interactive_mode=interactive_mode,
        )

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')  # type: ignore
        ax.legend()

        if title is not None:
            ax.set_title(title)

        plt.show()

    @classmethod
    def show_structures(
            cls,
            coordinates_list: list[np.ndarray],
            visual_params_list: list[StructureVisualParams],
            to_build_bonds_list: list[bool],
            title: str | None = None,
            show_coordinates: bool | None = None,
            show_indexes: bool | None = None,
            num_of_min_distances: int = 2,
            skip_first_distances: int = 0,
            interactive_mode: bool = False,
            custom_indices_list: list[list[int] | None] | None = None,
    ) -> None:
        """ Show 3D plot with multiple structures """

        if len(coordinates_list) != len(visual_params_list) != len(to_build_bonds_list):
            raise ValueError("len(coordinates_list) != len(visual_params_list) != len(to_build_bonds_list)")

        if custom_indices_list and len(custom_indices_list) != len(coordinates_list):
            raise ValueError("len(custom_indices_list) != len(coordinates_list)")

        # Prepare to visualize
        fig: Figure = plt.figure()
        ax: Axes = fig.add_subplot(111, projection='3d')

        for i, (coordinates, visual_params, to_build_bonds) in enumerate(zip(
                coordinates_list, visual_params_list, to_build_bonds_list)):
            custom_indices = custom_indices_list[i] if custom_indices_list else []

            # if i == 1:
            #     num_of_min_distances=3
            # elif i == 2:
            #     num_of_min_distances=3

            cls._plot_atoms_3d(
                fig=fig,
                ax=ax,
                coordinates=coordinates,
                visual_params=visual_params,
                to_build_bonds=to_build_bonds,
                num_of_min_distances=num_of_min_distances,
                skip_first_distances=skip_first_distances,
                show_coordinates=show_coordinates,
                show_indexes=show_indexes,
                interactive_mode=interactive_mode if visual_params.label != "Carbon" else False,
                custom_indexes=custom_indices if custom_indices else [],
            )

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')  # type: ignore
        ax.legend()

        if title is not None:
            ax.set_title(title)

        plt.show()

    @staticmethod
    def get_2d_plot(
        coordinates: np.ndarray,
        title: str | None = None,
        visual_params: StructureVisualParams = VisualizationParams.carbon,
        show_coordinates: bool | None = None,
        show_indexes: bool | None = None,
    ) -> Axes:
        # Prepare to visualize in 2D
        fig: Figure = plt.figure()
        ax: Axes = fig.add_subplot(111)  # No 3D projection here, just 2D

        # Plot points
        x: np.ndarray = coordinates[:, 0]
        y: np.ndarray = coordinates[:, 1]
        ax.scatter(
            x, y,
            color=visual_params.color_atoms,
            label='Points',
            alpha=visual_params.transparency,
        )

        if show_coordinates:
            # Show coordinates near each point
            for xx, yy in zip(x, y):
                ax.annotate(
                    f"({xx:.2f}, {yy:.2f})",
                    (xx, yy),
                    textcoords="offset points",
                    xytext=(5, 5),  # Offset from the point
                    fontsize=6,
                )

        if show_indexes:
            # Show coordinates near each point
            for xx, yy, i in zip(x, y, range(len(coordinates))):
                ax.annotate(
                    str(i),
                    (xx, yy),
                    textcoords="offset points",
                    xytext=(5, 5),  # Offset from the point
                    fontsize=6,
                )

        ax.set_xlabel('X')
        ax.set_ylabel('Y')

        if title is not None:
            ax.set_title(title)

        # ax.legend()
        plt.grid(True)

        return ax

    @classmethod
    def show_2d_graph(
        cls,
        coordinates: np.ndarray,
        title: str | None = None,
        visual_params: StructureVisualParams = VisualizationParams.carbon,
        show_coordinates: bool | None = None,
        show_indexes: bool | None = None,
    ) -> None:
        cls.get_2d_plot(
            coordinates,
            title,
            visual_params,
            show_coordinates,
            show_indexes
        )
        plt.show()

    @classmethod
    def _plot_atoms_3d(
            cls,
            fig: Figure,
            ax: Axes,
            coordinates: ndarray,
            visual_params: StructureVisualParams,
            set_equal_scale: bool | None = None,
            to_build_bonds: bool = True,
            num_of_min_distances: int = 2,
            skip_first_distances: int = 0,
            show_coordinates: bool | None = None,
            show_indexes: bool | None = None,
            interactive_mode: bool = False,
            custom_indexes: list[int] = [],
    ) -> PathCollection | None:
        if coordinates.size == 0:
            logger.warning(f"No points to plot for {visual_params.label}.")
            return

        x: ndarray = coordinates[:, 0]
        y: ndarray = coordinates[:, 1]
        z: ndarray = coordinates[:, 2]

        scatter: PathCollection = ax.scatter(
            x, y, z,
            c=visual_params.color_atoms,
            label=visual_params.label,
            s=visual_params.size,  # type: ignore
            alpha=visual_params.transparency,
            picker=True if interactive_mode else False,
        )

        if set_equal_scale is None:
            set_equal_scale = visual_params.set_equal_scale

        if set_equal_scale:
            cls._set_equal_scale(ax, x, y, z)

        if show_coordinates is True or (show_coordinates is None and visual_params.show_coordinates):
            # Show coordinates near each point
            for xx, yy, zz in zip(x, y, z):
                ax.text(
                    xx, yy, zz,
                    f"({xx:.2f}, {yy:.2f}, {zz:.2f})",  # type: ignore
                    fontsize=6,
                    color="black",
                    ha="center",
                    va="center",
                )

        if show_indexes is not False and visual_params.show_indexes:
            # Show coordinates near each point
            if len(custom_indexes) > 0:
                for i, (xx, yy, zz) in enumerate(coordinates):
                    ax.text(
                        xx, yy, zz,
                        str(custom_indexes[i]),  # type: ignore
                        fontsize=10,
                        color="black",
                        ha="center",
                        va="center",
                    )
            else:
                for i, (xx, yy, zz) in enumerate(coordinates):
                    ax.text(
                        xx, yy, zz,
                        str(i),  # type: ignore
                        fontsize=10,
                        color="black",
                        ha="center",
                        va="center",
                    )

        if to_build_bonds:
            # Carbon
            LinesBuilder.add_lines(
                coordinates=coordinates, ax=ax,
                visual_params=visual_params,
                num_of_min_distances=num_of_min_distances,
                skip_first_distances=skip_first_distances)

        if interactive_mode:
            def on_pick(event) -> None:
                try:
                    ind: int = event.ind[0]
                except Exception:
                    return

                while True:
                    try:
                        result: str = input(
                            f"Update coordinates for point with coordinates={coordinates[ind]}? (y/n): ")
                        if result == "y":
                            break
                        else:
                            return
                    except Exception:
                        return

                new_x = float(input(f"Enter new X ({coordinates[ind][0]:.2f}): ") or coordinates[ind][0])
                new_y = float(input(f"Enter new Y ({coordinates[ind][1]:.2f}): ") or coordinates[ind][1])
                new_z = float(input(f"Enter new Z ({coordinates[ind][2]:.2f}): ") or coordinates[ind][2])

                upd_point: list[float] = [new_x, new_y, new_z]

                coordinates[ind] = upd_point
                scatter._offsets3d = (coordinates[:, 0], coordinates[:, 1], coordinates[:, 2])  # type: ignore
                plt.draw()

                logger.info(f"Updated point coordinates: {coordinates[ind]}")

                # min_dists: np.ndarray = DistanceMeasure.calculate_min_distances(coordinates, np.array([upd_point]))
                # min_dist: float = np.min(min_dists[min_dists > 0])
                # logger.info(f"Distance to the nearest atom: {min_dist}")

            fig.canvas.mpl_connect('pick_event', on_pick)

        return scatter

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

        try:
            ax.set_zlim(z_mid - delta, z_mid + delta)  # type: ignore
        except Exception:
            pass
