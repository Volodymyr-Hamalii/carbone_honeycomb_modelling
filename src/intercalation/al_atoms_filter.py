import math
import numpy as np
from numpy import ndarray
from scipy.spatial.distance import cdist

from ..utils import Logger
from ..coordinates_actions import PlanesBuilder, CoordinatesFilter, StructureRotator
from ..data_preparation import StructureSettings, ChannelPoints

from .al_lattice_type import AlLatticeType

logger = Logger(__name__)


class AlAtomsFilter:

    @classmethod
    def find_max_filtered_atoms(
            cls,
            coordinates_carbone: ndarray,
            coordinates_al: ndarray,
            structure_settings: StructureSettings,
    ) -> ndarray:
        """
        Parallel move and filter Al atoms related carbone atoms
        to find the option with the maximum Al atoms after filtering.
        """

        max_atoms: int = 0
        xy_variance: float = np.inf

        coordinates_al_with_max_atoms: ndarray = coordinates_al.copy()

        al_lattice_parameter: float = structure_settings.al_lattice_parameter

        range_to_move: ndarray = np.arange(0, al_lattice_parameter, 0.1 * al_lattice_parameter)
        angle_range_to_rotate: ndarray = np.arange(- math.pi / 4, math.pi / 2, math.pi / 16)

        for step_x in range_to_move:
            moved_x_coordinates_al: ndarray = coordinates_al.copy()
            moved_x_coordinates_al[:, 0] += step_x

            for step_y in range_to_move:
                moved_xy_coordinates_al: ndarray = moved_x_coordinates_al.copy()
                moved_xy_coordinates_al[:, 1] += step_y

                for angle_x in angle_range_to_rotate:
                    x_rotaded_coordinates_al: ndarray = StructureRotator.rotate_on_angle_related_center(
                        moved_xy_coordinates_al.copy(), angle_x=angle_x)

                    for angle_y in angle_range_to_rotate:
                        xy_rotaded_coordinates_al: ndarray = StructureRotator.rotate_on_angle_related_center(
                            x_rotaded_coordinates_al.copy(), angle_y=angle_y)

                        for angle_z in angle_range_to_rotate:
                            xyz_rotaded_coordinates_al: ndarray = StructureRotator.rotate_on_angle_related_center(
                                xy_rotaded_coordinates_al.copy(), angle_z=angle_z)

                            coordinates_al_filtered: ndarray = cls.filter_al_atoms_related_carbone(
                                xyz_rotaded_coordinates_al, coordinates_carbone, structure_settings)

                            num_of_atoms: int = len(coordinates_al_filtered)

                            if num_of_atoms > max_atoms:
                                max_atoms = num_of_atoms
                                logger.info("max_atoms", max_atoms)

                                coordinates_al_with_max_atoms = coordinates_al_filtered

                                # Reset xy_variance
                                xy_variance = np.inf

                            elif num_of_atoms == max_atoms:
                                # Calculate the xy_variance to check if these X and Y coordinates
                                # are more equidistant (to have a more equilibrium structure)

                                # TODO: check case for complex structures
                                current_xy_variance = float(
                                    np.var(coordinates_al_filtered[:, 0]) + np.var(coordinates_al_filtered[:, 1]))

                                if current_xy_variance < xy_variance:
                                    coordinates_al_with_max_atoms = coordinates_al_filtered
                                    xy_variance = current_xy_variance

        return coordinates_al_with_max_atoms

    @classmethod
    def filter_al_atoms_related_carbone(
            cls,
            coordinates_al: ndarray,
            coordinates_carbone: ndarray,
            structure_settings: StructureSettings
    ) -> ndarray:
        """Filter Al atoms related planes and then related carbone atoms."""

        coordinates_al_filtered: ndarray = cls._filter_atoms_related_clannel_planes(
            coordinates=coordinates_al,
            points_to_set_channel_planes=structure_settings.points_to_set_channel_planes)

        coordinates_al_filtered = cls._filter_atoms_relates_carbone_atoms(
            coordinates_al=coordinates_al_filtered,
            coordinates_carbone=coordinates_carbone,
            max_distance_to_carbone_atoms=structure_settings.al_lattice_parameter / 2)

        return CoordinatesFilter.filter_by_min_max_z(
            coordinates_to_filter=coordinates_al_filtered,
            coordinates_with_min_max_z=coordinates_carbone,
            move_align_z=True)

    @staticmethod
    def _filter_atoms_related_clannel_planes(
            coordinates: ndarray,
            points_to_set_channel_planes: list[ChannelPoints],
            distance_from_plane: float = 0,
    ) -> ndarray:
        """
        Filter points from coordinates array by planes
        (remove atoms that are outside channel and atoms inside that are closer than distance_from_plane param).
        """

        filtered_coordinates: ndarray = coordinates

        for plane_data in points_to_set_channel_planes:
            # Build plane parameters
            A, B, C, D = PlanesBuilder.build_plane_parameters(
                p1=plane_data.points[0],
                p2=plane_data.points[1],
                p3=plane_data.points[2])

            filtered_coordinates = CoordinatesFilter.filter_coordinates_related_to_plane(
                filtered_coordinates,
                A, B, C, D,
                min_distance=distance_from_plane,
                direction=plane_data.direction)

        return filtered_coordinates

    @staticmethod
    def _filter_atoms_relates_carbone_atoms(
        coordinates_al: ndarray,
        coordinates_carbone: ndarray,
        max_distance_to_carbone_atoms: float,
    ) -> ndarray:
        """
        Filter points from the coordinates_al array
        by the max distance (max_distance_to_carbone_atoms param)
        to the points in the coordinates_carbone array.
        """

        # Calculate the distance matrix
        distances_matrix: ndarray = cdist(coordinates_al, coordinates_carbone)

        # Find the minimum distance for each atom in coordinates_al to any atom in coordinates_carbone
        min_distances: ndarray = np.min(distances_matrix, axis=1)

        # Filter the atoms in coordinates_al based on the max distance to carbone atoms
        filtered_coordinates: ndarray = coordinates_al[min_distances >= max_distance_to_carbone_atoms]

        return filtered_coordinates
