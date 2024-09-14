import numpy as np
from numpy import ndarray
from scipy.spatial.distance import cdist

from ..utils import Logger
from ..coordinates_actions import PlanesBuilder, CoordinatesFilter
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
        coordinates_al_with_max_atoms: ndarray = coordinates_al.copy()

        al_lattice_parameter: float = structure_settings.al_lattice_parameter

        logger.info(structure_settings)

        range_to_move: ndarray = np.arange(0, al_lattice_parameter, 0.05 * al_lattice_parameter)

        for step_x in range_to_move:
            moved_x_coordinates_al: ndarray = coordinates_al.copy()
            moved_x_coordinates_al[:, 0] += step_x

            for step_y in range_to_move:
                moved_xy_coordinates_al: ndarray = moved_x_coordinates_al.copy()
                moved_xy_coordinates_al[:, 1] += step_y

                coordinates_al_filtered: ndarray = cls.filter_al_atoms_related_carbone(
                    moved_xy_coordinates_al, coordinates_carbone, structure_settings)

                num_of_atoms: int = len(coordinates_al_filtered)

                if num_of_atoms > max_atoms:
                    logger.info("max_atoms", max_atoms)

                    max_atoms = num_of_atoms
                    coordinates_al_with_max_atoms = coordinates_al_filtered

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

        return cls._filter_atoms_relates_carbone_atoms(
            coordinates_al=coordinates_al_filtered,
            coordinates_carbone=coordinates_carbone,
            max_distance_to_carbone_atoms=structure_settings.al_lattice_parameter / 2)

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
