import numpy as np
from numpy import ndarray
from scipy.spatial.distance import cdist

from ..utils import PathBuilder, FileReader, Logger
from ..structure_visualizer import AtomsUniverseBuilder
from ..coordinates_actions import StructureTranslator, PlanesBuilder, CoordinatesFilter
from ..data_preparation import StructureSettingsManager, ChannelLimits, StructureSettings, ChannelPoints

from .al_lattice_type import AlLatticeType


logger = Logger(__name__)


class IntercalatedChannelBuilder:
    @staticmethod
    def build_carbone_coordinates(structure_folder: str, structure_settings: None | StructureSettings) -> ndarray:
        path_to_pdb_file: str = PathBuilder.build_path_to_result_data_file(structure_folder)

        return AtomsUniverseBuilder.builds_atoms_coordinates(
            path_to_pdb_file,
            channel_coordinate_limits=structure_settings.channel_limits if structure_settings else None)

    @staticmethod
    def build_al_coordinates_for_cell(
            to_translate_al: bool = True,
            al_file: str = "al.pdb",
            structure_settings: None | StructureSettings = None,
    ) -> ndarray:
        path_to_al_pdb_file: str = PathBuilder.build_path_to_init_data_file(file=al_file)
        coordinates_al: ndarray = AtomsUniverseBuilder.builds_atoms_coordinates(path_to_al_pdb_file)

        if to_translate_al:
            if structure_settings is None:
                raise ValueError(
                    "To translate Al please provide structure_settings.json file with channel_limits.")

            return StructureTranslator.translate_cell(
                cell_coordinates=coordinates_al, translation_limits=structure_settings.channel_limits)

        return coordinates_al

    @staticmethod
    def build_al_coordinates_for_close_packed(
            al_lattice_type: AlLatticeType,
            structure_settings: None | StructureSettings = None,
            to_translate_al: bool = True,
    ) -> ndarray:

        if structure_settings is None or structure_settings.al_lattice_parameter == 0:
            raise ValueError(
                "To translate Al please provide structure_settings.json file with channel_limits and al_lattice_parameter.")

        return AtomsUniverseBuilder.build_close_packed_structure(
            lattice_parameter=structure_settings.al_lattice_parameter,
            lattice_type=al_lattice_type,
            channel_coordinate_limits=structure_settings.channel_limits,
            to_translate_al=to_translate_al)

    @classmethod
    def build_al_in_carbone(
            cls,
            coordinates_carbone: ndarray,
            coordinates_al: ndarray,
            structure_settings: None | StructureSettings,
            to_filter_al_atoms: bool = True,
    ) -> tuple[ndarray, ndarray]:
        """ Return coordinates_carbone, coordinates_al """

        if to_filter_al_atoms:
            if structure_settings is None:
                logger.error("Not able to filter AL atoms without structure_settings.json")
            else:
                # distance_from_plane: float = structure_settings["distance_from_plane"]
                coordinates_al_filtered: ndarray = cls._filter_atoms_related_clannel_planes(
                    coordinates=coordinates_al,
                    points_to_set_channel_planes=structure_settings.points_to_set_channel_planes)

                coordinates_al_filtered: ndarray = cls._filter_atoms_relates_carbone_atoms(
                    coordinates_al=coordinates_al_filtered,
                    coordinates_carbone=coordinates_carbone,
                    max_distance_to_carbone_atoms=structure_settings.max_distance_to_carbone_atoms)

                return coordinates_carbone, coordinates_al_filtered

        return coordinates_carbone, coordinates_al

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
