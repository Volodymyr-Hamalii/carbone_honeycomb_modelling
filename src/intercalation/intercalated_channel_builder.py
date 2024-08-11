import numpy as np
from numpy import ndarray
from scipy.spatial.distance import cdist

from ..utils import PathBuilder, FileReader, Logger
from ..structure_visualizer import AtomsUniverseBuilder
from ..coordinates_actions import StructureTranslator, CoordinateLimits, PlanesBuilder, CoordinatesFilter


logger = Logger(__name__)


class IntercalatedChannelBuilder:
    @classmethod
    def build_al_in_carbone(cls, structure_folder: str, filter_al_atoms: bool = True) -> tuple[ndarray, ndarray]:
        """ Return coordinates_carbone, coordinates_al """
        path_to_init_pdb_file: str = PathBuilder.build_path_to_result_data_file(structure_folder)
        structure_settings: dict | None = FileReader.read_json_file(structure_folder)  # type: ignore

        if structure_settings is None:
            raise FileExistsError(f"Settings file for {structure_folder} not found")

        # Build one carbone channel
        channel_limits: dict[str, float] = structure_settings["channel_limits"]
        # TODO: to automate searching for these limits instead of providing them in the file
        channel_coordinate_limits = CoordinateLimits(
            x_min=channel_limits["x_min"], x_max=channel_limits["x_max"],
            y_min=channel_limits["y_min"], y_max=channel_limits["y_max"],
            z_min=channel_limits["z_min"], z_max=channel_limits["z_max"],
        )

        coordinates_carbone: ndarray = AtomsUniverseBuilder.builds_atoms_coordinates(
            path_to_init_pdb_file, channel_coordinate_limits)

        # Build AL structure

        path_to_al_pdb_file: str = PathBuilder.build_path_to_init_data_file(file="al.pdb")
        coordinates_al: ndarray = AtomsUniverseBuilder.builds_atoms_coordinates(path_to_al_pdb_file)

        coordinates_al_translated: ndarray = StructureTranslator.translate(
            coordinates=coordinates_al, translation_limits=channel_coordinate_limits
        )

        if filter_al_atoms:

            # distance_from_plane: float = structure_settings["distance_from_plane"]
            coordinates_al_filtered: ndarray = cls._filter_atoms_related_clannel_planes(
                coordinates=coordinates_al_translated,
                points_to_set_channel_planes=structure_settings["points_to_set_channel_planes"],
                # distance_from_plane=distance_from_plane,
            )

            max_distance_to_carbone_atoms: float = structure_settings["max_distance_to_carbone_atoms"]

            coordinates_al_filtered: ndarray = cls._filter_atoms_relates_carbone_atoms(
                coordinates_al=coordinates_al_filtered,
                coordinates_carbone=coordinates_carbone,
                max_distance_to_carbone_atoms=max_distance_to_carbone_atoms,
            )

            return coordinates_carbone, coordinates_al_filtered
        return coordinates_carbone, coordinates_al

    @staticmethod
    def _filter_atoms_related_clannel_planes(
            coordinates: ndarray,
            points_to_set_channel_planes: list[dict],
            distance_from_plane: float = 0,
    ) -> ndarray:
        """
        Filter points from coordinates array by planes
        (remove atoms that are outside channel and atoms inside that are closer than distance_from_plane param).
        """

        filtered_coordinates: ndarray = coordinates

        for plane_data in points_to_set_channel_planes:
            plane_points: list[list[float]] = plane_data["points"]

            # Build plane parameters
            A, B, C, D = PlanesBuilder.build_plane_parameters(
                p1=plane_points[0],
                p2=plane_points[1],
                p3=plane_points[2],
            )

            direction: int = plane_data["direction"]

            filtered_coordinates = CoordinatesFilter.filter_coordinates_related_to_plane(
                filtered_coordinates,
                A, B, C, D,
                min_distance=distance_from_plane,
                direction=direction,
            )

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
