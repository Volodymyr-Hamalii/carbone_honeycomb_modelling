from numpy import ndarray

from ..utils import PathBuilder, FileReader, Logger
from ..structure_visualizer import AtomsUniverseBuilder
from ..coordinates_actions import StructureTranslator, CoordinateLimits, PlanesBuilder, CoordinatesFilter


logger = Logger(__name__)


class IntercalatedChannelBuilder:
    @classmethod
    def build_al_in_carbone(cls, structure_folder: str) -> tuple[ndarray, ndarray]:
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

        distance_from_plane: float = structure_settings["distance_from_plane"]  # type: ignore
        coordinates_al_filtered: ndarray = cls._filter_atoms_related_clannel_planes(
            coordinates=coordinates_al_translated,
            points_to_set_channel_planes=structure_settings["points_to_set_channel_planes"],
            distance_from_plane=distance_from_plane,
        )

        return coordinates_carbone, coordinates_al_filtered

    @staticmethod
    def _filter_atoms_related_clannel_planes(
            coordinates: ndarray,
            points_to_set_channel_planes: list[dict],
            distance_from_plane: float = 0,
    ) -> ndarray:

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
