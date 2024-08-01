from numpy import ndarray

from .utils import PathBuilder, FilesConverter, Logger
from .structure_visualizer import StructureVisualizer, AtomsUniverseBuilder
from .coordinates_actions import CoordinateLimits
from .intercalation import IntercalatedChannelBuilder

logger = Logger(__name__)


class Actions:
    @staticmethod
    def convert_init_dat_to_pdb(structure_folder: str) -> None:

        FilesConverter.dat_to_pdb(
            structure_folder=structure_folder,
        )

    @staticmethod
    def show_init_structure(structure_folder: str) -> None:
        path_to_init_pdb_file: str = PathBuilder.build_path_to_result_data_file(structure_folder)
        coordinates: ndarray = AtomsUniverseBuilder.builds_atoms_coordinates(path_to_init_pdb_file)
        StructureVisualizer.show_structure(coordinates)

    @staticmethod
    def show_one_channel_structure(structure_folder: str) -> None:
        path_to_init_pdb_file: str = PathBuilder.build_path_to_result_data_file(structure_folder)

        channel_coordinate_limits = CoordinateLimits(
            x_min=-2.00,
            x_max=4.50,
            y_min=-0.40,
            y_max=4.8,
        )

        coordinates: ndarray = AtomsUniverseBuilder.builds_atoms_coordinates(
            path_to_init_pdb_file, channel_coordinate_limits)
        StructureVisualizer.show_structure(coordinates, to_build_bonds=True)

    @staticmethod
    def show_al_in_one_channel_structure(structure_folder: str) -> None:
        coordinates: tuple[ndarray, ndarray] = IntercalatedChannelBuilder.build_al_in_carbone(
            structure_folder=structure_folder
        )
        coordinates_carbone: ndarray = coordinates[0]
        coordinates_al: ndarray = coordinates[1]

        logger.info("Number of carbone atoms:", len(coordinates_carbone))
        logger.info("Number of al atoms:", len(coordinates_al))

        StructureVisualizer.show_two_structures(
            coordinates_first=coordinates_carbone,
            coordinates_second=coordinates_al,
            coordinates_second_transparency=0.15,
            to_build_bonds=False)

    @classmethod
    def full_flow(cls, structure_folder: str) -> None:
        # cls.convert_init_dat_to_pdb(structure_folder)
        cls.show_init_structure(structure_folder)
        cls.show_one_channel_structure(structure_folder)
        cls.show_al_in_one_channel_structure(structure_folder)
