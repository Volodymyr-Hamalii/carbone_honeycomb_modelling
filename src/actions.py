from numpy import ndarray

from .utils import PathBuilder, FilesConverter, Logger
from .structure_visualizer import StructureVisualizer, AtomsUniverseBuilder, VisualizationParameters
from .data_preparation import StructureSettings, StructureSettingsManager, ChannelLimits
from .intercalation import IntercalatedChannelBuilder

logger = Logger(__name__)


class Actions:
    @classmethod
    def help(cls, structure_folder: str) -> None:
        """ Print all available methods. """

        # Get all attributes of the class
        attributes: list[str] = dir(cls)

        # Filter out only methods
        methods: list[str] = [attr for attr in attributes if callable(getattr(cls, attr)) and not attr.startswith("__")]

        print()
        logger.info("Available actions:")
        for method in methods:
            print(method)

        print()
        logger.info("To run specific action provide also structure_folder as a parameter. For example:")
        print(f"python3 main.py show_init_structure {structure_folder}\n")

    @staticmethod
    def convert_init_dat_to_pdb(structure_folder: str) -> None:
        """
        Convert init_data/{structure_folder}/ljout.dat into result_data/{structure_folder}/ljout-result.pdb
        Also create result_data/{structure_folder}/structure_settings.json template if it didn't exist.
        """

        FilesConverter.dat_to_pdb(structure_folder=structure_folder)

        # Create template for coordinates if it doesn't exists
        StructureSettingsManager.create_structure_settings_template(structure_folder=structure_folder)

    @staticmethod
    def show_init_structure(structure_folder: str) -> None:
        """ Show 3D model of result_data/{structure_folder}/ljout-result.pdb """

        path_to_init_pdb_file: str = PathBuilder.build_path_to_result_data_file(structure_folder)
        coordinates: ndarray = AtomsUniverseBuilder.builds_atoms_coordinates(path_to_init_pdb_file)
        StructureVisualizer.show_structure(coordinates, to_build_bonds=True)

    @staticmethod
    def show_init_al_structure() -> None:
        """ Show 3D model of init_data/al.pdb """

        path_to_al_pdb_file: str = PathBuilder.build_path_to_init_data_file(file="al.pdb")
        coordinates: ndarray = AtomsUniverseBuilder.builds_atoms_coordinates(path_to_al_pdb_file)

        StructureVisualizer.show_structure(
            coordinates=coordinates,
            to_build_bonds=True,
            visual_parameters=VisualizationParameters.al,
            num_of_min_distances=1,
            skip_first_distances=2)

    @staticmethod
    def show_one_channel_structure(structure_folder: str) -> None:
        """
        Build one channel model from result_data/{structure_folder}/ljout-result.pdb atoms
        based on result_data/{structure_folder}/structure_settings.json channel limits.

        Write result to result_data/{structure_folder}/ljout-result-one-channel.pdb if it didn't exist.
        """

        path_to_init_pdb_file: str = PathBuilder.build_path_to_result_data_file(structure_folder)

        structure_settings: None | StructureSettings = StructureSettingsManager.read_file(
            structure_folder=structure_folder)

        channel_limits: ChannelLimits | None = structure_settings.channel_limits if structure_settings else None

        coordinates: ndarray = AtomsUniverseBuilder.builds_atoms_coordinates(
            path_to_init_pdb_file, channel_limits)
        StructureVisualizer.show_structure(coordinates, to_build_bonds=True)

    @staticmethod
    def show_al_in_one_channel_structure(structure_folder: str) -> None:
        """
        Build one channel model from result_data/{structure_folder}/ljout-result.pdb atoms
        based on result_data/{structure_folder}/structure_settings.json channel limits,
        filled with translated Al structure from init_data/al.pdb
        """

        coordinates: tuple[ndarray, ndarray] = IntercalatedChannelBuilder.build_al_in_carbone(
            structure_folder=structure_folder, filter_al_atoms=True)

        coordinates_carbone: ndarray = coordinates[0]
        coordinates_al: ndarray = coordinates[1]

        logger.info("Number of carbone atoms:", len(coordinates_carbone))
        logger.info("Number of al atoms:", len(coordinates_al))

        StructureVisualizer.show_two_structures(
            coordinates_first=coordinates_carbone,
            coordinates_second=coordinates_al,
            to_build_bonds=False)

    @classmethod
    def full_flow(cls, structure_folder: str) -> None:
        """ Run all actions """

        cls.convert_init_dat_to_pdb(structure_folder)
        cls.show_init_structure(structure_folder)
        cls.show_init_al_structure()
        cls.show_one_channel_structure(structure_folder)
        cls.show_al_in_one_channel_structure(structure_folder)
