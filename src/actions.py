from numpy import ndarray

from .utils import PathBuilder, FilesConverter, Logger, Inputs
from .structure_visualizer import StructureVisualizer, AtomsUniverseBuilder, VisualizationParameters
from .data_preparation import StructureSettings, StructureSettingsManager, ChannelLimits
from .intercalation import IntercalatedChannelBuilder, AlLatticeType

logger = Logger(__name__)


class Actions:
    @classmethod
    def help(cls, structure_folder: str) -> None:
        """ Print all available methods and general info. """

        # Get all attributes of the class
        attributes: list[str] = dir(cls)

        # Filter out only methods
        methods: list[str] = [attr for attr in attributes if callable(getattr(cls, attr)) and not attr.startswith("__")]

        message: str = "Available actions:\n"
        for method in methods:
            message = f"{message}* {method}\n"

        message = message + "\nmain.py takes the following parameters: action, structure_folder, set.\n" + \
            "   1) action -- any action from the list above,\n" + \
            "   2) structure_folder -- structure to process from 'init_data' or 'result_data' folders,\n" + \
            "   3) set -- parameters custumization (to build bonds, to filter atoms etc.); just put 'set' to customize building,\n" + \
            "   4) optional arguments -- specific parameters if some method needs them.\n" + \
            f"For example:\npython3 main.py show_init_structure {structure_folder} set\n" + \
            "If you don't want to specify some argument -- just provide '_' on this place."

        logger.info(message)

    @staticmethod
    def convert_init_dat_to_pdb(structure_folder: str, to_set: bool) -> None:
        """
        Convert init_data/{structure_folder}/ljout.dat into result_data/{structure_folder}/ljout-result.pdb
        Also create result_data/{structure_folder}/structure_settings.json template if it didn't exist.
        """

        FilesConverter.dat_to_pdb(structure_folder=structure_folder)

        # Create template for coordinates if it doesn't exists
        StructureSettingsManager.create_structure_settings_template(structure_folder=structure_folder)

    @staticmethod
    def show_init_structure(structure_folder: str, to_set: bool) -> None:
        """ Show 3D model of result_data/{structure_folder}/ljout-result.pdb """

        path_to_init_pdb_file: str = PathBuilder.build_path_to_result_data_file(structure_folder)
        coordinates: ndarray = AtomsUniverseBuilder.builds_atoms_coordinates(path_to_init_pdb_file)

        to_build_bonds: bool = Inputs.bool_input(to_set, default_value=True, text="To build bonds between atoms")
        StructureVisualizer.show_structure(coordinates, to_build_bonds=to_build_bonds)

    @staticmethod
    def show_init_al_structure(_, to_set: bool) -> None:
        """ Show 3D model of init_data/al.pdb """

        al_file: str = Inputs.text_input(to_set, default_value="al.pdb", text="Init AL file")
        path_to_al_pdb_file: str = PathBuilder.build_path_to_init_data_file(file=al_file)
        coordinates: ndarray = AtomsUniverseBuilder.builds_atoms_coordinates(path_to_al_pdb_file)

        to_build_bonds: bool = Inputs.bool_input(to_set, default_value=True, text="To build bonds between atoms")
        StructureVisualizer.show_structure(
            coordinates=coordinates,
            to_build_bonds=to_build_bonds,
            visual_parameters=VisualizationParameters.al,
            num_of_min_distances=1,
            skip_first_distances=2)

    @staticmethod
    def show_one_channel_structure(structure_folder: str, to_set: bool) -> None:
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

        to_build_bonds: bool = Inputs.bool_input(to_set, default_value=True, text="To build bonds between atoms")
        StructureVisualizer.show_structure(coordinates, to_build_bonds=to_build_bonds)

    @staticmethod
    def show_al_in_one_channel_structure(structure_folder: str, to_set: bool) -> None:
        """
        Build one channel model from result_data/{structure_folder}/ljout-result.pdb atoms
        based on result_data/{structure_folder}/structure_settings.json channel limits,
        filled with translated Al structure from init_data/al.pdb
        """

        structure_settings: None | StructureSettings = StructureSettingsManager.read_file(
            structure_folder=structure_folder)
        
        ### Collect data to process

        # Carbone
        coordinates_carbone: ndarray = IntercalatedChannelBuilder.build_carbone_coordinates(
            structure_folder=structure_folder, structure_settings=structure_settings)

        # Aluminium

        to_filter_al_atoms: bool = Inputs.bool_input(
            to_set, default_value=True, text="To filter AL atomes relative honeycomd bondaries")
        to_translate_al: bool = Inputs.bool_input(
            to_set, default_value=True, text="To translate AL atomes to fill full volume")

        al_lattice_type_str: str = Inputs.text_input(
            to_set, default_value="cell", text=AlLatticeType.get_info())
        al_lattice_type = AlLatticeType(al_lattice_type_str)

        if al_lattice_type.is_cell:
            al_file: str = Inputs.text_input(to_set, default_value="al.pdb", text="Init AL file")
            coordinates_al: ndarray = IntercalatedChannelBuilder.build_al_coordinates_for_cell(
                to_translate_al=to_translate_al,
                al_file=al_file,
                structure_settings=structure_settings)

        else:
            # Fill the volume with aluminium for close-packed lattice
            coordinates_al: ndarray = IntercalatedChannelBuilder.build_al_coordinates_for_close_packed(
                al_lattice_type=al_lattice_type,
                structure_settings=structure_settings,
                to_translate_al=to_translate_al)

        ### Process data

        coordinates: tuple[ndarray, ndarray] = IntercalatedChannelBuilder.build_al_in_carbone(
            coordinates_carbone=coordinates_carbone,
            coordinates_al=coordinates_al,
            structure_settings=structure_settings,
            to_filter_al_atoms=to_filter_al_atoms)

        processed_coordinates_carbone: ndarray = coordinates[0]
        processed_coordinates_al: ndarray = coordinates[1]

        logger.info("Number of carbone atoms:", len(processed_coordinates_carbone))
        logger.info("Number of al atoms:", len(processed_coordinates_al))

        to_build_bonds: bool = Inputs.bool_input(to_set, default_value=True, text="To build bonds between atoms")

        StructureVisualizer.show_two_structures(
            coordinates_first=processed_coordinates_carbone,
            coordinates_second=processed_coordinates_al,
            to_build_bonds=to_build_bonds)

    @classmethod
    def full_flow(cls, structure_folder: str, to_set: bool) -> None:
        """ Run all actions """

        cls.convert_init_dat_to_pdb(structure_folder, to_set)
        cls.show_init_structure(structure_folder, to_set)
        cls.show_init_al_structure(None, to_set)
        cls.show_one_channel_structure(structure_folder, to_set)
        cls.show_al_in_one_channel_structure(structure_folder, to_set)
