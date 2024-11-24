from pathlib import Path
from numpy import ndarray

from src.utils import Constants, PathBuilder, FileReader, FileWriter, Logger, Inputs
from src.structure_visualizer import StructureVisualizer
from src.data_preparation import StructureSettings, StructureSettingsManager
from src.projects.intercalation_and_sorption import IntercalatedChannelBuilder, AlAtomsFilter
from src.base_structure_classes import AlLatticeType

from .init_data_parsing import AppActionsInitDataParsing
from .show_init_data import AppActionsShowInitData

logger = Logger("Actions")


class AppActionsIntercalationAndSorption:
    @classmethod
    def show_al_in_one_channel_structure(cls, structure_folder: str, to_set: bool) -> None:
        """
        Build one channel model from result_data/{structure_folder}/ljout-from-init-dat.pdb atoms
        based on result_data/{structure_folder}/structure_settings.json channel limits,
        filled with translated Al structure
        """

        structure_settings: None | StructureSettings = StructureSettingsManager.get_structure_settings(
            structure_folder=structure_folder)

        # Collect data to process

        # Carbon
        coordinates_carbon: ndarray = IntercalatedChannelBuilder.build_carbon_coordinates(
            structure_folder=structure_folder, structure_settings=structure_settings)

        # Aluminium
        to_open_calculated_atoms: bool = Inputs.bool_input(
            to_set, default_value=True, text="To open previously calculated atoms (if they exist)",
            env_id="open_calculated_atoms")

        if to_open_calculated_atoms:
            # Try to load previously calculated points from the file
            folder_path: Path = PathBuilder.build_path_to_result_data_dir()
            try:
                processed_coordinates_al: ndarray = FileReader.read_dat_file(
                    structure_folder=structure_folder, folder_path=folder_path)

            except FileNotFoundError:
                logger.warning(f"Calculated Al points for {structure_folder} not found.")
                processed_coordinates_al: ndarray = cls._calculate_al_points(
                    to_set, structure_settings, coordinates_carbon)

        else:
            processed_coordinates_al: ndarray = cls._calculate_al_points(
                to_set, structure_settings, coordinates_carbon)

        logger.info("Number of carbon atoms:", len(coordinates_carbon))
        logger.info("Number of al atoms:", len(processed_coordinates_al))

        FileWriter.write_dat_file(processed_coordinates_al, structure_folder=structure_folder)
        to_build_bonds: bool = Inputs.bool_input(to_set, default_value=True, text="To build bonds between atoms")

        StructureVisualizer.show_two_structures(
            coordinates_first=coordinates_carbon,
            coordinates_second=processed_coordinates_al,
            to_build_bonds=to_build_bonds,
            title=structure_folder)

    @classmethod
    def show_filtered_al_one_channel_structure(cls, structure_folder: str, to_set: bool) -> None:
        """
        Build one channel model from result_data/{structure_folder}/ljout-from-init-dat.pdb atoms
        based on result_data/{structure_folder}/structure_settings.json channel limits,
        filled with Al structure
        """

        structure_settings: None | StructureSettings = StructureSettingsManager.get_structure_settings(
            structure_folder=structure_folder)

        if structure_settings is None:
            logger.error("Please, provide structure_settings.")
            return

        # Carbon
        coordinates_carbon: ndarray = IntercalatedChannelBuilder.build_carbon_coordinates(
            structure_folder=structure_folder, structure_settings=structure_settings)

        # Aluminium
        coordinates_al: ndarray = cls._build_al_atoms(to_set, structure_settings, coordinates_carbon)
        coordinates_al = AlAtomsFilter.filter_al_atoms_related_carbon(
            coordinates_al, coordinates_carbon, structure_settings)

        to_build_bonds: bool = Inputs.bool_input(to_set, default_value=True, text="To build bonds between atoms")
        StructureVisualizer.show_two_structures(
            coordinates_first=coordinates_carbon,
            coordinates_second=coordinates_al,
            to_build_bonds=to_build_bonds,
            title=structure_folder)

    @staticmethod
    def _build_al_atoms(
            to_set: bool, structure_settings: None | StructureSettings, coordinates_carbon: ndarray
    ) -> ndarray:
        to_translate_al: bool = Inputs.bool_input(
            to_set, default_value=True, text="To translate AL atomes to fill full volume")

        al_lattice_type_str: str = Inputs.text_input(
            to_set,
            default_value="FCC",
            text=AlLatticeType.get_info(),
            available_values=AlLatticeType.get_available_types(),
            env_id="al_lattice_type")
        al_lattice_type = AlLatticeType(al_lattice_type_str)

        if al_lattice_type.is_cell:
            al_file: str = Inputs.text_input(to_set, default_value=Constants.filenames.AL_FILE, text="Init AL file")
            return IntercalatedChannelBuilder.build_al_coordinates_for_cell(
                to_translate_al=to_translate_al,
                al_file=al_file,
                structure_settings=structure_settings)

        else:
            # Fill the volume with aluminium for close-packed lattice
            return IntercalatedChannelBuilder.build_al_coordinates_for_close_packed(
                al_lattice_type=al_lattice_type,
                structure_settings=structure_settings,
                to_translate_al=to_translate_al)

    @classmethod
    def _calculate_al_points(
            cls, to_set: bool, structure_settings: None | StructureSettings, coordinates_carbon: ndarray
    ) -> ndarray:
        """ Calculate Al points inside channel from coordinates_carbon. """

        coordinates_al: ndarray = cls._build_al_atoms(to_set, structure_settings, coordinates_carbon)

        # Process data

        to_filter_al_atoms: bool = Inputs.bool_input(
            to_set, default_value=True, text="To filter AL atomes relative honeycomd bondaries")

        equidistant_al_points: bool = Inputs.bool_input(
            to_set=to_set, default_value=True, text="Set Al atoms maximally equidistant from the channel atoms",
            env_id="set_equidistant")

        return IntercalatedChannelBuilder.build_al_in_carbon(
            coordinates_carbon=coordinates_carbon,
            coordinates_al=coordinates_al,
            structure_settings=structure_settings,
            to_filter_al_atoms=to_filter_al_atoms,
            equidistant_al_points=equidistant_al_points)

    @classmethod
    def full_flow(cls, structure_folder: str, to_set: bool) -> None:
        """ Run all actions to intercalate (sorp) Al into carbon channel. """

        AppActionsInitDataParsing.convert_init_dat_to_pdb(structure_folder, to_set)
        AppActionsShowInitData.show_init_structure(structure_folder, to_set)
        AppActionsShowInitData.show_one_channel_structure(structure_folder, to_set)
        cls.show_filtered_al_one_channel_structure(structure_folder, to_set)
        cls.show_al_in_one_channel_structure(structure_folder, to_set)
