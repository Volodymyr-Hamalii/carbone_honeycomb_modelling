from pathlib import Path
import numpy as np

from src.utils import Constants, PathBuilder, FileReader, FileWriter, Logger, Inputs
from src.base_structure_classes import AlLatticeType, Points, CoordinateLimits
from src.structure_visualizer import StructureVisualizer
from src.data_preparation import StructureSettings, StructureSettingsManager
from src.projects import (
    IntercalatedChannelBuilder,
    # AlAtomsFilter,
    CarbonHoneycombChannel,
    CarbonHoneycombActions,
    CoordinatesTableManager,
    AtomsParser,
    AlAtomsTranslator,
)

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

        structure_settings: StructureSettings = StructureSettingsManager.get_structure_settings(
            structure_folder=structure_folder)

        # Collect data to process

        # Carbon
        coordinates_carbon: Points = IntercalatedChannelBuilder.build_carbon_coordinates(
            structure_folder=structure_folder)

        carbon_channels: list[CarbonHoneycombChannel] = CarbonHoneycombActions.split_init_structure_into_separate_channels(
            coordinates_carbon=coordinates_carbon)
        carbon_channel: CarbonHoneycombChannel = carbon_channels[0]

        # Aluminium
        to_open_calculated_atoms: bool = Inputs.bool_input(
            to_set, default_value=True, text="To open previously calculated atoms (if they exist)",
            env_id="open_calculated_atoms")

        if to_open_calculated_atoms:
            # Try to load previously calculated points from the file
            folder_path: Path = PathBuilder.build_path_to_result_data_dir()
            try:
                al_points_data: np.ndarray = FileReader.read_dat_file(
                    structure_folder=structure_folder, folder_path=folder_path)

                al_points = Points(points=al_points_data)

            except FileNotFoundError:
                logger.warning(f"Calculated Al points for {structure_folder} not found.")
                al_points: Points = cls._calculate_al_points(
                    to_set, structure_settings, carbon_channel)

        else:
            al_points: Points = cls._calculate_al_points(
                to_set, structure_settings, carbon_channel)

        logger.info("Number of carbon atoms:", len(carbon_channel))
        logger.info("Number of al atoms:", len(al_points))

        FileWriter.write_dat_file(al_points.points, structure_folder=structure_folder)
        to_build_bonds: bool = Inputs.bool_input(
            to_set,
            default_value=True,
            text="To build bonds between atoms",
            env_id="to_build_bonds",
        )

        StructureVisualizer.show_two_structures(
            coordinates_first=carbon_channel.planes[0].points,
            coordinates_second=al_points.points,
            to_build_bonds=to_build_bonds,
            title=structure_folder,
            show_coordinates=True,
        )

    @staticmethod
    def fill_all_channels(structure_folder: str, to_set: bool) -> None:
        """
        Reads .dat files with carbon and built Al in one channel coordinates,
        translates these Al coordinates through all channels, shows and writes results.
        """

        folder_path: Path = PathBuilder.build_path_to_result_data_dir()

        al_points: np.ndarray = FileReader.read_dat_file(
            structure_folder=structure_folder,
            folder_path=folder_path,
        )

        coordinates_carbon: Points = IntercalatedChannelBuilder.build_carbon_coordinates(
            structure_folder=structure_folder)

        carbon_channels: list[CarbonHoneycombChannel] = CarbonHoneycombActions.split_init_structure_into_separate_channels(
            coordinates_carbon=coordinates_carbon)

        upd_al_points: np.ndarray = IntercalatedChannelBuilder.translate_al_points_through_channels(
            al_points, carbon_channels)

        to_build_bonds: bool = Inputs.bool_input(
            to_set,
            default_value=True,
            text="To build bonds between atoms",
            env_id="to_build_bonds",
        )

        StructureVisualizer.show_two_structures(
            coordinates_first=coordinates_carbon.points,
            coordinates_second=upd_al_points,
            to_build_bonds=to_build_bonds,
            title=structure_folder)

        FileWriter.write_dat_file(upd_al_points, structure_folder=structure_folder, filename="al_in_all_channels.dat")

    @staticmethod
    def update_al_coordinates_tbl(structure_folder: str, to_set: bool) -> None:
        """ 
        Run in loop:
        1) read Excel file to get Al atoms;
        2) filter empty lines (if some lines were removed);
        3) show current structure based on the coordinates from the file;
        4) after closing the plot go to 1).
        """
        carbon_channel: CarbonHoneycombChannel = AtomsParser.build_carbon_channel(structure_folder)
        to_build_bonds: bool = True

        number_of_planes: int = int(Inputs.text_input(
            to_set,
            default_value="1",
            text="Number of planes to translate",
            env_id="number_of_planes",
        ))

        while True:
            al_plane_coordinates: Points = AtomsParser.get_al_plane_coordinates(
                structure_folder, carbon_channel, number_of_planes)

            plane_points: np.ndarray = np.vstack([carbon_channel.planes[i].points for i in range(number_of_planes)])

            StructureVisualizer.show_two_structures(
                coordinates_first=plane_points,
                coordinates_second=al_plane_coordinates.points,
                to_build_bonds=to_build_bonds,
                title=structure_folder,
                # show_coordinates=False,
                # show_indexes=True,
            )

            try:
                CoordinatesTableManager.update_tbl_file(structure_folder, carbon_channel, number_of_planes)

            except IOError as e:
                logger.warning(e)
                file_is_closed: bool = Inputs.bool_input(
                    to_set,
                    default_value=True,
                    text="Did you save and close Excel file?",
                )

                if file_is_closed is False:
                    raise

    @staticmethod
    def translate_al_to_other_planes(structure_folder: str, to_set: bool) -> None:
        """ 
        Read Al coordinates from the Excel table and translate the structure to other planes.
        """
        carbon_channel: CarbonHoneycombChannel = AtomsParser.build_carbon_channel(structure_folder)
        to_build_bonds: bool = True

        al_plane_coordinates: Points = AtomsParser.get_al_channel_coordinates(structure_folder, carbon_channel)

        al_coordinates: Points = AlAtomsTranslator.translate_for_all_planes(carbon_channel, al_plane_coordinates)

        StructureVisualizer.show_two_structures(
            coordinates_first=carbon_channel.points,
            coordinates_second=al_coordinates.points,
            to_build_bonds=to_build_bonds,
            title=structure_folder,
            # show_coordinates=False,
            # show_indexes=True,
        )

        FileWriter.write_excel_file(
            df=al_coordinates.to_df(),
            structure_folder=structure_folder,
            sheet_name="Al atoms for the channel",
            file_name=Constants.filenames.AL_CHANNEL_COORDINATES_XLSX_FILE,
            is_init_data_dir=False,
        )

    @staticmethod
    def translate_al_to_all_channels(structure_folder: str, to_set: bool) -> None:
        """ 
        Read Al plane coordinates from the Excel table and translate the structure to other planes.
        """
        coordinates_carbon: Points = IntercalatedChannelBuilder.build_carbon_coordinates(
            structure_folder=structure_folder)

        carbon_channels: list[CarbonHoneycombChannel] = CarbonHoneycombActions.split_init_structure_into_separate_channels(
            coordinates_carbon=coordinates_carbon)

        al_channel_coordinates: Points = AtomsParser.get_al_channel_coordinates(
            structure_folder, carbon_channels.pop(0))
        al_coordinates: Points = AlAtomsTranslator.translate_for_all_channels(
            coordinates_carbon=coordinates_carbon,
            carbon_channels=carbon_channels,
            al_channel_coordinates=al_channel_coordinates,
        )

        to_build_bonds = True

        StructureVisualizer.show_two_structures(
            coordinates_first=coordinates_carbon.points,
            coordinates_second=al_coordinates.points,
            to_build_bonds=to_build_bonds,
            title=structure_folder,
            # show_coordinates=False,
            show_indexes=False,
        )

        FileWriter.write_excel_file(
            df=al_coordinates.to_df(columns=["i", "x_Al", "y_Al", "z_Al"]),
            structure_folder=structure_folder,
            sheet_name="Al atoms for the channel",
            file_name=Constants.filenames.AL_CHANNEL_COORDINATES_XLSX_FILE,
            is_init_data_dir=False,
        )

    @staticmethod
    def _build_al_atoms(
            to_set: bool, coordinate_limits: CoordinateLimits
    ) -> Points:
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
                coordinate_limits=coordinate_limits,
            )

        else:
            # Fill the volume with aluminium for close-packed lattice
            return IntercalatedChannelBuilder.build_al_coordinates_for_close_packed(
                al_lattice_type=al_lattice_type,
                coordinate_limits=coordinate_limits,
            )

    @classmethod
    def _calculate_al_points(
            cls, to_set: bool, structure_settings: StructureSettings, carbon_channel: CarbonHoneycombChannel
    ) -> Points:
        """ Calculate Al points inside channel from coordinates_carbon. """

        coordinates_al: Points = cls._build_al_atoms(to_set, carbon_channel.coordinate_limits)

        to_filter_al_atoms: bool = Inputs.bool_input(
            to_set, default_value=True, text="To filter AL atomes relative honeycomd bondaries")

        equidistant_al_points: bool = Inputs.bool_input(
            to_set=to_set, default_value=True, text="Set Al atoms maximally equidistant from the channel atoms",
            env_id="set_equidistant")

        return IntercalatedChannelBuilder.build_al_in_carbon(
            carbon_channel=carbon_channel,
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
        # cls.show_filtered_al_one_channel_structure(structure_folder, to_set)
        cls.show_al_in_one_channel_structure(structure_folder, to_set)
