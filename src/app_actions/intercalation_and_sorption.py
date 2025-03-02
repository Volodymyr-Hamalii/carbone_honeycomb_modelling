from pathlib import Path
import numpy as np
import pandas as pd

from src.utils import Constants, PathBuilder, FileReader, FileWriter, Logger, Inputs
from src.base_structure_classes import AlLatticeType, Points, CoordinateLimits
from src.structure_visualizer import StructureVisualizer, VisualizationParams
from src.data_preparation import StructureSettings
from src.coordinate_operations import DistanceMeasure
from src.projects import (
    IntercalatedChannelBuilder,
    # AlAtomsFilter,
    CarbonHoneycombChannel,
    CarbonHoneycombActions,
    CoordinatesTableManager,
    AtomsParser,
    AlAtomsTranslator,
    FullChannelBuilder,
)

from .init_data_parsing import AppActionsInitDataParsing
from .show_init_data import AppActionsShowInitData

logger = Logger("Actions")


class AppActionsIntercalationAndSorption:
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

        num_of_min_distances: int = int(Inputs.text_input(
            to_set,
            default_value="2",
            text="Number of min distances for bonds to show on plot",
            env_id="number_of_min_distances",
        ))

        StructureVisualizer.show_two_structures(
            coordinates_first=coordinates_carbon.points,
            coordinates_second=upd_al_points,
            to_build_bonds=to_build_bonds,
            title=structure_folder,
            num_of_min_distances=num_of_min_distances,
        )

        FileWriter.write_dat_file(upd_al_points, structure_folder=structure_folder, filename="al_in_all_channels.dat")

    @classmethod
    def update_al_coordinates_tbl(cls, structure_folder: str, to_set: bool) -> None:
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

        num_of_min_distances: int = int(Inputs.text_input(
            to_set,
            default_value="2",
            text="Number of min distances for bonds to show on plot",
            env_id="number_of_min_distances",
        ))

        show_al_layers: bool = Inputs.bool_input(
            to_set,
            default_value=False,
            text="Show AL layers in different colors",
            env_id="show_al_layers",
        )

        while True:
            al_plane_coordinates: Points = AtomsParser.get_al_plane_coordinates(
                structure_folder, carbon_channel, number_of_planes)

            plane_points: np.ndarray = np.vstack([carbon_channel.planes[i].points for i in range(number_of_planes)])

            cls._show_structures(
                carbon_channel_points=plane_points,
                al_points=al_plane_coordinates.points,
                to_build_bonds=to_build_bonds,
                title=structure_folder,
                num_of_min_distances=num_of_min_distances,
                show_al_layers=show_al_layers,
                # skip_first_distances=0,
                # show_coordinates=False,
                # show_indexes=True,
            )

            try:
                CoordinatesTableManager.update_plane_tbl_file(structure_folder, carbon_channel, number_of_planes)

            except IOError as e:
                logger.warning(e)
                file_is_closed: bool = Inputs.bool_input(
                    to_set,
                    default_value=True,
                    text="Did you save and close Excel file?",
                )

                if file_is_closed is False:
                    raise

    @classmethod
    def translate_al_to_other_planes(cls, structure_folder: str, to_set: bool) -> None:
        """
        Read Al coordinates from the Excel table and translate the structure to other planes.
        """
        carbon_channel: CarbonHoneycombChannel = AtomsParser.build_carbon_channel(structure_folder)
        to_build_bonds: bool = True

        number_of_planes: int = int(Inputs.text_input(
            to_set,
            default_value="1",
            text="Number of planes to translate",
            env_id="number_of_planes",
        ))

        num_of_min_distances: int = int(Inputs.text_input(
            to_set,
            default_value="2",
            text="Number of min distances for bonds to show on plot",
            env_id="number_of_min_distances",
        ))

        try_to_reflect_al_atoms: bool = Inputs.bool_input(
            to_set,
            default_value=True,
            text="Try to reflect aluminum atoms for other planes",
            env_id="try_to_reflect_al_atoms",
        )

        show_al_layers: bool = Inputs.bool_input(
            to_set,
            default_value=False,
            text="Show AL layers in different colors",
            env_id="show_al_layers",
        )

        al_coordinates: Points = AtomsParser.get_al_channel_coordinates(
            structure_folder, carbon_channel, number_of_planes, try_to_reflect_al_atoms)

        cls._show_structures(
            carbon_channel_points=carbon_channel.points,
            al_points=al_coordinates.points,
            to_build_bonds=to_build_bonds,
            title=structure_folder,
            num_of_min_distances=num_of_min_distances,
            show_al_layers=show_al_layers,
            # show_coordinates=False,
            # show_indexes=True,
        )

        FileWriter.write_excel_file(
            df=al_coordinates.to_df(columns=["i", "x_Al", "y_Al", "z_Al"]),
            structure_folder=structure_folder,
            sheet_name="Al atoms for the channel",
            file_name=Constants.filenames.AL_CHANNEL_COORDINATES_XLSX_FILE,
            is_init_data_dir=False,
        )

    @classmethod
    def update_al_full_channel_coordinates_tbl(cls, structure_folder: str, to_set: bool) -> None:
        """
        Run in loop:
        1) read Excel file to get Al atoms for the full channel;
        2) filter empty lines (if some lines were removed);
        3) show current structure based on the coordinates from the file;
        4) after closing the plot go to 1).
        """
        carbon_channel: CarbonHoneycombChannel = AtomsParser.build_carbon_channel(structure_folder)
        to_build_bonds: bool = True

        num_of_min_distances: int = int(Inputs.text_input(
            to_set,
            default_value="2",
            text="Number of min distances for bonds to show on plot",
            env_id="number_of_min_distances",
        ))

        show_al_layers: bool = Inputs.bool_input(
            to_set,
            default_value=False,
            text="Show AL layers in different colors",
            env_id="show_al_layers",
        )

        while True:
            al_full_channel_coordinates_df: pd.DataFrame | None = FileReader.read_excel_file(
                structure_folder=structure_folder,
                file_name=Constants.filenames.AL_FULL_CHANNEL_COORDINATES_XLSX_FILE,
                is_init_data_dir=False,
            )

            if al_full_channel_coordinates_df is not None:
                al_full_channel_coordinates: Points = AtomsParser._parse_al_coordinates_df(
                    al_full_channel_coordinates_df)
            else:
                logger.warning(
                    f"Excel file with Al atoms for the full channel not found in {structure_folder}; building full channel.")

                al_bulk: Points = cls._build_al_atoms(
                    to_set, carbon_channel.coordinate_limits)

                # number_of_planes: int = int(Inputs.text_input(
                #     to_set,
                #     default_value="1",
                #     text="Number of planes to translate",
                #     env_id="number_of_planes",
                # ))

                try_to_reflect_al_atoms: bool = Inputs.bool_input(
                    to_set,
                    default_value=True,
                    text="Try to reflect aluminum atoms for other planes",
                    env_id="try_to_reflect_al_atoms",
                )

                al_channel_planes_coordinates: Points = AtomsParser.get_al_channel_coordinates(
                    structure_folder, carbon_channel, number_of_planes=1,
                    try_to_reflect_al_atoms=try_to_reflect_al_atoms)

                al_full_channel_coordinates: Points = FullChannelBuilder.build_full_channel(
                    carbon_channel, al_channel_planes_coordinates, al_bulk)

                FileWriter.write_excel_file(
                    df=al_full_channel_coordinates.to_df(columns=["i", "x_Al", "y_Al", "z_Al"]),
                    structure_folder=structure_folder,
                    sheet_name="Al atoms for the full channel",
                    file_name=Constants.filenames.AL_FULL_CHANNEL_COORDINATES_XLSX_FILE,
                    is_init_data_dir=False,
                )

            cls._show_structures(
                carbon_channel_points=carbon_channel.points,
                al_points=al_full_channel_coordinates.points,
                to_build_bonds=to_build_bonds,
                title=structure_folder,
                num_of_min_distances=num_of_min_distances,
                show_al_layers=show_al_layers,
                # show_coordinates=False,
                # show_indexes=True,
            )

            try:
                CoordinatesTableManager.update_full_channel_tbl_file(structure_folder)

            except IOError as e:
                logger.warning(e)
                file_is_closed: bool = Inputs.bool_input(
                    to_set,
                    default_value=True,
                    text="Did you save and close Excel file?",
                )

                if file_is_closed is False:
                    raise

    @classmethod
    def translate_al_to_all_channels(cls, structure_folder: str, to_set: bool) -> None:
        """
        Read Al plane coordinates from the Excel table and translate the structure to other planes.
        """
        coordinates_carbon: Points = IntercalatedChannelBuilder.build_carbon_coordinates(
            structure_folder=structure_folder)

        carbon_channels: list[CarbonHoneycombChannel] = CarbonHoneycombActions.split_init_structure_into_separate_channels(
            coordinates_carbon=coordinates_carbon)

        number_of_planes: int = int(Inputs.text_input(
            to_set,
            default_value="1",
            text="Number of planes to translate",
            env_id="number_of_planes",
        ))

        num_of_min_distances: int = int(Inputs.text_input(
            to_set,
            default_value="2",
            text="Number of min distances for bonds to show on plot",
            env_id="number_of_min_distances",
        ))

        try_to_reflect_al_atoms: bool = Inputs.bool_input(
            to_set,
            default_value=True,
            text="Try to reflect aluminum atoms for other planes",
            env_id="try_to_reflect_al_atoms",
        )

        show_al_layers: bool = Inputs.bool_input(
            to_set,
            default_value=False,
            text="Show AL layers in different colors",
            env_id="show_al_layers",
        )

        al_channel_coordinates: Points = AtomsParser.get_al_channel_coordinates(
            structure_folder, carbon_channels.pop(0), number_of_planes, try_to_reflect_al_atoms)

        al_coordinates: Points = AlAtomsTranslator.translate_for_all_channels(
            coordinates_carbon=coordinates_carbon,
            carbon_channels=carbon_channels,
            al_channel_coordinates=al_channel_coordinates,
        )

        to_build_bonds = True

        cls._show_structures(
            carbon_channel_points=coordinates_carbon.points,
            al_points=al_coordinates.points,
            to_build_bonds=to_build_bonds,
            title=structure_folder,
            num_of_min_distances=num_of_min_distances,
            show_al_layers=show_al_layers,
            # show_coordinates=False,
            # show_indexes=False,
        )

        FileWriter.write_excel_file(
            df=al_coordinates.to_df(columns=["i", "x_Al", "y_Al", "z_Al"]),
            structure_folder=structure_folder,
            sheet_name="Al atoms for the channel",
            file_name=Constants.filenames.AL_ALL_CHANNELS_COORDINATES_XLSX_FILE,
            is_init_data_dir=False,
        )

        FileWriter.write_dat_file(
            al_coordinates.points,
            structure_folder=structure_folder,
            filename=Constants.filenames.AL_ALL_CHANNELS_COORDINATES_DAT_FILE,
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

    @staticmethod
    def _show_structures(
        carbon_channel_points: np.ndarray,
        al_points: np.ndarray,
        to_build_bonds: bool = False,
        title: str | None = None,
        show_coordinates: bool | None = None,
        show_indexes: bool | None = None,
        show_al_layers: bool | None = None,
        num_of_min_distances: int = 2,
        skip_first_distances: int = 0,
    ) -> None:

        if show_al_layers:
            # Split al points into layers by the min dists to Al
            min_dists: np.ndarray = DistanceMeasure.calculate_min_distances(al_points, carbon_channel_points)

            min_dist = np.min(min_dists)
            min_dist_with_threshold = min_dist * 1.2

            first_layer_points: np.ndarray = al_points[min_dists <= min_dist_with_threshold]
            other_layers_points: np.ndarray = al_points[min_dists > min_dist_with_threshold]

            StructureVisualizer.show_structures(
                coordinates_list=[
                    carbon_channel_points,
                    first_layer_points,
                    other_layers_points,
                ],
                visual_params_list=[
                    VisualizationParams.carbon,
                    VisualizationParams.al,
                    VisualizationParams.al_2,
                ],
                to_build_bonds_list=[
                    to_build_bonds,
                    False,
                    False,
                ],
                title=title,
                num_of_min_distances=num_of_min_distances,
                skip_first_distances=skip_first_distances,
                show_coordinates=show_coordinates,
                show_indexes=show_indexes,
            )

        else:
            StructureVisualizer.show_two_structures(
                coordinates_first=carbon_channel_points,
                coordinates_second=al_points,
                to_build_bonds=to_build_bonds,
                title=title,
                num_of_min_distances=num_of_min_distances,
                skip_first_distances=skip_first_distances,
                show_coordinates=show_coordinates,
                show_indexes=show_indexes,
            )

    @classmethod
    def full_flow(cls, structure_folder: str, to_set: bool) -> None:
        """ Run all actions to intercalate (sorp) Al into carbon channel. """

        AppActionsInitDataParsing.convert_init_dat_to_pdb(structure_folder, to_set)
        AppActionsShowInitData.show_init_structure(structure_folder, to_set)
        AppActionsShowInitData.show_one_channel_structure(structure_folder, to_set)
        # cls.show_filtered_al_one_channel_structure(structure_folder, to_set)
