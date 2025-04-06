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
from .view_model_params_setter import VMParamsSetter


logger = Logger("Actions")


class VMIntercalationAndSorption(VMParamsSetter):
    # def fill_all_channels(self, structure_folder: str) -> None:
    #     """
    #     Reads .dat files with carbon and built Al in one channel coordinates,
    #     translates these Al coordinates through all channels, shows and writes results.
    #     """

    #     folder_path: Path = PathBuilder.build_path_to_result_data_dir()

    #     al_points: np.ndarray = FileReader.read_dat_file(
    #         structure_folder=structure_folder,
    #         folder_path=folder_path,
    #     )

    #     coordinates_carbon: Points = IntercalatedChannelBuilder.build_carbon_coordinates(
    #         structure_folder=structure_folder)

    #     carbon_channels: list[CarbonHoneycombChannel] = CarbonHoneycombActions.split_init_structure_into_separate_channels(
    #         coordinates_carbon=coordinates_carbon)

    #     upd_al_points: np.ndarray = IntercalatedChannelBuilder.translate_al_points_through_channels(
    #         al_points, carbon_channels)

    #     to_build_bonds: bool = Inputs.bool_input(
    #         to_set,
    #         default_value=True,
    #         text="To build bonds between atoms",
    #         env_id="to_build_bonds",
    #     )

    #     num_of_min_distances: int = int(Inputs.text_input(
    #         to_set,
    #         default_value="2",
    #         text="Number of min distances for bonds to show on plot",
    #         env_id="number_of_min_distances",
    #     ))

    #     StructureVisualizer.show_two_structures(
    #         coordinates_first=coordinates_carbon.points,
    #         coordinates_second=upd_al_points,
    #         to_build_bonds=to_build_bonds,
    #         title=structure_folder,
    #         num_of_min_distances=num_of_min_distances,
    #     )

    #     FileWriter.write_dat_file(upd_al_points, structure_folder=structure_folder, filename="al_in_all_channels.dat")

    def update_plane_tbl_excel_file(self, structure_folder: str) -> None:
        carbon_channel: CarbonHoneycombChannel = AtomsParser.build_carbon_channel(structure_folder)

        try:
            CoordinatesTableManager.update_plane_tbl_file(structure_folder, carbon_channel, self.number_of_planes)

        except IOError as e:
            logger.warning(e)
            raise IOError(
                f"Error updating plane table file for {structure_folder} (possible reason: Excel file is blocked)"
            ) from e

    def update_al_coordinates_tbl(self, structure_folder: str) -> None:
        """
        Run in loop:
        1) read Excel file to get Al atoms;
        2) filter empty lines (if some lines were removed);
        3) show current structure based on the coordinates from the file;
        4) after closing the plot go to 1).
        """
        carbon_channel: CarbonHoneycombChannel = AtomsParser.build_carbon_channel(structure_folder)

        al_plane_coordinates: Points = AtomsParser.get_al_plane_coordinates(
            structure_folder, carbon_channel, self.number_of_planes, self.file_name)

        plane_points: np.ndarray = np.vstack(
            [carbon_channel.planes[i].points for i in range(self.number_of_planes)])

        self._show_structures(
            carbon_channel_points=plane_points,
            al_points=al_plane_coordinates.points,
            title=structure_folder,
        )

    def translate_al_to_other_planes(self, structure_folder: str) -> Path:
        """
        Read Al coordinates from the Excel table and translate the structure to other planes.
        """
        carbon_channel: CarbonHoneycombChannel = AtomsParser.build_carbon_channel(structure_folder)

        al_coordinates: Points = AtomsParser.get_al_channel_coordinates(
            structure_folder, carbon_channel, self.number_of_planes, self.to_try_to_reflect_al_atoms)

        if self.number_of_planes > 1:
            # Build only specified planes
            carbon_channel_points: np.ndarray = np.vstack(
                [carbon_channel.planes[i].points for i in range(self.number_of_planes)])
        else:
            # Build all planes
            carbon_channel_points: np.ndarray = carbon_channel.points

        # min_dists: np.ndarray = DistanceMeasure.calculate_min_distances(al_coordinates.points, carbon_channel_points)
        # DistanceMeasure.calculate_min_distances_between_points()

        self._show_structures(
            carbon_channel_points=carbon_channel_points,
            al_points=al_coordinates.points,
            title=structure_folder,
        )

        path_to_file: Path | None = FileWriter.write_excel_file(
            df=al_coordinates.to_df(columns=["i", "x_Al", "y_Al", "z_Al"]),
            structure_folder=structure_folder,
            sheet_name="Al atoms for the channel",
            file_name=Constants.filenames.AL_CHANNEL_COORDINATES_XLSX_FILE,
            is_init_data_dir=False,
        )

        if path_to_file is None:
            raise ValueError("Failed to write Excel file")

        return path_to_file

    def update_al_full_channel_coordinates_tbl(self, structure_folder: str) -> None:
        """
        Run in loop:
        1) read Excel file to get Al atoms for the full channel;
        2) filter empty lines (if some lines were removed);
        3) show current structure based on the coordinates from the file;
        4) after closing the plot go to 1).
        """
        return
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

        interactive_mode: bool = Inputs.bool_input(
            to_set,
            default_value=False,
            text="Interactive mode to update point coordinates",
            env_id="interactive_mode",
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

                al_bulk: Points = self._build_al_atoms(
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

            self._show_structures(
                carbon_channel_points=carbon_channel.points,
                al_points=al_full_channel_coordinates.points,
                to_build_bonds=to_build_bonds,
                title=structure_folder,
                num_of_min_distances=num_of_min_distances,
                show_al_layers=show_al_layers,
                interactive_mode=interactive_mode,
                # show_coordinates=False,
                # show_indexes=True,
                to_set=to_set,
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

    def get_al_in_channel_details(self, structure_folder: str) -> None:
        """
        Get details of Al atoms in the channel.
        """
        return
        carbon_channel: CarbonHoneycombChannel = AtomsParser.build_carbon_channel(structure_folder)

        number_of_planes: int = int(Inputs.text_input(
            to_set,
            default_value="1",
            text="Number of planes to translate",
            env_id="number_of_planes",
        ))

        al_coordinates: Points = AtomsParser.get_al_channel_coordinates(
            structure_folder, carbon_channel, number_of_planes=number_of_planes, try_to_reflect_al_atoms=True)

        # Prepare data for DataFrame
        data: list[dict] = []

        for index, al_coordinate in enumerate(al_coordinates.points):
            min_dist_to_plane: float = float("inf")
            min_dist_to_carbon: float = np.min(DistanceMeasure.calculate_min_distances(
                np.array([al_coordinate]), carbon_channel.points))

            for plane in carbon_channel.planes:
                dist: float = DistanceMeasure.calculate_distance_from_plane(
                    np.array([al_coordinate]), plane.plane_params)

                if dist < min_dist_to_plane:
                    min_dist_to_plane = dist

            dists_to_al: np.ndarray = DistanceMeasure.calculate_min_distances(
                al_coordinates.points, np.array([al_coordinate]))
            min_dist_to_al: float = np.min(dists_to_al[dists_to_al > 0])  # Exclude self-distance

            # Collect data for each Al coordinate
            data.append({
                # ("Al_atom", "index"): index,
                ("Al_atom", "X"): np.round(al_coordinate[0], 2),
                ("Al_atom", "Y"): np.round(al_coordinate[1], 2),
                ("Al_atom", "Z"): np.round(al_coordinate[2], 2),
                ("min_dist", "to_plane"): np.round(min_dist_to_plane, 2),
                ("min_dist", "to_C"): np.round(min_dist_to_carbon, 2),
                ("min_dist", "to_Al"): np.round(min_dist_to_al, 2),
                **{("dists_to_Al", f"{i}"): np.round(dist, 2) for i, dist in enumerate(dists_to_al)}
            })

        # Create DataFrame with multi-level columns
        df: pd.DataFrame = pd.DataFrame(data)

        # Set multi-level columns
        df.columns = pd.MultiIndex.from_tuples(df.columns)

        file_name: str = f"{structure_folder}_{Constants.filenames.AL_CHANNEL_DETAILS_XLSX_FILE}"

        # Write DataFrame to Excel file
        FileWriter.write_excel_file(
            df=df,
            structure_folder=structure_folder,
            sheet_name="Al atoms in channel details",
            file_name=file_name,
            is_init_data_dir=False,
        )

    def translate_al_to_all_channels(self, structure_folder: str) -> None:
        """
        Read Al plane coordinates from the Excel table and translate the structure to other planes.
        """
        return
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

        interactive_mode: bool = Inputs.bool_input(
            to_set,
            default_value=False,
            text="Interactive mode to update point coordinates",
            env_id="interactive_mode",
        )

        al_channel_coordinates: Points = AtomsParser.get_al_channel_coordinates(
            structure_folder, carbon_channels.pop(0), number_of_planes, try_to_reflect_al_atoms)

        al_coordinates: Points = AlAtomsTranslator.translate_for_all_channels(
            coordinates_carbon=coordinates_carbon,
            carbon_channels=carbon_channels,
            al_channel_coordinates=al_channel_coordinates,
        )

        to_build_bonds = True

        self._show_structures(
            carbon_channel_points=coordinates_carbon.points,
            al_points=al_coordinates.points,
            to_build_bonds=to_build_bonds,
            title=structure_folder,
            num_of_min_distances=num_of_min_distances,
            show_al_layers=show_al_layers,
            interactive_mode=interactive_mode,
            # show_coordinates=False,
            show_indexes=False,
            to_set=to_set,
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

        FileWriter.write_dat_file(
            coordinates_carbon.points,
            structure_folder=structure_folder,
            filename=Constants.filenames.C_ALL_CHANNELS_COORDINATES_DAT_FILE,
        )

    def _build_al_atoms(
            self,
            coordinate_limits: CoordinateLimits,
    ) -> Points:
        al_lattice_type = AlLatticeType(self.al_lattice_type)

        if al_lattice_type.is_cell:
            return IntercalatedChannelBuilder.build_al_coordinates_for_cell(
                to_translate_al=self.to_translate_al,
                al_file=self.file_name,
                coordinate_limits=coordinate_limits,
            )

        else:
            # Fill the volume with aluminium for close-packed lattice
            return IntercalatedChannelBuilder.build_al_coordinates_for_close_packed(
                al_lattice_type=al_lattice_type,
                coordinate_limits=coordinate_limits,
            )

    def _show_structures(
        self,
        carbon_channel_points: np.ndarray,
        al_points: np.ndarray,
        title: str | None = None,
    ) -> None:
        # TODO: refactor

        if self.num_of_al_layers == 1:
            StructureVisualizer.show_two_structures(
                coordinates_first=carbon_channel_points,
                coordinates_second=al_points,
                title=title,
                to_build_bonds=self.to_build_bonds,
                num_of_min_distances=self.bonds_num_of_min_distances,
                skip_first_distances=self.bonds_skip_first_distances,
                to_show_coordinates=self.to_show_coordinates,
                to_show_indexes=self.to_show_indexes,
                # is_interactive_mode=self.interactive_mode,
            )

        elif self.num_of_al_layers == 2:
            # Split the al_points by layers along Oz (by rounded z coordinate)
            al_groups_with_indices: list[tuple[np.ndarray, np.ndarray]] = self._split_atoms_along_z_axis(al_points)

            a_layer_indices: list[int] = []
            b_layer_indices: list[int] = []

            for i, (group, indices) in enumerate(al_groups_with_indices):
                if i % self.num_of_al_layers == 0:
                    a_layer_indices.extend(indices)
                else:
                    b_layer_indices.extend(indices)

            StructureVisualizer.show_structures(
                coordinates_list=[
                    carbon_channel_points,
                    al_points[a_layer_indices],
                    al_points[b_layer_indices],
                ],
                visual_params_list=[
                    VisualizationParams.carbon,
                    VisualizationParams.al,
                    VisualizationParams.al_2,
                ],
                to_build_bonds_list=[
                    self.to_build_bonds,
                    False,
                    False,
                ],
                title=title,
                num_of_min_distances=self.bonds_num_of_min_distances,
                skip_first_distances=self.bonds_skip_first_distances,
                to_show_coordinates=self.to_show_coordinates,
                to_show_indexes=self.to_show_indexes,
                # is_interactive_mode=self.interactive_mode,
                custom_indices_list=[None, a_layer_indices, b_layer_indices],
            )

        elif self.num_of_al_layers == 3:
            al_groups_with_indices: list[tuple[np.ndarray, np.ndarray]] = self._split_atoms_along_z_axis(al_points)

            a_layer_indices: list[int] = []
            b_layer_indices: list[int] = []
            c_layer_indices: list[int] = []

            for i, (group, indices) in enumerate(al_groups_with_indices):
                if i % self.num_of_al_layers == 0:
                    a_layer_indices.extend(indices)
                elif i % self.num_of_al_layers == 1:
                    b_layer_indices.extend(indices)
                else:
                    c_layer_indices.extend(indices)

            StructureVisualizer.show_structures(
                coordinates_list=[
                    carbon_channel_points,
                    al_points[a_layer_indices],
                    al_points[b_layer_indices],
                    al_points[c_layer_indices],
                ],
                visual_params_list=[
                    VisualizationParams.carbon,
                    VisualizationParams.al,
                    VisualizationParams.al_2,
                    VisualizationParams.al_3,
                ],
                to_build_bonds_list=[
                    self.to_build_bonds,
                    False,
                    False,
                    False,
                ],
                title=title,
                num_of_min_distances=self.bonds_num_of_min_distances,
                skip_first_distances=self.bonds_skip_first_distances,
                to_show_coordinates=self.to_show_coordinates,
                to_show_indexes=self.to_show_indexes,
                # is_interactive_mode=self.interactive_mode,
                custom_indices_list=[None, a_layer_indices, b_layer_indices, c_layer_indices],
            )
        else:
            raise NotImplementedError(f"Number of layers {self.num_of_al_layers} is not implemented")

    @staticmethod
    def _split_atoms_along_z_axis(coordinates: np.ndarray) -> list[tuple[np.ndarray, np.ndarray]]:
        """ Returns grouped coordinates with their original indices, grouped by rounded Z coordinate. """
        # Round Z values to the nearest integer or a desired precision
        rounded_z_values: np.ndarray = np.round(coordinates[:, 2], decimals=1)

        # Get unique rounded Z values
        unique_z_values: np.ndarray = np.unique(rounded_z_values)

        # Group points and their indices by their rounded Z coordinate
        grouped_coordinates: list[tuple[np.ndarray, np.ndarray]] = [
            (coordinates[rounded_z_values == z], np.where(rounded_z_values == z)[0])
            for z in unique_z_values
        ]

        return grouped_coordinates
