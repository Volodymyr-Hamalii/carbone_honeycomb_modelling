import numpy as np
import pandas as pd

from src.projects.carbon_honeycomb_actions.channel.planes.carbon_honeycomb_plane import CarbonHoneycombPlane
from src.utils import Constants, Logger, FileReader
from src.base_structure_classes import Points
from src.projects.carbon_honeycomb_actions import (
    CarbonHoneycombChannel,
    CarbonHoneycombActions,
)

from ..intercalated_channel_builder import (
    IntercalatedChannelBuilder,
)
from ..based_on_planes_configs import (
    AtomsBuilder,
    AtomsFilter,
)


logger = Logger("AtomsBuilder")


class AtomsParser:
    @classmethod
    def get_al_plane_coordinates(
            cls,
            structure_folder: str,
            carbon_channel: CarbonHoneycombChannel,
    ) -> Points:
        """
        Read Al coordinates from the file
        or build them if there is no Excel file.
        """

        al_plane_coordinates_df: pd.DataFrame | None = FileReader.read_excel_file(
            structure_folder=structure_folder,
            file_name=Constants.filenames.PLANE_COORDINATS_XLSX_FILE,
            is_init_data_dir=False,
        )

        if al_plane_coordinates_df is not None:
            return cls._parse_al_plane_coordinates_df(al_plane_coordinates_df)

        # logger.warning(f"Excel table with Al atoms for {structure_folder} structure not found. Al atoms builder.")

        # Build atoms
        # carbon_channel: CarbonHoneycombChannel = cls.build_carbon_channel(structure_folder)
        coordinates_al: Points = cls._build_al_plane_coordinates(
            carbon_channel, num_of_planes=1)

        return coordinates_al

    @staticmethod
    def build_carbon_channel(structure_folder: str) -> CarbonHoneycombChannel:
        coordinates_carbon: Points = IntercalatedChannelBuilder.build_carbon_coordinates(
            structure_folder=structure_folder)

        carbon_channels: list[CarbonHoneycombChannel] = CarbonHoneycombActions.split_init_structure_into_separate_channels(
            coordinates_carbon=coordinates_carbon)
        return carbon_channels[0]

    @staticmethod
    def _build_al_plane_coordinates(
            carbon_channel: CarbonHoneycombChannel,
            num_of_planes: int,
    ) -> Points:
        """ Build Al atoms for one plane """
        coordinates_al: Points = AtomsBuilder._build_al_atoms_near_planes(carbon_channel, planes_limit=num_of_planes)
        coordinates_al = AtomsFilter.replace_nearby_atoms_with_one_atom(coordinates_al)
        coordinates_al = AtomsFilter.remove_too_close_atoms(coordinates_al)
        return coordinates_al

    @staticmethod
    def _parse_al_plane_coordinates_df(al_plane_coordinates_df: pd.DataFrame) -> Points:
        """
        Parse al_plane_coordinates_df DataFrame with columns
        i, x_Al, y_Al, z_Al, min_dist_to_Al, Al_1, Al_2, Al_3 ...
        into Points with x_Al, y_Al, z_Al coordinates.

        The points with x_Al, y_Al, z_Al that equals NaN is ignored.
        """
        # Extract the x_Al, y_Al, z_Al columns
        required_columns = ["x_Al", "y_Al", "z_Al"]
        if not all(col in al_plane_coordinates_df.columns for col in required_columns):
            raise ValueError(f"DataFrame must contain columns: {required_columns}")

        # Filter out rows where any of the required coordinates are NaN
        filtered_df: pd.DataFrame = al_plane_coordinates_df.dropna(subset=required_columns)

        # Extract the coordinates as a numpy array
        points_array: np.ndarray = filtered_df[required_columns].to_numpy()

        # Return the Points instance
        return Points(points=points_array)
