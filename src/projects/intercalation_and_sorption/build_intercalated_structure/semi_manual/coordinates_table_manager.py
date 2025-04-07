from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from src.utils import Constants, Logger, FileReader, FileWriter
from src.base_structure_classes import Points
from src.projects.carbon_honeycomb_actions import (
    CarbonHoneycombChannel,
)
from .atoms_parser import AtomsParser


logger = Logger("CoordinatesTableManager")


class CoordinatesTableManager:
    @classmethod
    def update_plane_tbl_file(
            cls,
            structure_folder: str,
            carbon_channel: CarbonHoneycombChannel,
            number_of_planes: int,
    ) -> Path:
        al_plane_coordinates: Points = AtomsParser.get_al_plane_coordinates(
            structure_folder, carbon_channel, number_of_planes)
        df: pd.DataFrame = cls._build_updated_df(al_plane_coordinates)

        path_to_file: Path | None = FileWriter.write_excel_file(
            df=df,
            structure_folder=structure_folder,
            sheet_name="Al atoms for the plane",
            file_name=Constants.filenames.AL_PLANE_COORDINATES_XLSX_FILE,
            is_init_data_dir=False,
        )

        if path_to_file is None:
            raise IOError(f"Failed to write {Constants.filenames.AL_PLANE_COORDINATES_XLSX_FILE} file.")

        return path_to_file

    @classmethod
    def update_full_channel_tbl_file(cls, structure_folder: str,) -> Path:
        al_channel_coordinates_df: pd.DataFrame | None = FileReader.read_excel_file(
            structure_folder=structure_folder,
            file_name=Constants.filenames.AL_FULL_CHANNEL_COORDINATES_XLSX_FILE,
            is_init_data_dir=False,
        )

        if al_channel_coordinates_df is None:
            raise FileNotFoundError(
                f"Excel file with Al atoms for the full channel not found in {structure_folder}.")

        al_channel_coordinates: Points = AtomsParser._parse_al_coordinates_df(al_channel_coordinates_df)
        df: pd.DataFrame = cls._build_updated_df(al_channel_coordinates)

        path_to_file: Path | None = FileWriter.write_excel_file(
            df=df,
            structure_folder=structure_folder,
            sheet_name="Al atoms for the full channel",
            file_name=Constants.filenames.AL_FULL_CHANNEL_COORDINATES_XLSX_FILE,
            is_init_data_dir=False,
        )

        if path_to_file is None:
            raise IOError(f"Failed to write {Constants.filenames.AL_FULL_CHANNEL_COORDINATES_XLSX_FILE} file.")

        return path_to_file

    @staticmethod
    def _build_updated_df(al_points: Points) -> pd.DataFrame:
        """
        Build pd.DataFrame based on the points with the colums
        i, x_Al, y_Al, z_Al, min_dist_to_Al, Al_1, Al_2, Al_3 ...

        where
        i - index of the point in the al_points.points,
        x_Al, y_Al, z_Al - coordinates of the point,
        Al_1, Al_2, Al_3, ... - distance from the point to other points.
        """

        # Extract the array of points
        points: np.ndarray = al_points.points

        # Validate the shape of the points
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("Points must be a 2D array with shape (N, 3).")

        # Compute the distance matrix
        dists: np.ndarray = cdist(points, points)

        # Compute the minimum distance to each point (ignoring the self-distance)
        min_dist_to_Al: np.ndarray = np.min(dists + np.diag([np.inf] * len(points)), axis=1)

        # Prepare data for the DataFrame
        data: dict = {
            "i": np.arange(len(points)),  # Point index
            "x_Al": points[:, 0],        # X-coordinates
            "y_Al": points[:, 1],        # Y-coordinates
            "z_Al": points[:, 2],        # Z-coordinates
            "min_dist_to_Al": min_dist_to_Al,  # Minimum distance to other points
        }

        # Add columns for distances to other points
        for j in range(dists.shape[1]):
            data[f"Al_{j}"] = dists[:, j]

        # Build the DataFrame
        df = pd.DataFrame(data)

        return df
