import json
from typing import Any
from pathlib import Path

import pandas as pd
import numpy as np

from ..constants import Constants
from ..logger import Logger
from .path_builder import PathBuilder


logger = Logger("FileReader")


class FileReader:
    @staticmethod
    def read_json_file(
            folder_path: Path | str,
            file_name: str = Constants.filenames.STRUCTURE_SETTINGS_FILE,
    ) -> Any:
        """
        Read JSON file (by default 'structure_settings.json').
        If folder_path=None -- uses path to 'result_data' folder.
        """

        path_to_file: Path = Path(folder_path) / file_name

        if not path_to_file.exists():
            # logger.warning(f"File {path_to_file} not exists.")
            return None

        data_json: str = Path(path_to_file).read_text(encoding="utf-8")
        return json.loads(data_json)

    @staticmethod
    def read_dat_file(
            structure_folder: str,
            folder_path: Path | str | None = None,
            file_name: str = Constants.filenames.INIT_DAT_FILE,
    ) -> np.ndarray:

        if folder_path is None:
            path_to_file: Path = PathBuilder.build_path_to_init_data_file(
                structure_folder=structure_folder, file=file_name)
        else:
            path_to_file: Path = Path(folder_path) / structure_folder / file_name

        atom_data: list[list[float]] = []

        with Path(path_to_file).open("r") as dat_file:
            # Skip the first and second lines
            dat_file.readline()
            dat_file.readline()

            for line in dat_file:
                if line.strip():  # Skip empty liness
                    coords: list[str] = line.split()
                    if len(coords) == 3:
                        atom_data.append([float(coord) for coord in coords])

        return np.array(atom_data)

    @staticmethod
    def read_excel_file(
            structure_folder: str,
            file_name: str,
            folder_path: Path | str | None = None,
            sheet_name: str | int = 0,
    ) -> pd.DataFrame | None:
        """
        Read an Excel file (by default defined by Constants.filenames.PLANE_COORDINATES_XLSX_FILE).

        Parameters:
        - structure_folder: str, the name of the structure folder.
        - folder_path: Path | str | None, the base folder path. If None, uses default from PathBuilder.
        - file_name: str, the Excel file name to read.
        - sheet_name: str | int, the sheet name or index to read (default is the first sheet).

        Returns:
        - pd.DataFrame: A pandas DataFrame containing the data from the specified Excel file and sheet.

        Raises:
        - FileNotFoundError: If the file does not exist.
        - ValueError: If the file cannot be read as an Excel file.
        """

        if folder_path is None:
            path_to_file: Path = PathBuilder.build_path_to_init_data_file(
                structure_folder=structure_folder, file=file_name)
        else:
            path_to_file: Path = Path(folder_path) / structure_folder / file_name

        # Check if the file exists
        if not path_to_file.exists():
            logger.warning(f"File not found at {path_to_file}")

        try:
            # Read the Excel file into a pandas DataFrame
            df: pd.DataFrame = pd.read_excel(
                path_to_file, sheet_name=sheet_name, engine='openpyxl'
            )
            return df

        except Exception as e:
            logger.error(f"Failed to read file {path_to_file}: {e}")
