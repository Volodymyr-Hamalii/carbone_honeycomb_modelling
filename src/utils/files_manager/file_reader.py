import json
from typing import Any
from pathlib import Path

import pandas as pd
import numpy as np

from ..constants import Constants
from ..logger import Logger
from .file_manager import FileManager


logger = Logger("FileReader")


class FileReader(FileManager):
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

    @classmethod
    def read_dat_file(
            cls,
            structure_folder: str,
            folder_path: Path | str | None = None,
            file_name: str = Constants.filenames.INIT_DAT_FILE,
            is_init_data_dir: bool = True,
    ) -> np.ndarray:

        path_to_file: Path = cls._get_path_to_file(structure_folder, file_name, folder_path, is_init_data_dir)

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

    @classmethod
    def read_excel_file(
            cls,
            structure_folder: str,
            file_name: str,
            folder_path: Path | str | None = None,
            sheet_name: str | int = 0,
            is_init_data_dir: bool = True,
    ) -> pd.DataFrame | None:
        """
        Read an Excel file (by default defined by Constants.filenames.PLANE_COORDINATES_XLSX_FILE).

        Parameters:
        - structure_folder: str, the name of the structure folder.
        - folder_path: Path | str | None, the base folder path. If None, uses default from PathBuilder.
        - file_name: str, the Excel file name to read.
        - sheet_name: str | int, the sheet name or index to read (default is the first sheet).
        - is_init_data_dir: bool | None: to build path to the specific dir. If it's False - builds path to result data.

        Returns:
        - pd.DataFrame: A pandas DataFrame containing the data from the specified Excel file and sheet.
        """

        path_to_file: Path = cls._get_path_to_file(structure_folder, file_name, folder_path, is_init_data_dir)

        # Check if the file exists
        if not path_to_file.exists():
            logger.warning(f"File not found at {path_to_file}")

        try:
            # Read the Excel file into a pandas DataFrame
            df: pd.DataFrame = pd.read_excel(
                path_to_file, sheet_name=sheet_name, engine='openpyxl'
            )
            return df

        except FileNotFoundError:
            logger.warning(f"File {path_to_file} not exists.")

        except Exception as e:
            logger.error(f"Failed to read file {path_to_file}: {e}")
