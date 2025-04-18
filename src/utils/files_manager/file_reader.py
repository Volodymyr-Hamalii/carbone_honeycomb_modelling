import json
from typing import Any
from pathlib import Path

import pandas as pd
import numpy as np
import MDAnalysis as mda

from ..constants import Constants
from ..logger import Logger
from .file_manager import FileManager


logger = Logger("FileReader")


class FileReader(FileManager):
    @staticmethod
    def read_list_of_dirs(
            folder_path: Path | str = Constants.path.INIT_DATA_PATH,
    ) -> list[str]:
        """ Read a list of directories in the given folder path. By default uses 'init_data' folder. """
        try:
            return sorted(
                [
                    dir.name for dir in Path(folder_path).iterdir()
                    if dir.is_dir()
                ]
            )
        except FileNotFoundError:
            logger.error(f"Folder {folder_path} not found.")
            return []
        except Exception as e:
            logger.error(f"Failed to read list of directories in {folder_path}: {e}")
            return []

    @staticmethod
    def read_list_of_files(
            folder_path: Path | str,
            format: str | None = None,
            to_include_nested_files: bool = False,
            to_append_parent_dir: bool = False,
    ) -> list[str]:
        """
        Read a list of files in the given folder path. By default uses '.xlsx' format.
        If to_include_nested_files is True, it will include files from nested folders.
        """
        try:
            file_names: list[str] = []

            for file in Path(folder_path).iterdir():
                if file.is_file() and (
                    file.name.endswith(format) if format else True
                ):
                    if to_append_parent_dir:
                        file_names.append(file.parent.name + "/" + file.name)
                    else:
                        file_names.append(file.name)

                if to_include_nested_files and file.is_dir():
                    file_names.extend(
                        FileReader.read_list_of_files(
                            file,
                            format=format,
                            to_include_nested_files=to_include_nested_files,
                            to_append_parent_dir=True,
                        )
                    )

            return sorted(file_names)

        except FileNotFoundError:
            logger.error(f"Folder {folder_path} not found.")
            return []
        except Exception as e:
            logger.error(f"Failed to read list of files in {folder_path}: {e}")
            return []

    @staticmethod
    def read_json_file(
            folder_path: Path | str,
            file_name: str = Constants.file_names.STRUCTURE_SETTINGS_FILE,
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
            file_name: str = Constants.file_names.INIT_DAT_FILE,
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

    @staticmethod
    def read_pdb_file(
            structure_folder: str,
            folder_path: Path | str | None = None,
            file_name: str = Constants.file_names.INIT_PDB_FILE,
            is_init_data_dir: bool = True,
    ) -> np.ndarray:
        """
        Read a PDB file and return its atomic coordinates as a NumPy array.

        Parameters:
        - structure_folder: str, the name of the structure folder.
        - folder_path: Path | str | None, the base folder path. If None, uses default from PathBuilder.
        - file_name: str, the PDB file name to read.
        - is_init_data_dir: bool | None: to build path to the specific dir. If it's False - builds path to result data.

        Returns:
        - np.ndarray: A NumPy array containing the atomic coordinates from the PDB file.
        """
        path_to_file: Path = FileManager._get_path_to_file(structure_folder, file_name, folder_path, is_init_data_dir)

        # atom_data: list[list[float]] = []

        try:
            # with path_to_file.open("r") as pdb_file:
            #     for line in pdb_file:
            #         if line.startswith(("ATOM", "HETATM")):
            #             x: float = float(line[30:38].strip())
            #             y: float = float(line[38:46].strip())
            #             z: float = float(line[46:54].strip())
            #             atom_data.append([x, y, z])

            # Load a structure from a file
            u = mda.Universe(path_to_file)

            # Access atoms and coordinates
            atoms = u.atoms
            points: np.ndarray = atoms.positions  # type: ignore

            return np.array(points)

        except FileNotFoundError:
            logger.error(f"PDB file not found at {path_to_file}")
            return np.array([])

        except Exception as e:
            logger.error(f"Failed to read PDB file {path_to_file}: {e}")
            return np.array([])

    @classmethod
    def read_excel_file(
            cls,
            structure_folder: str,
            file_name: str,
            folder_path: Path | str | None = None,
            sheet_name: str | int = 0,
            is_init_data_dir: bool = True,
            to_print_warning: bool = True,
    ) -> pd.DataFrame | None:
        """
        Read an Excel file.

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
            if to_print_warning:
                logger.warning(f"File not found at {path_to_file}")
            return None

        try:
            # Read the Excel file into a pandas DataFrame
            df: pd.DataFrame = pd.read_excel(
                path_to_file, sheet_name=sheet_name, engine='openpyxl'
            )
            return df

        except FileNotFoundError:
            if to_print_warning:
                logger.warning(f"File {path_to_file} not exists.")

        except Exception as e:
            logger.error(f"Failed to read file {path_to_file}: {e}")

    @classmethod
    def read_init_data_file(cls, structure_folder: str, file_name: str) -> np.ndarray:
        file_format: str = file_name.split(".")[-1].lower()

        if file_format == "pdb":
            return cls.read_pdb_file(structure_folder, file_name=file_name)
        elif file_format == "dat":
            return cls.read_dat_file(structure_folder, file_name=file_name)
        elif file_format == "xlsx":
            carbon_points_df: pd.DataFrame | None = cls.read_excel_file(structure_folder, file_name=file_name)
            if carbon_points_df is None:
                raise IOError(f"Failed to read file: {file_name}")
            return carbon_points_df.to_numpy()
        else:
            raise ValueError(f"Unsupported file format: {file_format}")
