from pathlib import Path
from numpy import ndarray

from ..constants import Constants
from ..logger import Logger
from .path_builder import PathBuilder


logger = Logger(__name__)


class FileWriter:
    dat_file_first_lines: str = "        210         150           3  1.000000000000000E-002       10000\n" + \
                                "   5.00000000000000             3023   1.00000000000000\n"

    @classmethod
    def write_dat_file(
        cls,
        data_lines: list[str] | ndarray,
        path_to_file: Path | None = None,
        structure_folder: str | None = None,
        overwrite: bool = True,
    ) -> None:

        try:
            if len(data_lines) == 0:
                logger.warning("No data for .dat file.")
                return

            if path_to_file is None:
                if structure_folder is None:
                    raise ValueError(
                        "Provide either path_to_file or structure_folder param for FileWriter.write_dat_file")

                path_to_file = PathBuilder.build_path_to_result_data_file(
                    structure_folder,
                    file=Constants.filenames.INIT_DAT_FILE)

            if overwrite is False and path_to_file.exists():
                # Don't overwrite existing file
                return

            # Convert ndarray to list[str]
            if isinstance(data_lines, ndarray):
                data_lines = [f"{i[0]} {i[1]} {i[2]}" for i in data_lines]

            with Path(path_to_file).open("w") as dat_file:
                dat_file.write(cls.dat_file_first_lines)

                for line in data_lines:
                    dat_file.write(line + "\n")

            logger.info(f"File saved: {path_to_file}")

        except Exception as e:
            logger.error(f".dat file not saved: {e}")
