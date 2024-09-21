from pathlib import Path

from ..constants import Constants
from ..logger import Logger

from .path_builder import PathBuilder
from .pdb_file_builder import PdbFileBuilder
from .file_writer import FileWriter


logger = Logger(__name__)


class FileConverter:
    @classmethod
    def dat_to_pdb(
        cls,
        structure_folder: str,
        dat_file_name: str = Constants.filenames.INIT_DAT_FILE,
        pdb_file_name: str = Constants.filenames.INIT_PDB_FILE,
    ) -> None:

        # Path to init ljout.dat
        dat_file_path: Path = PathBuilder.build_path_to_init_data_file(
            structure_folder=structure_folder, file=dat_file_name)

        # Convert to ljout-from-init-dat.pdb
        atom_data: list[str] = cls._read_dat_file_as_pdb_file(dat_file_path)

        # Write ljout-from-init-dat.pdb
        pdb_file_path: Path = PathBuilder.build_path_to_result_data_file(
            structure_folder=structure_folder, file=pdb_file_name)
        cls._write_pdb_file(pdb_file_path, atom_data)

        logger.info(f"File saved {pdb_file_path}")

    @classmethod
    def _read_dat_file_as_pdb_file(cls, dat_file_path: Path | str) -> list[str]:
        atom_data: list[str] = []

        with Path(dat_file_path).open("r") as dat_file:
            # Skip the first and second lines
            dat_file.readline()
            dat_file.readline()

            # Read and parse the remaining lines
            atom_id = 1
            for line in dat_file:
                if line.strip():  # Skip empty lines

                    coords: list[str] = line.split()
                    if len(coords) == 3:
                        pdb_line: str = PdbFileBuilder.build_pdb_line(coords, atom_id)
                        atom_data.append(pdb_line)
                        atom_id += 1

        return atom_data
