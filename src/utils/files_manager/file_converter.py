# from pathlib import Path
# import numpy as np
# from numpy import ndarray

# from ..constants import Constants
# from ..logger import Logger

# from .path_builder import PathBuilder
# from .pdb_file_builder import PdbFileBuilder
# from .file_reader import FileReader
# from .file_writer import FileWriter


# logger = Logger("FileConverter")


# class FileConverter:
#     @classmethod
#     def dat_to_pdb(
#         cls,
#         structure_folder: str,
#         dat_file_name: str = Constants.file_names.INIT_DAT_FILE,
#         pdb_file_name: str = Constants.file_names.INIT_PDB_FILE,
#     ) -> None:

#         # Read and parse the file
#         atom_data: ndarray = FileReader.read_dat_file(
#             structure_folder=structure_folder, file_name=dat_file_name)

#         # Convert to ljout-from-init-dat.pdb
#         pdb_lines: list[str] = cls._convert_ndarray_to_pdb_lines(atom_data)

#         # Write ljout-from-init-dat.pdb
#         pdb_file_path: Path = PathBuilder.build_path_to_result_data_file(
#             structure_folder=structure_folder, file=pdb_file_name)

#         FileWriter.write_pdb_file(pdb_lines, pdb_file_path)
#         logger.info(f"File saved {pdb_file_path}")

#     @staticmethod
#     def _convert_ndarray_to_pdb_lines(atom_data: ndarray) -> list[str]:
#         pdb_lines: list[str] = []

#         for i, line in enumerate(atom_data):
#             pdb_line: str = PdbFileBuilder.build_pdb_line(line, i+1)
#             pdb_lines.append(pdb_line)

#         return pdb_lines
