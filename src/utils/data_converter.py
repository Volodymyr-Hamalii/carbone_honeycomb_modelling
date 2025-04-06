import numpy as np
import pandas as pd
from .logger import Logger


logger = Logger("DataConverter")


class DataConverter:
    _dat_file_first_lines: str = "        210         150           3  1.000000000000000E-002       10000\n" + \
        "   5.00000000000000             3023   1.00000000000000\n"

    @classmethod
    def convert_df_to_dat(cls, df: pd.DataFrame) -> list[str]:
        """ Return list of strings with .dat file lines. """
        dat_lines: list[str] = [cls._dat_file_first_lines]
        for _, row in df.iterrows():
            dat_lines.append(f"{row[0]:.6f}\t{row[1]:.6f}\t{row[2]:.6f}\n")
        return dat_lines

    @classmethod
    def convert_df_to_pdb(cls, df: pd.DataFrame) -> list[str]:
        """ Return list of strings with .pdb file lines. """
        pdb_lines: list[str] = []
        for i, row in df.iterrows():
            pdb_line: str = cls._build_pdb_line(row.to_numpy(), i + 1)  # type: ignore
            pdb_lines.append(pdb_line)
        return pdb_lines

    @classmethod
    def convert_ndarray_to_dat(cls, coordinates: np.ndarray) -> list[str]:
        """ Return list of strings with .dat file lines. """
        dat_lines: list[str] = [cls._dat_file_first_lines]
        for row in coordinates:
            dat_lines.append(f"{row[0]:.6f}\t{row[1]:.6f}\t{row[2]:.6f}")
        return dat_lines

    @classmethod
    def convert_ndarray_to_pdb(cls, coordinates: np.ndarray) -> list[str]:
        """ Return list of strings with .pdb file lines. """
        pdb_lines: list[str] = []
        for i, row in enumerate(coordinates):
            pdb_line: str = cls._build_pdb_line(row, i + 1)
            pdb_lines.append(pdb_line)
        return pdb_lines

    @staticmethod
    def _build_pdb_line(row: np.ndarray, index: int) -> str:
        """ Build a PDB line from a row of data. """
        return f"ATOM  {index:5d}  CA  ALA A   1    " \
            f"{row[0]:8.3f}{row[1]:8.3f}{row[2]:8.3f}  1.00  0.00           C\n"
