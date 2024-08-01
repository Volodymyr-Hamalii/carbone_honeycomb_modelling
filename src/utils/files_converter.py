from .constants import Constants
from .path_builder import PathBuilder
from .logger import Logger


logger = Logger(__name__)


class FilesConverter:
    @classmethod
    def dat_to_pdb(
        cls,
        structure_folder: str,
        dat_file_name: str = Constants.path.INIT_DAT_FILE,
        pdb_file_name: str = Constants.path.INIT_PDB_FILE,
    ):

        # Path to init ljout.dat
        dat_file_path: str = PathBuilder.build_path_to_init_data_file(
            structure_folder=structure_folder, file=dat_file_name)

        # Convert to ljout-result.pdb
        atom_data: list[str] = cls._build_pdb_file(dat_file_path)

        # Write ljout-result.pdb
        pdb_file_path: str = PathBuilder.build_path_to_result_data_file(
            structure_folder=structure_folder, file=pdb_file_name)
        cls._write_pdb_file(pdb_file_path, atom_data)

        logger.info(f"Generated {pdb_file_path}")

    @classmethod
    def _build_pdb_file(cls, dat_file_path: str) -> list[str]:
        atom_data: list[str] = []

        with open(dat_file_path, 'r') as dat_file:
            # Skip the first and second lines
            dat_file.readline()
            dat_file.readline()

            # Read and parse the remaining lines
            atom_id = 1
            for line in dat_file:
                if line.strip():  # Skip empty lines

                    coords: list[str] = line.split()
                    if len(coords) == 3:
                        pdb_line = cls._build_pdb_line(coords, atom_id)
                        atom_data.append(pdb_line)
                        atom_id += 1

        return atom_data

    @staticmethod
    def _write_pdb_file(pdb_file_path: str, atom_data: list[str]):
        with open(pdb_file_path, 'w') as pdb_file:
            # Write the PDB file header
            pdb_file.write("COMPND      BENS NIEUWE KRISTALLEN\n")
            pdb_file.write("AUTHOR      BWVANDEWAAL    27 04 00\n")

            # Write the atom data
            for atom_line in atom_data:
                pdb_file.write(atom_line)

            # Write the PDB file footer
            pdb_file.write(f"TER      {len(atom_data)}\n")
            pdb_file.write("END\n")

    @staticmethod
    def _build_pdb_line(coords: list[str], atom_id: int) -> str:
        x, y, z = map(float, coords)
        return f"ATOM  {atom_id:>5} C            1    {x:>8.3f}{y:>8.3f}{z:>8.3f}   1.000   0.000\n"
