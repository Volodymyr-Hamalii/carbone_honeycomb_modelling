import MDAnalysis as mda


class StructureUtils:
    @staticmethod
    def write_pdb_from_mda(output_pdb_file: str, atoms) -> None:
        with mda.Writer(output_pdb_file, n_atoms=atoms.n_atoms) as PDB:
            PDB.write(atoms)
