class PdbFileBuilder:
    first_lines: str = "COMPND      BENS NIEUWE KRISTALLEN\n" + \
        "AUTHOR      BWVANDEWAAL    27 04 00\n"

    @staticmethod
    def get_end_lines(num_of_lines: int) -> str:
        return f"TER      {num_of_lines}\n" + \
            "END\n"

    @staticmethod
    def build_pdb_line(coords: list[str], atom_id: int) -> str:
        x, y, z = map(float, coords)
        return f"ATOM  {atom_id:>5} C            1    {x:>8.3f}{y:>8.3f}{z:>8.3f}   1.000   0.000\n"
