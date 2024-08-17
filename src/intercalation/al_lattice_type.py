class AlLatticeType:
    def __init__(self, al_structure: str = "cell") -> None:
        """ AL structure to fill (can be 'cell' for cubic cell, 'FCC' or 'HCP' for planes) """
        self.is_cell: bool = al_structure == "cell"  # Cubic cell
        self.is_fcc: bool = al_structure == "FCC"  # Face-centered cubic
        self.is_hcp: bool = al_structure == "HCP"  # Hexagonal close-packed

        if not (self.is_cell or self.is_fcc or self.is_hcp):
            raise ValueError("Pleae, provide available Al structure ('cell', 'FCC' or 'HCP').")

    @staticmethod
    def get_info() -> str:
        return "AL structure to fill (can be 'cell' for cubic cell, 'FCC' or 'HCP' for planes)"
