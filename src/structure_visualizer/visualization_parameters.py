class StructureVisualParameters:
    def __init__(
            self,
            color_atoms: str,
            color_bonds: str,
            size: int,
            transparency: float,
    ) -> None:
        self.color_atoms: str = color_atoms
        self.color_bonds: str = color_bonds
        self.size: int = size
        self.transparency: float = transparency


class VisualizationParameters:
    carbone = StructureVisualParameters(
        color_atoms="#0500a4",
        color_bonds="#00065f",
        size=50,
        transparency=0.75,
    )

    al = StructureVisualParameters(
        color_atoms="#e00000",
        color_bonds="#500000",
        size=100,
        transparency=0.5,
    )
