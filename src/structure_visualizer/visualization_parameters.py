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
        color_atoms="#3d83fe",
        color_bonds="#a5c6ff",
        size=50,
        transparency=0.75,
    )

    al = StructureVisualParameters(
        color_atoms="#C70039",
        color_bonds="#ffa5a5",
        size=100,
        transparency=0.25,
    )
