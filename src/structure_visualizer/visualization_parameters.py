class StructureVisualParameters:
    def __init__(
            self,
            color_atoms: str,
            color_bonds: str,
            size: int,
            transparency: float,
            label: str,
    ) -> None:
        self.color_atoms: str = color_atoms
        self.color_bonds: str = color_bonds
        self.size: int = size
        self.transparency: float = transparency
        self.label: str = label


class VisualizationParameters:
    carbone = StructureVisualParameters(
        color_atoms="#0500a4",
        color_bonds="#00065f",
        size=200,
        transparency=0.35,
        label="Carbon",
    )

    al = StructureVisualParameters(
        color_atoms="#e00000",
        color_bonds="#500000",
        size=400,
        transparency=0.5,
        label="Aluminum",
    )
