from dataclasses import dataclass


@dataclass
class StructureVisualParameters:
    color_atoms: str
    color_bonds: str
    size: int
    transparency: float
    label: str


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
