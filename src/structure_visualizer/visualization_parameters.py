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
        label="Carbon",
        color_atoms="#0500a4",
        color_bonds="#00065f",
        transparency=0.25,
        size=200,
    )

    al = StructureVisualParameters(
        label="Aluminum",
        color_atoms="#e00000",
        color_bonds="#500000",
        transparency=0.5,
        size=400,
    )
