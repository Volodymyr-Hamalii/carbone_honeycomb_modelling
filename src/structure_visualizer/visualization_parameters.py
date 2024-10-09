from dataclasses import dataclass


@dataclass
class StructureVisualParameters:
    color_atoms: str
    color_bonds: str
    size: int
    transparency: float
    set_equal_scale: bool
    label: str


class VisualizationParameters:
    carbon = StructureVisualParameters(
        label="Carbon",
        color_atoms="#0500a4",
        color_bonds="#00065f",
        transparency=0.25,
        set_equal_scale=True,
        size=200,
    )

    al = StructureVisualParameters(
        label="Aluminum",
        color_atoms="#e00000",
        color_bonds="#500000",
        transparency=0.5,
        set_equal_scale=False,
        size=400,
    )
