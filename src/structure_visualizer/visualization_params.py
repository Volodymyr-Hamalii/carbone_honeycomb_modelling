from dataclasses import dataclass


@dataclass
class StructureVisualParams:
    color_atoms: str
    color_bonds: str
    size: int
    transparency: float
    set_equal_scale: bool
    label: str
    show_coordinates: bool


class VisualizationParams:
    carbon = StructureVisualParams(
        label="Carbon",
        color_atoms="#0500a4",
        color_bonds="#00065f",
        transparency=0.25,
        set_equal_scale=True,
        size=200,
        show_coordinates=False,
    )

    al = StructureVisualParams(
        label="Aluminum",
        color_atoms="#e00000",
        color_bonds="#500000",
        transparency=0.5,
        set_equal_scale=False,
        size=400,
        show_coordinates=True,
    )
