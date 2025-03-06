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
    show_indexes: bool


class ColorParams:
    carbon_atoms: str = "#0500a4"
    carbon_bonds: str = "#00065f"
    aluminum_atoms: str = "#e00000"
    aluminum_bonds: str = "#500000"
    aluminum_2_atoms: str = "#00d11d"
    aluminum_2_bonds: str = "#004309"


class VisualizationParams:
    carbon = StructureVisualParams(
        label="Carbon",
        color_atoms=ColorParams.carbon_atoms,
        color_bonds=ColorParams.carbon_bonds,
        transparency=0.25,
        set_equal_scale=True,
        size=200,
        show_coordinates=False,
        show_indexes=False,
    )

    al = StructureVisualParams(
        label="Aluminum",
        color_atoms=ColorParams.aluminum_atoms,
        color_bonds=ColorParams.aluminum_bonds,
        transparency=0.5,
        set_equal_scale=False,
        size=400,
        show_coordinates=False,
        show_indexes=True,
    )

    al_2 = StructureVisualParams(
        label="Aluminum",
        color_atoms=ColorParams.aluminum_2_atoms,
        color_bonds=ColorParams.aluminum_2_bonds,
        transparency=0.5,
        set_equal_scale=False,
        size=400,
        show_coordinates=False,
        show_indexes=True,
    )
