from dataclasses import dataclass


@dataclass
class StructureVisualParams:
    color_atoms: str
    color_bonds: str
    size: int
    bonds_width: float
    transparency: float
    transparency_bonds: float
    set_equal_scale: bool
    label: str
    show_coordinates: bool
    show_indexes: bool


class Colors:
    carbon_atoms: str = "#0500a4"
    carbon_bonds: str = "#00065f"

    aluminum_1_atoms: str = "#e00000"
    aluminum_1_bonds: str = "#500000"
    aluminum_2_atoms: str = "#00d11d"
    aluminum_2_bonds: str = "#004309"
    aluminum_3_atoms: str = "#1b9afa"
    aluminum_3_bonds: str = "##0062ad"

    black: str = "#000000"
    gray100: str = "#454545"
    gray200: str = "#6a6a6a"
    gray300: str = "#919191"
    gray400: str = "#afafaf"

    # # Uncomment to use B&W colors
    # carbon_atoms: str = black
    # carbon_bonds: str = black
    # aluminum_atoms: str = gray100
    # aluminum_bonds: str = gray100
    # aluminum_2_atoms: str = gray300
    # aluminum_2_bonds: str = gray300


class VisualizationParams:
    carbon = StructureVisualParams(
        label="Carbon",

        color_atoms=Colors.carbon_atoms,
        transparency=0.2,
        size=100,

        color_bonds=Colors.carbon_bonds,
        transparency_bonds=1,
        bonds_width=0.5,

        set_equal_scale=True,
        show_coordinates=False,
        show_indexes=False,
    )

    al = StructureVisualParams(
        label="Aluminum",

        color_atoms=Colors.aluminum_1_atoms,
        transparency=0.5,
        size=400,

        color_bonds=Colors.aluminum_1_bonds,
        transparency_bonds=1,
        bonds_width=1,

        set_equal_scale=False,
        show_coordinates=False,
        show_indexes=True,
    )

    al_2 = StructureVisualParams(
        label="Aluminum",

        color_atoms=Colors.aluminum_2_atoms,
        transparency=0.5,
        size=400,

        color_bonds=Colors.aluminum_2_bonds,
        transparency_bonds=1,
        bonds_width=1,

        set_equal_scale=False,
        show_coordinates=False,
        show_indexes=True,
    )

    al_3 = StructureVisualParams(
        label="Aluminum",

        color_atoms=Colors.aluminum_3_atoms,
        transparency=0.5,
        size=400,

        color_bonds=Colors.aluminum_3_bonds,
        transparency_bonds=1,
        bonds_width=1,

        set_equal_scale=False,
        show_coordinates=False,
        show_indexes=True,
    )
