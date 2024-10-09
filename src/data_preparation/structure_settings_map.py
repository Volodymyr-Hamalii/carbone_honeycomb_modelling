from dataclasses import dataclass


@dataclass
class ChannelLimits:
    x_min: float
    x_max: float
    y_min: float
    y_max: float
    z_min: float | None = None
    z_max: float | None = None


@dataclass
class ChannelPoints:
    points: list[list[float]]
    direction: int


@dataclass
class StructureSettings:
    channel_limits: ChannelLimits
    points_to_set_channel_planes: list[ChannelPoints]
    distance_from_plane: float = 0
    max_distance_to_carbon_atoms: float = 0
    al_lattice_parameter: float = 0
