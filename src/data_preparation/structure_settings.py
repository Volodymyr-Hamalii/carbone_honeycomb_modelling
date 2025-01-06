from dataclasses import dataclass


@dataclass
class StructureSettings:
    distance_from_plane: float = 0
    max_distance_to_carbon_atoms: float = 0
    direction_related_center: bool = False  # TODO[25.01.07]: set this logic and remove this
