import numpy as np
from dataclasses import dataclass

from src.base_structure_classes import Points

from .carbon_honeycomb_plane_hexagon_actions import CarbonHoneycombHexagonActions


@dataclass
class CarbonHoneycombHexagon(Points):
    @property
    def center(self) -> np.ndarray:
        """ The coordinates of the hexagon center. """
        return CarbonHoneycombHexagonActions.define_center_coordinates(self.points)
