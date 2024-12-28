import numpy as np
from dataclasses import dataclass

from src.base_structure_classes import Points

from .carbon_honeycomb_plane_polygon_actions import CarbonHoneycombPolygonActions


@dataclass
class CarbonHoneycombPolygon(Points):
    @property
    def center(self) -> np.ndarray:
        """ The coordinates of the hexagon center. """
        return CarbonHoneycombPolygonActions.define_polygon_center_coordinates(self.points)


@dataclass
class CarbonHoneycombHexagon(CarbonHoneycombPolygon):
    pass


@dataclass
class CarbonHoneycombPentagon(CarbonHoneycombPolygon):
    pass
