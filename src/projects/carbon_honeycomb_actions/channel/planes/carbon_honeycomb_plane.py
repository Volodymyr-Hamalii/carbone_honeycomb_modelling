from dataclasses import dataclass

from src.base_structure_classes import Points

from .carbon_honeycomb_plane_actions import CarbonHoneycombPlaneActions
from .plane_polygons import CarbonHoneycombPentagon, CarbonHoneycombHexagon


@dataclass
class CarbonHoneycombPlane(Points):
    # Define direction related plane for honeycomb volume:
    # inside - True
    # outside - False
    direction: bool

    @property
    def pentagons(self) -> list[CarbonHoneycombPentagon]:
        """ List of the plane pentagons. """
        return list(CarbonHoneycombPlaneActions.define_plane_pentagons(self.points))

    @property
    def hexagons(self) -> list[CarbonHoneycombHexagon]:
        """ List of the plane hexagons. """
        return list(CarbonHoneycombPlaneActions.define_plane_hexagons(self.points))
