from dataclasses import dataclass
from functools import cached_property

from src.base_structure_classes import FlatFigure

from .carbon_honeycomb_plane_actions import CarbonHoneycombPlaneActions
from .plane_polygons import CarbonHoneycombPentagon, CarbonHoneycombHexagon


@dataclass
class CarbonHoneycombPlane(FlatFigure):
    # Define direction related plane for honeycomb volume:
    # inside - True
    # outside - False
    direction: bool

    @cached_property
    def pentagons(self) -> list[CarbonHoneycombPentagon]:
        """ List of the plane pentagons. """
        return list(CarbonHoneycombPlaneActions.define_plane_pentagons(self.points))

    @cached_property
    def hexagons(self) -> list[CarbonHoneycombHexagon]:
        """ List of the plane hexagons. """
        return list(CarbonHoneycombPlaneActions.define_plane_hexagons(self.points))

    def copy(self) -> "CarbonHoneycombPlane":
        """ Returns a new CarbonHoneycombPlane instance. """
        return CarbonHoneycombPlane(
            points=self.points.copy(),
            direction=self.direction,
        )
