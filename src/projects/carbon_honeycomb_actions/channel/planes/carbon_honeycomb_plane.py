import numpy as np
from dataclasses import dataclass

from src.base_structure_classes import Points

from .carbon_honeycomb_plane_actions import CarbonHoneycombPlaneActions
from .carbon_honeycomb_plane_hexagon import CarbonHoneycombHexagon


@dataclass
class CarbonHoneycombPlane(Points):
    # Define direction related plane for honeycomb volume:
    # inside - True
    # outside - False
    direction: bool

    @property
    def hexagons(self) -> np.ndarray:
        """ List of the plane hexagons. """
        return CarbonHoneycombPlaneActions.define_plane_hexagons(self)
