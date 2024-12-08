import numpy as np
from dataclasses import dataclass

from .interfaces import CarbonHoneycombPlaneProtocol
from .carbon_honeycomb_plane_actions import CarbonHoneycombPlaneActions


@dataclass
class CarbonHoneycombPlane(CarbonHoneycombPlaneProtocol):
    points: np.ndarray
    direction: bool  # To define direction related honeycomb volume

    @property
    def hexagons(self) -> np.ndarray:
        """ List of the plane hexagons. """
        return CarbonHoneycombPlaneActions.define_plane_hexagons(self)
