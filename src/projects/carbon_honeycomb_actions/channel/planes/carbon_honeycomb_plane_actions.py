import numpy as np

from .interfaces import CarbonHoneycombPlaneProtocol


class CarbonHoneycombPlaneActions:
    @classmethod
    def define_plane_hexagons(cls, carbon_honeycomb_plane: CarbonHoneycombPlaneProtocol) -> np.ndarray:
        ...
