import numpy as np
from dataclasses import dataclass

from src.base_structure_classes import FlatFigure

from .carbon_honeycomb_plane_polygon_actions import CarbonHoneycombPolygonActions


@dataclass
class CarbonHoneycombPolygon(FlatFigure):
    pass


@dataclass
class CarbonHoneycombHexagon(CarbonHoneycombPolygon):
    pass


@dataclass
class CarbonHoneycombPentagon(CarbonHoneycombPolygon):
    pass
