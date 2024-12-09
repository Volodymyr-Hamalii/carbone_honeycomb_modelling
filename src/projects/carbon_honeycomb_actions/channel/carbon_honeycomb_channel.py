import numpy as np
from dataclasses import dataclass

from .planes import CarbonHoneycombPlane


@dataclass
class CarbonHoneycombChannel:
    points: np.ndarray

    @property
    def planes(self) -> list[CarbonHoneycombPlane]:
        ...
        # planes: list[CarbonHoneycombPlane]
