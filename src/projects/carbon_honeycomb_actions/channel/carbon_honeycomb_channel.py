import numpy as np
from dataclasses import dataclass

from .planes import CarbonHoneycombPlane
from .carbon_honeycomb_channel_actions import CarbonHoneycombChannelActions


@dataclass
class CarbonHoneycombChannel:
    points: np.ndarray

    # def __post_init__(self) -> None:
    #     self.planes: list[CarbonHoneycombPlane] = self._build_planes()

    @property
    def planes(self) -> list[CarbonHoneycombPlane]:
        """
        Returns a list of CarbonHoneycombPlane objects
        representing planes in the honeycomb channel.
        """
        return CarbonHoneycombChannelActions.build_planes(self.points)
