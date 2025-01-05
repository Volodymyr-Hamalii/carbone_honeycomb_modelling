import numpy as np
from dataclasses import dataclass

from src.base_structure_classes import Points

from .planes import CarbonHoneycombPlane
from .carbon_honeycomb_channel_actions import CarbonHoneycombChannelActions


@dataclass
class CarbonHoneycombChannel(Points):
    # def __post_init__(self) -> None:
    #     self.planes: list[CarbonHoneycombPlane] = self._build_planes()

    @property
    def planes(self) -> list[CarbonHoneycombPlane]:
        """
        A list of CarbonHoneycombPlane objects
        representing planes in the honeycomb channel.
        """
        return CarbonHoneycombChannelActions.build_planes(self.points)
