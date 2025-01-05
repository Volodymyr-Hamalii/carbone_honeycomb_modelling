from dataclasses import dataclass
from functools import cached_property

from src.base_structure_classes import Points

from .planes import CarbonHoneycombPlane
from .carbon_honeycomb_channel_actions import CarbonHoneycombChannelActions


@dataclass
class CarbonHoneycombChannel(Points):
    @cached_property
    def planes(self) -> list[CarbonHoneycombPlane]:
        """
        A list of CarbonHoneycombPlane objects
        representing planes in the honeycomb channel.
        """
        return CarbonHoneycombChannelActions.build_planes(self.points)
