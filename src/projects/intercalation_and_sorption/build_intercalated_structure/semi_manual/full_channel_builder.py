import numpy as np
import pandas as pd

from src.projects.carbon_honeycomb_actions.channel.planes.carbon_honeycomb_plane import CarbonHoneycombPlane
from src.utils import Constants, Logger, FileReader
from src.base_structure_classes import Points
from src.projects.carbon_honeycomb_actions import (
    CarbonHoneycombChannel,
    CarbonHoneycombActions,
)

from ..intercalated_channel_builder import (
    IntercalatedChannelBuilder,
)
from ..based_on_planes_configs import (
    AtomsBuilder,
    AtomsFilter,
)
from .al_atoms_translator import AlAtomsTranslator


class FullChannelBuilder:
    @staticmethod
    def build_full_channel(carbon_channel: CarbonHoneycombChannel) -> Points:
        ...
