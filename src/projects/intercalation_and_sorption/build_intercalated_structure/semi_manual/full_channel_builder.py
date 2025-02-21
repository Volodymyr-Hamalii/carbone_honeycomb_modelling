import numpy as np

from src.utils import Constants, Logger
from src.base_structure_classes import Points
from src.coordinate_operations import PointsFilter
from src.projects.carbon_honeycomb_actions import (
    CarbonHoneycombChannel,
)
from ..by_variance import (
    AlAtomsFilter,
)
from .al_atoms_translator import AlAtomsTranslator


logger = Logger(__name__)


class FullChannelBuilder:
    @classmethod
    def build_full_channel(
            cls,
            carbon_channel: CarbonHoneycombChannel,
            al_channel_planes_coordinates: Points,
            al_crystal: Points
    ) -> Points:

        al_bulk_coordinates: Points = cls._filter_al_bulk_coordinates(
            carbon_channel=carbon_channel,
            al_channel_planes_coordinates=al_channel_planes_coordinates,
            al_crystal=al_crystal,
        )

        return Points(
            points=np.vstack([al_channel_planes_coordinates.points, al_bulk_coordinates.points])
        )

    @staticmethod
    def _filter_al_bulk_coordinates(
        carbon_channel: CarbonHoneycombChannel,
        al_channel_planes_coordinates: Points,
        al_crystal: Points
    ) -> Points:

        al_bulk_coordinates: Points = AlAtomsFilter.filter_atoms_related_clannel_planes(
            coordinates_al=al_crystal,
            carbon_channel=carbon_channel,
            distance_from_plane=Constants.phys.al.MIN_RECOMENDED_DIST_BETWEEN_ATOMS,
        )

        return al_bulk_coordinates
