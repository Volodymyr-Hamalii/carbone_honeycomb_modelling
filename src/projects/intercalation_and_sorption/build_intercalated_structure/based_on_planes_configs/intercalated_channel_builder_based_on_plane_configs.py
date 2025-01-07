import math
from pathlib import Path

import numpy as np
# from scipy.spatial.distance import cdist
# from scipy.optimize import minimize

from src.utils import Constants, PathBuilder, Logger, execution_time_logger
from src.structure_visualizer import AtomsUniverseBuilder
from src.data_preparation import StructureSettings
from src.base_structure_classes import Points, AlLatticeType, CoordinateLimits
from src.coordinate_operations import PointsRotator
from src.projects.carbon_honeycomb_actions import CarbonHoneycombChannel, CarbonHoneycombPlane

from ...structure_operations import StructureTranslator



logger = Logger("IntercalatedChannelBuilder")


class IntercalatedChannelBuilderBasedOnPlaneConfigs:
    @classmethod
    def build_al_in_carbon(
            cls,
            carbon_channel: CarbonHoneycombChannel,
            structure_settings: StructureSettings,
            equidistant_al_points: bool = True,
    ) -> Points:
        """ Returns coordinates_al """

        coordinates_al: Points = cls.build_coordinates_al(
            carbon_channel=carbon_channel,
            structure_settings=structure_settings)

        # if equidistant_al_points:
        #     # Set Al atoms maximally equidistant from the channel atoms
        #     coordinates_al: Points = AlAtomsSetter.equidistant_points_sets_in_channel(
        #         carbon_channel, coordinates_al, structure_settings)

        return coordinates_al

    @classmethod
    def build_coordinates_al(
        cls,
        carbon_channel: CarbonHoneycombChannel,
        structure_settings: StructureSettings,
    ) -> Points:
        ...
