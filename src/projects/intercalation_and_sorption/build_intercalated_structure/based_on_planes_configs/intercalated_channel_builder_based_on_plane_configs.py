import numpy as np

from src.utils import Logger
from src.data_preparation import StructureSettings
from src.base_structure_classes import Points
from src.projects.carbon_honeycomb_actions import CarbonHoneycombChannel

from .atoms_builder import AtomsBuilder
from .atoms_filter import AtomsFilter


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

        coordinates_al: Points = AtomsBuilder._build_al_atoms_near_planes(carbon_channel)
        coordinates_al: Points = AtomsFilter.replace_nearby_atoms_with_one_atom(coordinates_al)
        coordinates_al: Points = AtomsFilter.remove_too_close_atoms(coordinates_al)

        return coordinates_al
