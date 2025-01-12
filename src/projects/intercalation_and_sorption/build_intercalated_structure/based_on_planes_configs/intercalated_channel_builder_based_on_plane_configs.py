import numpy as np

from src.utils import Logger
from src.data_preparation import StructureSettings
from src.base_structure_classes import Points
from src.projects.carbon_honeycomb_actions import CarbonHoneycombChannel

from .atoms_builder import AtomsBuilder
from .atoms_filter import AtomsFilter
from .atom_configurator import AtomConfigurator


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
        coordinates_al = AtomsFilter.replace_nearby_atoms_with_one_atom(coordinates_al)
        coordinates_al = AtomsFilter.remove_too_close_atoms(coordinates_al)
        coordinates_al = AtomConfigurator.space_al_atoms_equidistantly(coordinates_al, carbon_channel)

        return coordinates_al
