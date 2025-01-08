import math
from pathlib import Path

import numpy as np

from src.utils import Constants, PathBuilder, Logger, execution_time_logger
from src.structure_visualizer import AtomsUniverseBuilder
from src.data_preparation import StructureSettings
from src.base_structure_classes import Points, AlLatticeType, CoordinateLimits
from src.projects.carbon_honeycomb_actions import CarbonHoneycombChannel

from ..structure_operations import StructureTranslator
from .by_variance import IntercalatedChannelBuilderByVariance
from .based_on_planes_configs import IntercalatedChannelBuilderBasedOnPlaneConfigs


logger = Logger("IntercalatedChannelBuilder")


class IntercalatedChannelBuilder:
    @staticmethod
    def build_carbon_coordinates(structure_folder: str) -> Points:
        path_to_pdb_file: Path = PathBuilder.build_path_to_result_data_file(structure_folder)

        return AtomsUniverseBuilder.builds_atoms_coordinates(path_to_pdb_file)

    @staticmethod
    def build_al_coordinates_for_cell(
            to_translate_al: bool,
            al_file: str,
            coordinate_limits: CoordinateLimits = CoordinateLimits(),
    ) -> Points:
        path_to_al_pdb_file: Path = PathBuilder.build_path_to_init_data_file(file=al_file)
        coordinates_al: Points = AtomsUniverseBuilder.builds_atoms_coordinates(path_to_al_pdb_file)

        if to_translate_al:
            return StructureTranslator.translate_cell(
                cell_coordinates=coordinates_al,
                translation_limits=coordinate_limits)

        return coordinates_al

    @staticmethod
    def build_al_coordinates_for_close_packed(
            al_lattice_type: AlLatticeType,
            coordinate_limits: CoordinateLimits,
            # to_translate_al: bool,  # TODO
    ) -> Points:

        if al_lattice_type.is_fcc:
            return AtomsUniverseBuilder.build_fcc_lattice_type(
                lattice_parameter=Constants.phys.AL_LATTICE_PARAM,
                coordinate_limits=coordinate_limits,
            )

        if al_lattice_type.is_hcp:
            return AtomsUniverseBuilder.build_hcp_lattice_type(
                lattice_parameter=Constants.phys.AL_LATTICE_PARAM,
                coordinate_limits=coordinate_limits,
            )

        logger.warning("No Al atoms built for closed packed structure.")
        return Points(points=np.array([]))

    @classmethod
    def build_al_in_carbon_by_variance(
        cls,
        carbon_channel: CarbonHoneycombChannel,
        coordinates_al: Points,
        structure_settings: StructureSettings,
        to_filter_al_atoms: bool = True,
        equidistant_al_points: bool = True
    ) -> Points:
        return IntercalatedChannelBuilderByVariance.build_al_in_carbon(
            carbon_channel,
            coordinates_al,
            structure_settings,
            to_filter_al_atoms,
            equidistant_al_points,
        )

    @classmethod
    def build_al_in_carbon(
        cls,
        carbon_channel: CarbonHoneycombChannel,
        coordinates_al: Points,
        structure_settings: StructureSettings,
        to_filter_al_atoms: bool = True,
        equidistant_al_points: bool = True
    ) -> Points:
        # return IntercalatedChannelBuilderByVariance.build_al_in_carbon(
        #     carbon_channel,
        #     coordinates_al,
        #     structure_settings,
        #     to_filter_al_atoms,
        #     equidistant_al_points,
        # )
        return IntercalatedChannelBuilderBasedOnPlaneConfigs.build_al_in_carbon(
            carbon_channel,
            structure_settings,
            equidistant_al_points,
        )
