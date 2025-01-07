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

from .al_atoms_filter import AlAtomsFilter
from .al_atoms_setter import AlAtomsSetter


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
    def build_al_in_carbon(
            cls,
            carbon_channel: CarbonHoneycombChannel,
            coordinates_al: Points,
            structure_settings: StructureSettings,
            to_filter_al_atoms: bool = True,
            equidistant_al_points: bool = True
    ) -> Points:
        """ Return coordinates_carbon, coordinates_al """

        if to_filter_al_atoms:
            coordinates_al_filtered: Points = cls.find_max_filtered_atoms(
                carbon_channel=carbon_channel,
                coordinates_al=coordinates_al,
                structure_settings=structure_settings)

            if equidistant_al_points:
                # Set Al atoms maximally equidistant from the channel atoms
                coordinates_al_filtered: Points = AlAtomsSetter.equidistant_points_sets_in_channel(
                    carbon_channel, coordinates_al_filtered, structure_settings)

            return coordinates_al_filtered

        logger.info("Build AL in the channel without Al atoms filtering.")
        return coordinates_al

    @classmethod
    @execution_time_logger
    def find_max_filtered_atoms(
            cls,
            carbon_channel: CarbonHoneycombChannel,
            coordinates_al: Points,
            structure_settings: StructureSettings,
    ) -> Points:
        """
        Parallel move, rotate and filter Al atoms related carbon atoms
        to find the option with the maximum Al atoms after filtering.
        """

        max_atoms: int = 0
        min_dist_between_al_sum: float = np.inf
        dist_and_rotation_variance: float | np.floating = 0

        coordinates_al_result: Points = Points(points=np.array([]))

        al_param: float = Constants.phys.AL_LATTICE_PARAM
        step_to_move: float = al_param / 25
        range_to_move: np.ndarray = np.arange(0, al_param, step_to_move)

        step_to_rotate: float = math.pi / 45
        angle_range_to_rotate: np.ndarray = np.arange(0, math.pi / 3, step_to_rotate)

        for step_x in range_to_move:
            moved_x_coordinates_al: Points = coordinates_al.copy()
            moved_x_coordinates_al.points[:, 0] += step_x

            for step_y in range_to_move:
                moved_xy_coordinates_al: Points = moved_x_coordinates_al.copy()
                moved_xy_coordinates_al.points[:, 1] += step_y

                for angle_x in angle_range_to_rotate:
                    x_rotaded_coordinates_al: Points = PointsRotator.rotate_on_angle_related_center(
                        moved_xy_coordinates_al.copy(), angle_x=angle_x)

                    for angle_y in angle_range_to_rotate:
                        xy_rotaded_coordinates_al: Points = PointsRotator.rotate_on_angle_related_center(
                            x_rotaded_coordinates_al.copy(), angle_y=angle_y)

                        for angle_z in angle_range_to_rotate:
                            xyz_rotaded_coordinates_al: Points = PointsRotator.rotate_on_angle_related_center(
                                xy_rotaded_coordinates_al.copy(), angle_z=angle_z)

                            result: tuple = AlAtomsFilter.get_filtered_al_atoms(
                                carbon_channel=carbon_channel,
                                coordinates_al=xyz_rotaded_coordinates_al,
                                coordinates_al_prev=coordinates_al_result,
                                structure_settings=structure_settings,
                                max_atoms=max_atoms,
                                min_dist_between_al_sum_prev=min_dist_between_al_sum,
                                dist_and_rotation_variance_prev=dist_and_rotation_variance,
                            )

                            coordinates_al_result = result[0]
                            min_dist_between_al_sum = result[1]
                            dist_and_rotation_variance = result[2]
                            max_atoms = result[3]

        return coordinates_al_result
