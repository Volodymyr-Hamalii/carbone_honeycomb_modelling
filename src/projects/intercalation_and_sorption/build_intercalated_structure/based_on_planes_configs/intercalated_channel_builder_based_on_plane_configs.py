from typing import Any, Generator
import numpy as np

from src.utils import Logger, Constants, execution_time_logger
from src.data_preparation import StructureSettings
from src.base_structure_classes import Points
from src.projects.carbon_honeycomb_actions import CarbonHoneycombChannel
from src.coordinate_operations import DistanceMeasure
from src.structure_visualizer import StructureVisualizer

from .atoms_builder import AtomsBuilder
from .atoms_filter import AtomsFilter
from .atom_configurator import AtomConfigurator
from .intercalation_exceptions import NotOptimalStructureError, AlBuildError


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

        # StructureVisualizer.show_two_structures(
        #     coordinates_first=carbon_channel.points, coordinates_second=coordinates_al.points, to_build_bonds=True)

        coordinates_al = AtomConfigurator.reorganize_al_atoms(coordinates_al, carbon_channel)
        cls._print_statistics(coordinates_al, carbon_channel)

        return coordinates_al

    @staticmethod
    def _iterate_params_for_reorganizing_al_atoms() -> Generator[tuple[float, int], Any, None]:
        # max_points_to_move_before_reset_coef: list[float] = [1.1, 1.3, 1.5, 1.7, 2]
        max_points_to_move_before_reset_coef: list[float] = [1, 0.9, 0.75, 0.5]
        num_of_points_to_skip_list: list[int] = [0, 1, 2, 3]

        for coef in max_points_to_move_before_reset_coef:
            for num_of_points_to_skip in num_of_points_to_skip_list:
                yield (coef, num_of_points_to_skip)

    @classmethod
    @execution_time_logger
    def _reorganize_al_atoms(cls, coordinates_al: Points, carbon_channel: CarbonHoneycombChannel) -> Points:
        """
        Move Al atoms to have the correct minimal distance between them
        (or remove them if we can't achive optimal configuration just by moving).
        """

        logger.info(f"Number of atoms before reorganizing: {len(coordinates_al.points)}")

        # # Take a step as 1% of the distance between atoms
        # step: float = Constants.phys.al.DIST_BETWEEN_ATOMS * 0.05

        # start: float = Constants.phys.al.DIST_BETWEEN_ATOMS * 1.2
        # end: float = Constants.phys.al.DIST_BETWEEN_ATOMS * 1.21

        # best_result: Points = coordinates_al
        # best_result_len: int = 0

        # for min_dist in np.arange(start, end, step):
        #     logger.info(f"Min dist: {min_dist}")

        #     coordinates_al_upd: Points = AtomsFilter.remove_some_close_atoms(coordinates_al, min_dist)

        #     if len(coordinates_al_upd.points) > best_result_len:
        #         best_result = coordinates_al_upd
        #         best_result_len = len(coordinates_al_upd.points)

        #     elif len(coordinates_al_upd.points) == len(coordinates_al.points):
        #         break

        # Filter and move Al atoms to have the correct minimal distance between them
        dist_between_al_atoms: float = Constants.phys.al.DIST_BETWEEN_ATOMS
        min_dist: float = dist_between_al_atoms * 1.2
        step: float = dist_between_al_atoms * 0.025

        max_neighbours_start: int = -1
        max_neighbours: int = -1

        # Bool to rotate between increasing min_dist and decreasing max_neighbours
        # to_change_step: bool = False

        best_result: Points = Points(points=np.array([]))
        min_distances: np.ndarray = np.array([])

        for max_points_to_move_before_reset_coef, num_of_points_to_skip in cls._iterate_params_for_reorganizing_al_atoms():
            logger.info(f"Max points to move before reset coef: {max_points_to_move_before_reset_coef}")
            while True:
                try:
                    logger.info(f"Min dist: {min_dist}; max neighbours: {max_neighbours}")

                    result: tuple[Points, int] = AtomsFilter.remove_some_close_atoms(
                        coordinates_al, min_dist, max_neighbours, num_of_points_to_skip)

                    coordinates_al_upd: Points = result[0]
                    max_neighbours = result[1]
                    logger.info(f"Number of atoms after filtering to reorganize: {len(coordinates_al_upd.points)}")

                    if max_neighbours > max_neighbours_start:
                        max_neighbours_start = max_neighbours

                    # For debugging purposes
                    # if len(coordinates_al_upd) <= 30:
                    #     pass
                    max_points_to_move_before_reset: int = round(
                        len(coordinates_al_upd) * max_points_to_move_before_reset_coef)

                    coordinates_al_upd = AtomConfigurator.space_al_atoms_equidistantly(
                        coordinates_al_upd,
                        carbon_channel,
                        max_points_to_move_before_reset,
                        num_of_points_to_skip=0,
                    )

                    # logger.info(f"Number of atoms after reorganizing: {len(coordinates_al.points)}")
                    if len(coordinates_al_upd.points) > len(best_result.points):
                        best_result = coordinates_al_upd
                        min_distances = DistanceMeasure.calculate_min_distances(
                            carbon_channel.points, coordinates_al_upd.points)

                    elif len(coordinates_al.points) == len(coordinates_al_upd.points):
                        # Get the set with the minimal sum of distances to carbon atoms
                        min_distances_upd: np.ndarray = DistanceMeasure.calculate_min_distances(
                            carbon_channel.points, coordinates_al_upd.points)

                        if np.sum(min_distances_upd) < np.sum(min_distances):
                            best_result = coordinates_al_upd
                            min_distances = min_distances_upd

                except NotOptimalStructureError as e:
                    logger.warning("Not optimal structure:", e)

                if min_dist < dist_between_al_atoms:
                    # msg: str = "Min dist is less than the distance between atoms, Al building failed."
                    # raise AlBuildError(msg)
                    break

                if max_neighbours <= max_neighbours_start // 2:
                    # Reset max_neighbours
                    max_neighbours = -1
                    max_neighbours_start = -1
                    min_dist -= step
                else:
                    max_neighbours -= 1

        return best_result
