import numpy as np
from scipy.spatial.distance import cdist

from src.utils import Logger, Constants, execution_time_logger
from src.data_preparation import StructureSettings
from src.base_structure_classes import Points
from src.projects.carbon_honeycomb_actions import CarbonHoneycombChannel
from src.coordinate_operations import DistanceMeasure

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

        coordinates_al = cls._reorganize_al_atoms(coordinates_al, carbon_channel)

        DistanceMeasure.calculate_min_distances_between_points(coordinates_al.points)

        return coordinates_al

    @staticmethod
    @execution_time_logger
    def _reorganize_al_atoms(coordinates_al: Points, carbon_channel: CarbonHoneycombChannel) -> Points:
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
        min_dist: float = dist_between_al_atoms * 1.15
        step: float = dist_between_al_atoms * 0.025

        max_neighbours_start: int = -1
        max_neighbours: int = -1

        # Bool to rotate between increasing min_dist and decreasing max_neighbours
        # to_change_step: bool = False

        while True:
            try:
                logger.info(f"Min dist: {min_dist}; max neighbours: {max_neighbours}")

                result: tuple[Points, int] = AtomsFilter.remove_some_close_atoms(
                    coordinates_al, min_dist, max_neighbours)

                coordinates_al_upd: Points = result[0]
                max_neighbours = result[1]
                logger.info(f"Number of atoms after filtering to reorganize: {len(coordinates_al_upd.points)}")

                if max_neighbours > max_neighbours_start:
                    max_neighbours_start = max_neighbours

                # For debugging purposes
                # if len(coordinates_al_upd) <= 30:
                #     pass

                coordinates_al = AtomConfigurator.space_al_atoms_equidistantly(coordinates_al_upd, carbon_channel)

                # logger.info(f"Number of atoms after reorganizing: {len(coordinates_al.points)}")
                return coordinates_al

            except NotOptimalStructureError as e:
                logger.warning("Not optimal structure:", e)

                if min_dist < dist_between_al_atoms:
                    raise AlBuildError("Min dist is less than the distance between atoms, Al building failed.")

                if max_neighbours <= max_neighbours_start // 2:
                    # Reset max_neighbours
                    max_neighbours = -1
                    max_neighbours_start = -1
                    min_dist -= step
                else:
                    max_neighbours -= 1
