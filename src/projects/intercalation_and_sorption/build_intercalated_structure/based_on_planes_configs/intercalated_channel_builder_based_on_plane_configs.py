import numpy as np

from src.utils import Logger, Constants, execution_time_logger
from src.data_preparation import StructureSettings
from src.base_structure_classes import Points
from src.projects.carbon_honeycomb_actions import CarbonHoneycombChannel
from src.structure_visualizer import StructureVisualizer
from src.coordinate_operations import DistanceMeasure

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

        # StructureVisualizer.show_two_structures(
        #     coordinates_first=carbon_channel.points, coordinates_second=coordinates_al.points, to_build_bonds=True)

        # coordinates_al = AtomConfigurator.reorganize_al_atoms(coordinates_al, carbon_channel)
        cls._print_statistics(coordinates_al, carbon_channel)

        return coordinates_al

    @staticmethod
    def _print_statistics(
        coordinates_al: Points,
        carbon_channel: CarbonHoneycombChannel,
    ):

        al: np.ndarray = coordinates_al.points
        carbon: np.ndarray = carbon_channel.points

        min_between_al: np.float64 = np.min(DistanceMeasure.calculate_min_distances_between_points(al))
        max_between_al: np.float64 = np.max(DistanceMeasure.calculate_min_distances_between_points(al))
        mean_between_al: np.float64 = np.mean(DistanceMeasure.calculate_min_distances_between_points(al))

        min_between_c_and_al: np.float64 = np.min(DistanceMeasure.calculate_min_distances(al, carbon))
        max_between_c_and_al: np.float64 = np.max(DistanceMeasure.calculate_min_distances(al, carbon))
        mean_between_c_and_al: np.float64 = np.mean(DistanceMeasure.calculate_min_distances(al, carbon))

        logger.info(
            "Final stats:\n",
            f"min_between_al: {round(min_between_al, 3)}\n",
            f"max_between_al: {round(max_between_al, 3)}\n",
            # f"mean_between_al: {round(mean_between_al, 3)}\n",
            f"min_between_c_and_al: {round(min_between_c_and_al, 3)}\n",
            f"max_between_c_and_al: {round(max_between_c_and_al, 3)}\n",
            # f"mean_between_c_and_al: {round(mean_between_c_and_al, 3)}\n",
        )
