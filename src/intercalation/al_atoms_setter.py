import math

import numpy as np
from numpy import ndarray, floating
from scipy.optimize import minimize

from src.utils import Logger, execution_time_logger

from ..coordinates_actions import PointsOrganizer
from ..calculators import DistanceCalculator, VarianceCalculator
from ..data_preparation import StructureSettings
# from ..structure_visualizer import StructureVisualizer

from .al_atoms_filter import AlAtomsFilter


logger = Logger("equidistant_points_sets_in_channel")


class AlAtomsSetter:
    @classmethod
    def equidistant_points_sets_in_channel(
            cls,
            channel_points: ndarray,
            inner_points: ndarray,
            structure_settings: StructureSettings,
            to_rotate: bool = True,
    ) -> ndarray:
        """
        Move the points inside the channel (from inner_points set) to occupy equilibrium positions,
        i.e., maximally equidistant from the channel atoms.

        If to_rotate=True -> rotate init points to have the min deviations.

        Returns updated inner_points.
        """

        if inner_points.size == 0:
            logger.warning("No points provided to equidistant_points_sets_in_channel.")
            return inner_points

        return cls.rotate_and_translate_points(
            channel_points=channel_points,
            inner_points=inner_points,
            structure_settings=structure_settings
        )

    @classmethod
    def rotate_and_translate_points(
        cls, channel_points: ndarray, inner_points: ndarray, structure_settings: StructureSettings
    ) -> ndarray:
        # translation coordinates (X and Y) and rotation angles (along Ox and Oy)
        initial_vectors: ndarray = np.array([0.0, 0.0, 0.0, 0.0])

        PointsOrganizer.align_inner_points_along_channel_oz(
            channel_points=channel_points, inner_points=inner_points
        )

        # Use optimization to find the best translation that minimizes the variance
        result = minimize(
            cls.calculate_distance_and_rotation_variance,
            initial_vectors,
            args=(inner_points, channel_points, structure_settings),
            method="BFGS",
            options={"disp": True}
        )

        # Optimal translation vector found
        optimal_vectors: ndarray = result.x

        moved_points: ndarray = PointsOrganizer.move_and_rotate_related_xy(
            vectors=optimal_vectors, points=inner_points)

        return moved_points

    @classmethod
    def calculate_distance_and_rotation_variance(
            cls,
            initial_vectors: ndarray,
            inner_points: ndarray,
            channel_points: ndarray,
            structure_settings: StructureSettings
    ) -> float | floating:

        num_of_atoms: int = len(inner_points)

        moved_points: ndarray = PointsOrganizer.move_and_rotate_related_xy(
            vectors=initial_vectors, points=inner_points)

        filtered_atoms: ndarray = AlAtomsFilter.filter_al_atoms_related_carbon(
            coordinates_al=inner_points,
            coordinates_carbon=channel_points,
            structure_settings=structure_settings)
        num_of_atoms_after_filter: int = len(filtered_atoms)

        if num_of_atoms_after_filter < num_of_atoms:
            return np.inf

        min_distance_sum: floating = DistanceCalculator.calculate_min_distance_sum(moved_points, channel_points)
        variance_related_channel: floating = VarianceCalculator.calculate_variance_related_channel(
            moved_points, channel_points)
        distance_and_rotation_variance: floating = variance_related_channel * len(inner_points) - min_distance_sum
        return distance_and_rotation_variance
