import math

import numpy as np
from numpy import ndarray, floating
from scipy.optimize import minimize

from src.utils import Logger, execution_time_logger
from src.data_preparation import StructureSettings
from src.coordinate_operations import (
    DistanceMeasure,
    PointsMover,
    PointsRotator,
)
from ...intercalated_coordinates_utils import IntercalatedCoordinatesUtils
from .variance_calculator import VarianceCalculator
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

        IntercalatedCoordinatesUtils.align_inner_points_along_channel_oz(
            channel_points=channel_points, intercaleted_points=inner_points
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

        moved_points: ndarray = cls.move_and_rotate_related_xy(
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

        moved_points: ndarray = cls.move_and_rotate_related_xy(
            vectors=initial_vectors, points=inner_points)

        filtered_atoms: ndarray = AlAtomsFilter.filter_al_atoms_related_carbon(
            coordinates_al=inner_points,
            coordinates_carbon=channel_points,
            structure_settings=structure_settings)
        num_of_atoms_after_filter: int = len(filtered_atoms)

        if num_of_atoms_after_filter < num_of_atoms:
            return np.inf

        min_distance_sum: floating = DistanceMeasure.calculate_min_distance_sum(moved_points, channel_points)
        variance_related_channel: floating = VarianceCalculator.calculate_variance_related_channel(
            moved_points, channel_points)
        distance_and_rotation_variance: floating = variance_related_channel * len(inner_points) - min_distance_sum
        return distance_and_rotation_variance

    @staticmethod
    def move_and_rotate_related_xy(vectors: ndarray, points: ndarray) -> ndarray:
        translation_x, translation_y = vectors[:2]
        move_vector = np.array([translation_x, translation_y, 0])
        angle_x, angle_y = vectors[2:]

        # Apply translation to inner points
        translated_inner_points: ndarray = points.copy()

        moved_points = PointsMover.move_on_vector(
            points=translated_inner_points,
            vector=np.array([vectors[0], vectors[1], 0]),
        )
        # Rotate inner points by angle_x and angle_y
        return PointsRotator.rotate_on_angle_related_center(
            translated_inner_points, angle_x=angle_x, angle_y=angle_y)
