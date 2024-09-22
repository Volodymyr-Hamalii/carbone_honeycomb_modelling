import math

import numpy as np
from numpy import ndarray, floating
from scipy.optimize import minimize
from scipy.spatial.distance import cdist

from src.utils import Logger

from .structure_rotator import StructureRotator
from ..structure_visualizer import StructureVisualizer


logger = Logger(__name__)


class PointsOrganizer:
    @staticmethod
    def _calculate_distance_variance(
            translation_vector: ndarray, channel_points: ndarray, inner_points: ndarray) -> floating:
        """
        Calculate the variance of the minimum distances between inner points
        and channel points after applying a translation.
        """

        # Apply translation to inner points
        translated_inner_points = inner_points + translation_vector

        # Calculate distances between each translated inner point and all channel points
        distances: ndarray = cdist(translated_inner_points, channel_points)

        # Get minimum distance from each inner point to any channel point
        min_distances: floating = np.min(distances, axis=1)

        # Calculate the variance of these minimum distances
        variance: floating = np.var(min_distances)

        return variance

    @staticmethod
    def _calculate_total_variance(
        translated_and_rotated_inner_points: ndarray, channel_points: ndarray
    ) -> floating:
        """
        Calculate variance of the minimum distances after translation and rotation.
        """
        distances = cdist(translated_and_rotated_inner_points, channel_points)
        min_distances = np.min(distances, axis=1)
        variance = np.var(min_distances)
        return variance

    @staticmethod
    def rotate_to_find_min_variance(channel_points: ndarray, inner_points: ndarray):
        """
        Rotate the points inside the channel (from inner_points set) to find equilibrium positions,
        i.e., maximally equidistant from the channel atoms.
        """
        min_variance = float("inf")  # To track the minimal variance found
        best_inner_points = inner_points.copy()

        # Define angle ranges to rotate over
        angle_range_to_rotate = np.arange(0, math.pi / 8, math.pi / 25)

        for angle_x in angle_range_to_rotate:
            # Rotate inner points around the x-axis
            x_rotated_points: ndarray = StructureRotator.rotate_on_angle_related_center(
                inner_points.copy(), angle_x=angle_x
            )

            for angle_y in angle_range_to_rotate:
                # Rotate inner points around the y-axis after the x-axis rotation
                xy_rotated_points: ndarray = StructureRotator.rotate_on_angle_related_center(
                    x_rotated_points.copy(), angle_y=angle_y
                )

                # Calculate the variance after this rotation
                variance = PointsOrganizer._calculate_total_variance(
                    translated_and_rotated_inner_points=xy_rotated_points, channel_points=channel_points
                )

                # Check if this is the best configuration so far
                if variance < min_variance:
                    min_variance = variance
                    best_inner_points: ndarray = xy_rotated_points

                    print("min_variance", min_variance)
                    StructureVisualizer.show_structure(coordinates=xy_rotated_points)

        return best_inner_points

    @classmethod
    def equidistant_points_sets_in_channel(
            cls, channel_points: ndarray, inner_points: ndarray, to_rotate: bool = True) -> ndarray:
        """
        Move the points inside the channel (from inner_points set) to occupy equilibrium positions,
        i.e., maximally equidistant from the channel atoms.

        If to_rotate=True -> rotate init points to have the min deviations.

        Returns updated inner_points.
        """
        # Initial translation vector (starting with no translation)
        initial_translation: ndarray = np.array([0.0, 0.0, 0.0])

        # Use optimization to find the best translation that minimizes the variance
        result = minimize(
            cls._calculate_distance_variance,
            initial_translation,
            args=(channel_points, inner_points),
            method="BFGS",
            options={"disp": True}
        )

        # Optimal translation vector found
        optimal_translation = result.x

        # Apply the optimal translation to the inner points
        optimized_inner_points: ndarray = inner_points + optimal_translation

        if to_rotate:
            # Return the best inner points configuration found (both translated and rotated)
            return cls.rotate_to_find_min_variance(channel_points, inner_points)

        return optimized_inner_points
