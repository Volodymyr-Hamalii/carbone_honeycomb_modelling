import numpy as np
from numpy import ndarray, floating
from scipy.spatial.distance import cdist

from ...coordinate_operations.distance_measurer import DistanceMeasure


class VarianceCalculator:
    # @classmethod
    # def _calculate_rotation_variance(
    #         cls, angles: ndarray, inner_points: ndarray, channel_points: ndarray) -> floating | float:
    #     """
    #     Objective function for the optimizer to minimize. It rotates the inner points using the angles and then
    #     calculates the variance of the minimum distances.
    #     """

        # Calculate variance
        # variance_related_channel: floating = cls.calculate_variance_related_channel(
        #     inner_points=rotated_points, channel_points=channel_points
        # )
        # variance_xy: floating = cls.calculate_xy_variance(rotated_points)
        # return min(variance_related_channel, variance_xy)
        # return -DistanceMeasure.calculate_min_distance_sum(rotated_points, channel_points)

    @classmethod
    def _calculate_distance_variance(
            cls, translation_vector: ndarray, channel_points: ndarray, inner_points: ndarray) -> floating | float:
        """
        Calculate the variance of the minimum distances between inner points
        and channel points after applying a translation.
        """

        # Apply translation to inner points
        translated_inner_points: ndarray = inner_points.copy()
        translated_inner_points[:, 0] += translation_vector[0]  # Along Ox
        translated_inner_points[:, 1] += translation_vector[1]  # Along Oy

        # Calculate distances between each translated inner point and all channel points
        # distances: ndarray = cdist(translated_inner_points, channel_points)

        # # Get minimum distance from each inner point to any channel point
        # min_distances: floating = np.min(distances, axis=1)

        # # Calculate the variance of these minimum distances
        # variance: floating = np.var(min_distances)

        # TO CHECK ValueError: The user-provided objective function must return a scalar value.
        return -DistanceMeasure.calculate_min_distance_sum(translated_inner_points, channel_points)

    @staticmethod
    def calculate_variance_related_channel(
        inner_points: ndarray, channel_points: ndarray
    ) -> floating:
        """ Calculate variance of the minimum distances after translation and rotation. """

        distances: ndarray = cdist(inner_points, channel_points)
        min_distances: ndarray = np.min(distances, axis=1)
        variance: floating = np.var(min_distances)
        return variance

    @staticmethod
    def calculate_xy_variance(points: ndarray) -> floating:
        """ Calculate variance of the x and y coordinates. """
        return np.var(points[:, 0]) + np.var(points[:, 1])
