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
            translation_vector: ndarray, channel_points: ndarray, inner_points: ndarray) -> floating | float:
        """
        Calculate the variance of the minimum distances between inner points
        and channel points after applying a translation.
        """

        # Apply translation to inner points
        translated_inner_points: ndarray = inner_points.copy()
        translated_inner_points[:, 0] += translation_vector[0]  # Along Ox
        translated_inner_points[:, 1] += translation_vector[1]  # Along Oy

        # Calculate distances between each translated inner point and all channel points
        distances: ndarray = cdist(translated_inner_points, channel_points)

        # Get minimum distance from each inner point to any channel point
        min_distances: floating = np.min(distances, axis=1)

        # Calculate the variance of these minimum distances
        variance: floating = np.var(min_distances)

        return variance

    @staticmethod
    def _calculate_rotation_variance(angles: ndarray, inner_points: ndarray, channel_points: ndarray) -> floating:
        """
        Objective function for the optimizer to minimize. It rotates the inner points using the angles and then
        calculates the variance of the minimum distances.
        """
        angle_x, angle_y = angles

        # Rotate inner points by angle_x and angle_y
        rotated_points_x: ndarray = StructureRotator.rotate_on_angle_related_center(
            inner_points.copy(), angle_x=angle_x
        )
        rotated_points_xy: ndarray = StructureRotator.rotate_on_angle_related_center(
            rotated_points_x.copy(), angle_y=angle_y
        )

        if (np.min(rotated_points_xy[:, 2]) < np.min(channel_points[:, 2])) or (
                np.max(rotated_points_xy[:, 2]) > np.max(channel_points[:, 2])):
            # Filter the cases when inner_points are out of the channel
            return np.inf

        # Calculate variance
        # variance: floating = PointsOrganizer._calculate_total_variance(
        #     translated_and_rotated_inner_points=rotated_points_xy, channel_points=channel_points
        # )
        variance: floating = np.var(rotated_points_xy[:, 0]) + np.var(rotated_points_xy[:, 1])

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

        cls.align_inner_points_along_channel_oz(
            channel_points=channel_points, inner_points=inner_points)

        # Initial guess for angles (both x and y) in radians
        initial_angles = np.array([0.0, 0.0])

        # Minimize the objective function with respect to the rotation angles
        result = minimize(
            PointsOrganizer._calculate_rotation_variance,
            initial_angles,
            args=(inner_points, channel_points),
            method="BFGS",
            options={"disp": True}
        )

        # Extract the optimal angles
        optimal_angle_x, optimal_angle_y = result.x

        # Apply the optimal rotations to the inner points
        rotated_points_x = StructureRotator.rotate_on_angle_related_center(
            inner_points.copy(), angle_x=optimal_angle_x
        )
        rotated_points_xy = StructureRotator.rotate_on_angle_related_center(
            rotated_points_x.copy(), angle_y=optimal_angle_y
        )

        return rotated_points_xy

    @classmethod
    def equidistant_points_sets_in_channel(
            cls, channel_points: ndarray, inner_points: ndarray, to_rotate: bool = True) -> ndarray:
        """
        Move the points inside the channel (from inner_points set) to occupy equilibrium positions,
        i.e., maximally equidistant from the channel atoms.

        If to_rotate=True -> rotate init points to have the min deviations.

        Returns updated inner_points.
        """

        if to_rotate:
            # Return the best inner points configuration found (both translated and rotated)
            inner_points = cls.rotate_to_find_min_variance(channel_points, inner_points)

        # Initial translation vector for x and y (starting with no translation)
        initial_translation: ndarray = np.array([0.0, 0.0])

        # Use optimization to find the best translation that minimizes the variance
        result = minimize(
            cls._calculate_distance_variance,
            initial_translation,
            args=(channel_points, inner_points),
            method="BFGS",
            options={"disp": True}
        )

        # Optimal translation vector found
        optimal_translation: ndarray = result.x

        # Apply the optimal translation to the inner points
        inner_points[:, 0] += optimal_translation[0]
        inner_points[:, 1] += optimal_translation[1]

        return inner_points
