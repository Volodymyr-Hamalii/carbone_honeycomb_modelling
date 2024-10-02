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

    @classmethod
    def _calculate_rotation_variance(
            cls, angles: ndarray, inner_points: ndarray, channel_points: ndarray) -> floating | float:
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
        variance_related_channel: floating = cls.calculate_variance_related_channel(
            inner_points=rotated_points_xy, channel_points=channel_points
        )
        variance_xy: floating = cls.calculate_xy_variance(rotated_points_xy)
        return min(variance_related_channel, variance_xy)

    @staticmethod
    def calculate_variance_related_channel(
        inner_points: ndarray, channel_points: ndarray
    ) -> floating:
        """ Calculate variance of the minimum distances after translation and rotation. """

        distances = cdist(inner_points, channel_points)
        min_distances = np.min(distances, axis=1)
        variance = np.var(min_distances)
        return variance

    @staticmethod
    def calculate_xy_variance(points: ndarray) -> floating:
        """ Calculate variance of the x and y coordinates. """
        return np.var(points[:, 0]) + np.var(points[:, 1])

    @classmethod
    def rotate_to_find_min_variance_by_minimize(cls, channel_points: ndarray, inner_points: ndarray) -> ndarray:
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
            cls._calculate_rotation_variance,
            initial_angles,
            args=(inner_points, channel_points),
            method="BFGS",
            options={"disp": True, "maxiter": 1000}
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
    def rotate_to_find_min_variance(cls, channel_points: ndarray, inner_points: ndarray) -> ndarray:
        """
        Rotate the points inside the channel (from inner_points set) to find equilibrium positions,
        i.e., maximally equidistant from the channel atoms.
        """

        cls.align_inner_points_along_channel_oz(
            channel_points=channel_points, inner_points=inner_points)

        min_variance: float | np.floating = np.inf  # To track the minimal variance found
        best_inner_points = inner_points.copy()

        # Define angle ranges to rotate over
        angle_range_to_rotate = np.arange(- math.pi / 16, math.pi / 16, math.pi / 64)

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

                cls.align_inner_points_along_channel_oz(
                    channel_points=channel_points, inner_points=inner_points)

                if min(inner_points[:, 2]) < min(channel_points[:, 2]):
                    # Skip if we are out the channel limits
                    continue

                # Calculate the variance after this rotation
                variance = cls.calculate_variance_related_channel(
                    inner_points=xy_rotated_points, channel_points=channel_points)

                # Check if this is the best configuration so far
                if variance < min_variance:
                    min_variance = variance
                    best_inner_points: ndarray = xy_rotated_points

        return best_inner_points

    @staticmethod
    def align_inner_points_along_channel_oz(channel_points: ndarray, inner_points: ndarray) -> None:
        """ Align the points in the middle relative to the Oz channel. """

        min_channel, max_channel = np.min(channel_points[:, 2]), np.max(channel_points[:, 2])
        min_inner, max_inner = np.min(inner_points[:, 2]), np.max(inner_points[:, 2])

        mid_channel: np.float64 = (max_channel + min_channel) / 2
        mid_inner: np.float64 = (max_inner + min_inner) / 2

        move_z_to: np.float64 = mid_channel - mid_inner
        inner_points[:, 2] += move_z_to

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
