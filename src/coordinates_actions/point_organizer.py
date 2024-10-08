import math

import numpy as np
from numpy import ndarray, floating
from scipy.optimize import minimize
from scipy.spatial.distance import cdist

from src.utils import Logger, execution_time_logger

from .structure_rotator import StructureRotator
from ..structure_visualizer import StructureVisualizer


logger = Logger("PointsOrganizer")


class PointsOrganizer:

    # Calculate distances

    @staticmethod
    def calculate_min_distance(points_set_1: ndarray, points_set_2: ndarray) -> floating:
        """ Returns min distance between 2 provided point sets. """
        distances: ndarray = cdist(points_set_1, points_set_2)
        return np.min(distances, axis=1)

    @staticmethod
    def calculate_min_distance_sum(points_set_1: ndarray, points_set_2: ndarray) -> floating:
        """ Returns sum of min distances between 2 provided point sets. """
        distances: ndarray = cdist(points_set_1, points_set_2)
        min_distances: ndarray = np.min(distances, axis=1)
        return np.sum(min_distances)

    # Calculate variance

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
        return -cls.calculate_min_distance_sum(translated_inner_points, channel_points)

    @classmethod
    def _calculate_rotation_variance(
            cls, angles: ndarray, inner_points: ndarray, channel_points: ndarray) -> floating | float:
        """
        Objective function for the optimizer to minimize. It rotates the inner points using the angles and then
        calculates the variance of the minimum distances.
        """
        angle_x, angle_y = angles

        # Rotate inner points by angle_x and angle_y
        rotated_points: ndarray = StructureRotator.rotate_on_angle_related_center(
            inner_points.copy(), angle_x=angle_x, angle_y=angle_y
        )

        if (np.min(rotated_points[:, 2]) < np.min(channel_points[:, 2])) or (
                np.max(rotated_points[:, 2]) > np.max(channel_points[:, 2])):
            # Filter the cases when inner_points are out of the channel
            return np.inf

        # Calculate variance
        # variance_related_channel: floating = cls.calculate_variance_related_channel(
        #     inner_points=rotated_points, channel_points=channel_points
        # )
        # variance_xy: floating = cls.calculate_xy_variance(rotated_points)
        # return min(variance_related_channel, variance_xy)
        return -cls.calculate_min_distance_sum(rotated_points, channel_points)

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

    @staticmethod
    def _move_and_rotate_related_xy(vectors: ndarray, points: ndarray) -> ndarray:
        translation_x, translation_y = vectors[:2]
        angle_x, angle_y = vectors[2:]

        # Apply translation to inner points
        translated_inner_points: ndarray = points.copy()
        translated_inner_points[:, 0] += translation_x  # Along Ox
        translated_inner_points[:, 1] += translation_y  # Along Oy

        # Rotate inner points by angle_x and angle_y
        return StructureRotator.rotate_on_angle_related_center(
            translated_inner_points, angle_x=angle_x, angle_y=angle_y)

    @classmethod
    def _calculate_distance_and_rotation_variance(
        cls, initial_vectors: ndarray, inner_points: ndarray, channel_points: ndarray
    ) -> floating:

        moved_points: ndarray = cls._move_and_rotate_related_xy(
            vectors=initial_vectors, points=inner_points)

        min_distance_sum: floating = cls.calculate_min_distance_sum(moved_points, channel_points)
        variance_related_channel: floating = cls.calculate_variance_related_channel(
            moved_points, channel_points) * len(inner_points)
        distance_and_rotation_variance: floating = variance_related_channel - min_distance_sum
        return distance_and_rotation_variance

    # Rotations

    @classmethod
    def rotate_to_find_min_variance_by_minimize(cls, channel_points: ndarray, inner_points: ndarray) -> ndarray:
        """
        Rotate the points inside the channel (from inner_points set) to find equilibrium positions,
        i.e., maximally equidistant from the channel atoms.
        """

        cls.align_inner_points_along_channel_oz(
            channel_points=channel_points, inner_points=inner_points)

        # Initial guess for angles (both x and y) in radians
        initial_angles: ndarray = np.array([0.0, 0.0])

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
        rotated_points: ndarray = StructureRotator.rotate_on_angle_related_center(
            inner_points.copy(), angle_x=optimal_angle_x, angle_y=optimal_angle_y
        )

        return rotated_points

    @classmethod
    @execution_time_logger
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
                inner_points.copy(), angle_x=angle_x)

            for angle_y in angle_range_to_rotate:
                # Rotate inner points around the y-axis after the x-axis rotation
                xy_rotated_points: ndarray = StructureRotator.rotate_on_angle_related_center(
                    x_rotated_points.copy(), angle_y=angle_y)

                cls.align_inner_points_along_channel_oz(
                    channel_points=channel_points, inner_points=xy_rotated_points)

                if min(inner_points[:, 2]) < min(channel_points[:, 2]):
                    # Skip if some inner_points are out the channel along Oz
                    continue

                # Calculate the variance after this rotation
                variance: floating = cls.calculate_variance_related_channel(
                    channel_points=channel_points, inner_points=xy_rotated_points)

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

    # Mixed operations

    @classmethod
    def rotate_and_translate_points(
        cls, channel_points: ndarray, inner_points: ndarray
    ) -> ndarray:
        # translation coordinates (X and Y) and rotation angles (along Ox and Oy)
        initial_vectors: ndarray = np.array([0.0, 0.0, 0.0, 0.0])

        # Use optimization to find the best translation that minimizes the variance
        result = minimize(
            cls._calculate_distance_and_rotation_variance,
            initial_vectors,
            args=(inner_points, channel_points),
            method="BFGS",
            options={"disp": True}
        )

        # Optimal translation vector found
        optimal_vectors: ndarray = result.x

        moved_points: ndarray = cls._move_and_rotate_related_xy(
            vectors=optimal_vectors, points=inner_points)

        return moved_points

    @classmethod
    @execution_time_logger
    def equidistant_points_sets_in_channel(
            cls, channel_points: ndarray, inner_points: ndarray, to_rotate: bool = True) -> ndarray:
        """
        Move the points inside the channel (from inner_points set) to occupy equilibrium positions,
        i.e., maximally equidistant from the channel atoms.

        If to_rotate=True -> rotate init points to have the min deviations.

        Returns updated inner_points.
        """

        if inner_points.size == 0:
            logger.warning("No points provided to PointsOrganizer.equidistant_points_sets_in_channel.")
            return inner_points

        # TODO: REFACTOR
        return cls.rotate_and_translate_points(
            channel_points=channel_points,
            inner_points=inner_points,
        )

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
