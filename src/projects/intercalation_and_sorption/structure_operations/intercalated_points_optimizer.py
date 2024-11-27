import math
from typing import Callable
import numpy as np
from numpy import ndarray, floating
from scipy.optimize import minimize
from scipy.spatial.distance import cdist

from src.utils import Logger, execution_time_logger
from src.coordinate_operations import (
    PointsRotator,
)
from ..build_intercalated_structure.by_variance.variance_calculator import VarianceCalculator
from ..intercalated_coordinates_utils import IntercalatedCoordinatesUtils


class IntercalatedPointsOptimizer:
    """ To find equilibrium positions of the intercalated (or sorbed) points. """

    @classmethod
    def rotate_to_find_min_variance_by_minimize(
        cls,
        inner_points: ndarray,
        variance_function: Callable,
        variance_function_args: tuple | None = None,
        minimize_method: str = "BFGS",  # TODO: check the method
    ) -> ndarray:
        """
        Rotate the points inside the channel (from inner_points set) to find equilibrium positions,
        i.e., maximally equidistant from the channel atoms.
        """

        # Initial guess for angles (both x and y) in radians
        initial_angles: ndarray = np.array([0.0, 0.0])

        # Minimize the objective function with respect to the rotation angles
        result = minimize(
            variance_function,
            initial_angles,
            args=variance_function_args,
            method=minimize_method,
            options={"disp": True, "maxiter": 1000}
        )

        # Extract the optimal angles
        optimal_angle_x, optimal_angle_y = result.x

        # Apply the optimal rotations to the inner points
        return PointsRotator.rotate_on_angle_related_center(
            inner_points.copy(), angle_x=optimal_angle_x, angle_y=optimal_angle_y
        )

    @classmethod
    @execution_time_logger
    def rotate_to_find_min_variance_by_iterations(cls, channel_points: ndarray, inner_points: ndarray) -> ndarray:
        """
        Rotate the points inside the channel (from inner_points set) to find equilibrium positions,
        i.e., maximally equidistant from the channel atoms.
        """

        min_variance: float | np.floating = np.inf  # To track the minimal variance found
        best_inner_points = inner_points.copy()

        # Define angle ranges to rotate over
        angle_range_to_rotate = np.arange(- math.pi / 16, math.pi / 16, math.pi / 64)

        for angle_x in angle_range_to_rotate:
            # Rotate inner points around the x-axis
            x_rotated_points: ndarray = PointsRotator.rotate_on_angle_related_center(
                inner_points.copy(), angle_x=angle_x)

            for angle_y in angle_range_to_rotate:
                # Rotate inner points around the y-axis after the x-axis rotation
                xy_rotated_points: ndarray = PointsRotator.rotate_on_angle_related_center(
                    x_rotated_points.copy(), angle_y=angle_y)

                IntercalatedCoordinatesUtils.align_inner_points_along_channel_oz(
                    channel_points=channel_points, intercaleted_points=xy_rotated_points)

                if min(inner_points[:, 2]) < min(channel_points[:, 2]):
                    # Skip if some inner_points are out the channel along Oz
                    continue

                # Calculate the variance after this rotation
                variance: floating = VarianceCalculator.calculate_variance_related_channel(
                    channel_points=channel_points, inner_points=xy_rotated_points)

                # Check if this is the best configuration so far
                if variance < min_variance:
                    min_variance = variance
                    best_inner_points: ndarray = xy_rotated_points

        return best_inner_points
