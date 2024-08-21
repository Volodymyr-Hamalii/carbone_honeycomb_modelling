import numpy as np
from numpy import ndarray, floating
from scipy.optimize import minimize
from scipy.spatial.distance import cdist


class PointsOrganizer:
    @staticmethod
    def _calculate_distance_variance(
            translation_vector: ndarray, channel_points: ndarray, inner_points: ndarray) -> floating:
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
    def equidistant_points_sets_in_channel(channel_points: ndarray, inner_points: ndarray) -> ndarray:
        """
        Move the points inside the channel (from inner_points set) to occupy equilibrium positions,
        i.e., maximally equidistant from the channel atoms.

        Returns updated inner_points.
        """
        # Initial translation vector (starting with no translation)
        initial_translation: ndarray = np.array([0.0, 0.0, 0.0])

        # Use optimization to find the best translation that minimizes the variance
        result = minimize(
            PointsOrganizer._calculate_distance_variance,
            initial_translation,
            args=(channel_points, inner_points),
            method="BFGS",
            options={"disp": True}
        )

        # Optimal translation vector found
        optimal_translation = result.x

        # Apply the optimal translation to the inner points
        optimized_inner_points = inner_points + optimal_translation

        return optimized_inner_points
