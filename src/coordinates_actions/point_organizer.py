import numpy as np
from numpy import ndarray
from scipy.spatial.distance import cdist
from scipy.optimize import minimize


class PointsOrganizer:
    @staticmethod
    def _objective_function(point: ndarray, channel_points: ndarray) -> float:
        """ Objective function to maximize the minimum distance from the point to the channel points """
        distances = cdist([point], channel_points, 'euclidean')[0]
        return -np.min(distances)  # Minimize the negative of the minimum distance

    @staticmethod
    def equidistant_points_sets_in_channel(channel_points: ndarray, inner_points: ndarray) -> tuple[ndarray, ndarray]:
        """
        Move the points inside the channel (from inner_points set) to occupy equilibrium positions,
        i.e., maximally equidistant from the channel atoms.

        Returns updated channel_points and inner_points.
        """
        # Define the number of iterations and tolerance
        max_iterations = 100
        tolerance = 1e-6

        # Initialize the result with the initial inner points
        optimized_inner_points = inner_points.copy()

        for i in range(len(optimized_inner_points)):
            point = optimized_inner_points[i]
            result = minimize(
                fun=PointsOrganizer._objective_function,
                x0=point,
                args=(channel_points,),
                method='L-BFGS-B',
                bounds=[(-np.inf, np.inf)] * point.shape[0],
                options={
                    'maxiter': max_iterations,
                    'ftol': tolerance
                }
            )

            optimized_inner_points[i] = result.x

        # Calculate the final distances from the optimized inner points to channel points
        distances_to_channel = cdist(optimized_inner_points, channel_points, 'euclidean').min(axis=1)

        return optimized_inner_points, distances_to_channel
