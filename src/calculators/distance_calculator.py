import numpy as np
from numpy import ndarray, floating
from scipy.spatial.distance import cdist


class DistanceCalculator:
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

    @staticmethod
    def calculate_min_distances_between_points(points: ndarray) -> ndarray:
        inf_diag_matrix: ndarray = np.diag([np.inf] * len(points))

        # Add inf_diag_matrix to distances tp remove zeros
        distances: ndarray = cdist(points, points) + inf_diag_matrix

        return np.min(distances, axis=1)
