
import numpy as np
from numpy import ndarray, floating
from scipy.spatial.distance import cdist

class DistanceMeasure:
    @staticmethod
    def calculate_distance_from_plane(points: ndarray, A: float, B: float, C: float, D: float) -> float:
        # Calculate the distance of each point from the plane
        numerator: ndarray = np.abs(A * points[:, 0] + B * points[:, 1] + C * points[:, 2] + D)
        denominator = np.sqrt(A**2 + B**2 + C**2)
        return numerator / denominator

    @staticmethod
    def calculate_signed_distance_from_plane(points: ndarray, A: float, B: float, C: float, D: float) -> float:
        # Calculate the signed distance of each point from the plane
        numerator: ndarray = A * points[:, 0] + B * points[:, 1] + C * points[:, 2] + D
        denominator = np.sqrt(A**2 + B**2 + C**2)
        return numerator / denominator

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
