from math import sqrt
import numpy as np
from numpy import ndarray, floating
from scipy.spatial.distance import cdist


class DistanceMeasure:
    @staticmethod
    def calculate_distance_between_2_points(p1: tuple | np.ndarray, p2: tuple | np.ndarray) -> float:
        """
        Compute distance between two points (x1, y1, z1) and (x2, y2, z2),
        where z1 and z2 are optional.
        """
        # Add z coordinate if it doesn't exist
        if len(p1) == 2:
            p1 = (p1[0], p1[1], 0)

        if len(p2) == 2:
            p2 = (p2[0], p2[1], 0)

        # Calculate the distance
        return sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

    @staticmethod
    def calculate_distance_from_plane(points: ndarray, A: float, B: float, C: float, D: float) -> float:
        """ Calculate the distance of each point from the plane. """
        numerator: ndarray = np.abs(A * points[:, 0] + B * points[:, 1] + C * points[:, 2] + D)
        denominator = np.sqrt(A**2 + B**2 + C**2)
        return numerator / denominator

    @staticmethod
    def calculate_signed_distance_from_plane(points: ndarray, A: float, B: float, C: float, D: float) -> float:
        """ Calculate the signed distance of each point from the plane. """
        numerator: ndarray = A * points[:, 0] + B * points[:, 1] + C * points[:, 2] + D
        denominator = np.sqrt(A**2 + B**2 + C**2)
        return numerator / denominator

    @staticmethod
    def calculate_dist_matrix(points: ndarray) -> ndarray:
        """ Calculate distance matrix between provided points (with inf in diagonal). """
        inf_diag_matrix: ndarray = np.diag([np.inf] * len(points))

        # Add inf_diag_matrix to distances tp remove zeros
        return cdist(points, points) + inf_diag_matrix

    @staticmethod
    def calculate_min_distances(points_1: ndarray, points_2: ndarray) -> ndarray:
        """ Returns min distance between 2 provided point sets. """
        distances: ndarray = cdist(points_1, points_2)
        return np.min(distances, axis=1)

    @classmethod
    def calculate_min_distance_sum(cls, points_1: ndarray, points_2: ndarray) -> floating:
        """ Returns sum of min distances between 2 provided point sets. """
        min_distances: ndarray = cls.calculate_min_distances(points_1, points_2)
        return np.sum(min_distances)

    @classmethod
    def calculate_min_distances_between_points(cls, points: ndarray) -> ndarray:
        distances: ndarray = cls.calculate_dist_matrix(points)
        return np.min(distances, axis=1)
