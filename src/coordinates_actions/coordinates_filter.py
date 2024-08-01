import numpy as np
from numpy import ndarray


class CoordinatesFilter:
    @staticmethod
    def _distance_from_plane(points: ndarray, A: float, B: float, C: float, D: float) -> float:
        # Calculate the distance of each point from the plane
        numerator: ndarray = np.abs(A * points[:, 0] + B * points[:, 1] + C * points[:, 2] + D)
        denominator = np.sqrt(A**2 + B**2 + C**2)
        return numerator / denominator

    @staticmethod
    def _signed_distance_from_plane(points: ndarray, A: float, B: float, C: float, D: float) -> float:
        # Calculate the signed distance of each point from the plane
        numerator: ndarray = A * points[:, 0] + B * points[:, 1] + C * points[:, 2] + D
        denominator = np.sqrt(A**2 + B**2 + C**2)
        return numerator / denominator

    @classmethod
    def filter_coordinates_related_to_plane(
            cls,
            points: ndarray,
            A: float, B: float, C: float, D: float,
            min_distance: float = 0,
            direction: int = 1
    ):
        signed_distances: float = cls._signed_distance_from_plane(points, A, B, C, D)
        if direction == 1:
            # Keep points above the plane at the minimum distance
            filtered_points = points[signed_distances >= min_distance]
        elif direction == -1:
            # Keep points below the plane at the minimum distance
            filtered_points = points[signed_distances <= -min_distance]
        else:
            raise ValueError("Direction should be 1 (above the plane) or -1 (below the plane)")
        return filtered_points
