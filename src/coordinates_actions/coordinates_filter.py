import numpy as np
from numpy import ndarray

from .planes_builder import PlanesBuilder


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
        """
        A, B, C, D are the parameters from the
        Ax + By + Cz + D = 0 plane equation.
        """

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

    @staticmethod
    def filter_by_min_max_z(
        coordinates_to_filter: ndarray,
        coordinates_with_min_max_z: ndarray,
        move_align_z: bool = False,
    ) -> ndarray:
        """Filter coordinates_to_filter by min and max z coordinate of coordinates_with_min_max_z."""

        if coordinates_to_filter.size == 0:
            return coordinates_to_filter

        z_min: np.float32 = np.min(coordinates_with_min_max_z[:, 2])
        z_max: np.float32 = np.max(coordinates_with_min_max_z[:, 2])

        filtered_by_min_z = coordinates_to_filter[coordinates_to_filter[:, 2] >= z_min]

        if filtered_by_min_z.size == 0:
            return filtered_by_min_z

        if move_align_z is True:
            # Move Al atoms along Oz down to align the lowest Al atom with the channel bottom
            move_to: np.float32 = np.min(filtered_by_min_z[:, 2]) - z_min
            filtered_by_min_z[:, 2] = filtered_by_min_z[:, 2] - move_to

        return coordinates_to_filter[coordinates_to_filter[:, 2] <= z_max]
