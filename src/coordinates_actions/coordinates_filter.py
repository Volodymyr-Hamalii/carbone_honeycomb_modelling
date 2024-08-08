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

    @classmethod
    def filter_by_max_min_z(cls, coordinates: ndarray, z_min: float, z_max: float) -> ndarray:
        """ Remove coordinates below z_min and above z_max limits """

        # Below
        min_plane_params: tuple[float, float, float, float] = PlanesBuilder.build_plane_parameters(
            p1=[0, 0, z_min], p2=[1, 1, z_min], p3=[1, 0, z_min])
        coordinates = cls.filter_coordinates_related_to_plane(
            coordinates, *min_plane_params, direction=-1
        )

        # Above
        min_plane_params: tuple[float, float, float, float] = PlanesBuilder.build_plane_parameters(
            p1=[0, 0, z_max], p2=[1, 1, z_max], p3=[1, 0, z_max])
        coordinates = cls.filter_coordinates_related_to_plane(
            coordinates, *min_plane_params, direction=1
        )

        return coordinates
