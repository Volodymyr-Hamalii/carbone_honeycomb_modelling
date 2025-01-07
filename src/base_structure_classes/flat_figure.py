from functools import cached_property
from dataclasses import dataclass
import numpy as np

from src.base_structure_classes import CoordinateLimits
from src.coordinate_operations import PlanesBuilder, LinesOperations
from .points import Points


@dataclass
class FlatFigure(Points):
    @cached_property
    def center(self) -> np.ndarray:
        """
        Given a set of points (N, 3) that form a flat figure in 3D space,
        returns the coordinates of the center (centroid).
        """
        # Ensure points is a 2D array of shape (N, 3)
        if len(self.points.shape) != 2 or self.points.shape[1] != 3:
            raise ValueError("self.points must be of shape (N, 3).")

        # Compute the centroid as the mean of the coordinates
        return self.points.mean(axis=0)

    @cached_property
    def plane_params(self) -> tuple[float, float, float, float]:
        """ Compute and cache result of get_plane_params() method. """
        return self.get_plane_params()

    def get_plane_params(self) -> tuple[float, float, float, float]:
        """
        Define the plane like Ax + By + Cz + D = 0
        using the three provided points.

        Takes 3 points as a parameters as lists with 3 coordinates.
        Returns A, B, C, D parameters from the equation above.
        """

        points: np.ndarray = self.points
        coordinate_limits: CoordinateLimits = self.coordinate_limits

        if points.shape[0] < 3:
            raise ValueError("At least 3 points are required to define a plane.")

        # 1) Pick first point p1
        p1 = points[0]

        # 2) Find a second point p2 that is not identical to p1
        p2 = None
        for i in range(1, len(points)):
            if not np.allclose(points[i], p1):
                p2 = points[i]
                break
        if p2 is None:
            raise ValueError("All points are identical. Cannot define a plane.")

        p3 = None
        for j in reversed(range(1, len(points))):
            c = points[j]
            # We already used p1, p2, but we might revisit them in the loop;
            # that's okay if they help us find a non-collinear point.
            if not LinesOperations.points_are_collinear(p1, p2, c):
                p3 = c
                break

        if p3 is None:
            raise ValueError("All points are collinear. Cannot define a plane.")

        return PlanesBuilder.build_plane_params(p1, p2, p3)
