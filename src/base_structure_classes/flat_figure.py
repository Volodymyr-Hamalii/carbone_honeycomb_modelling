import numpy as np
from dataclasses import dataclass

from src.coordinate_operations import PlanesBuilder
from .points import Points


@dataclass
class FlatFigure(Points):
    @property
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

    def get_plane_params(self) -> tuple[float, float, float, float]:
        """
        Define the plane like Ax + By + Cz + D = 0
        using the three provided points.

        Takes 3 points as a parameters as lists with 3 coordinates.
        Returns A, B, C, D parameters from the equation above.
        """

        plane_points: list[np.ndarray] = []
        used_x = set()
        used_y = set()
        used_z = set()

        for point in self.points:
            x, y, z = point

            if (x not in used_x) and (y not in used_y) and (z not in used_z):
                plane_points.append(point)
                used_x.add(x)
                used_y.add(y)
                used_z.add(z)

                if len(plane_points) == 3:
                    break

        if len(plane_points) < 3:
            raise ValueError("Not enough distinct points to define a plane.")

        return PlanesBuilder.build_plane_params(
            plane_points[0], plane_points[1], plane_points[2])
