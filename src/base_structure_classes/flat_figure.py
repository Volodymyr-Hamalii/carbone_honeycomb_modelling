import numpy as np
from dataclasses import dataclass

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

