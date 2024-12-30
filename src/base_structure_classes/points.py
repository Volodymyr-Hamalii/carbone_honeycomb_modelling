from dataclasses import dataclass
import numpy as np

from .coordinate_limits import CoordinateLimits


@dataclass
class Points:
    """ Template for any class with points array as a property. """
    points: np.ndarray

    def __len__(self) -> int:
        return len(self.points)

    def get_coordinate_limits(self) -> CoordinateLimits:
        """ Returns CoordinateLimits of self.points. """

        x_coords: np.ndarray = self.points[:, 0]
        y_coords: np.ndarray = self.points[:, 1]
        z_coords: np.ndarray = self.points[:, 2]

        return CoordinateLimits(
            x_min=np.min(x_coords),
            x_max=np.max(x_coords),

            y_min=np.min(y_coords),
            y_max=np.max(y_coords),

            z_min=np.min(z_coords),
            z_max=np.max(z_coords),
        )

    def copy(self) -> "Points":
        """ Returns a new Points instance. """
        return Points(points=self.points.copy())
