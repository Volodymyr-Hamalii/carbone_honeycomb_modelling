from dataclasses import dataclass
from functools import cached_property
import numpy as np

from .coordinate_limits import CoordinateLimits


@dataclass(frozen=True)
class Points:
    """ Template for any class with points array as a property. """
    points: np.ndarray

    def __len__(self) -> int:
        return len(self.points)

    @cached_property
    def coordinate_limits(self) -> CoordinateLimits:
        """ Returns CoordinateLimits of self.points. """

        if len(self.points) == 0:
            return CoordinateLimits()

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
