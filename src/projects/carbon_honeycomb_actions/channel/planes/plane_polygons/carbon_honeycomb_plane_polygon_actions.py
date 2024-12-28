import numpy as np

from src.coordinate_operations import PointsOrganizer, DistanceMeasure

from ....carbon_honeycomb_utils import CarbonHoneycombUtils


class CarbonHoneycombPolygonActions:

    @staticmethod
    def define_polygon_center_coordinates(points: np.ndarray) -> np.ndarray:
        """
        Given a set of points (N, 3) that form a flat polygon in 3D space,
        returns the coordinates of the center (centroid).
        """
        # Ensure points is a 2D array of shape (N, 3)
        if len(points.shape) != 2 or points.shape[1] != 3:
            raise ValueError("points must be of shape (N, 3).")

        # Compute the centroid as the mean of the coordinates
        return points.mean(axis=0)
