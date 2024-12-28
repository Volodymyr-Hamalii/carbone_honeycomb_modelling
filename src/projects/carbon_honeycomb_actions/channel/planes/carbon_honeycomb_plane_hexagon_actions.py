import numpy as np

from src.coordinate_operations import PointsOrganizer, DistanceMeasure

from ...carbon_honeycomb_utils import CarbonHoneycombUtils

class CarbonHoneycombHexagonActions:
    
    @staticmethod
    def define_center_coordinates(hexagon_points: np.ndarray) -> np.ndarray:
        """
        Given a set of points (N, 3) that form a hexagon in 3D space (all points lie in the same plane),
        returns the coordinates of the center (centroid).
        """
        # Ensure hexagon_points is a 2D array of shape (N, 3)
        if len(hexagon_points.shape) != 2 or hexagon_points.shape[1] != 3:
            raise ValueError("hexagon_points must be of shape (N, 3).")

        # Compute the centroid as the mean of the coordinates
        return hexagon_points.mean(axis=0)
