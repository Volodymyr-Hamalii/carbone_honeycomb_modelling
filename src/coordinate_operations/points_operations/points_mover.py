import numpy as np
from numpy import ndarray


class PointsMover:
    @staticmethod
    def move_on_vector(points: ndarray, vector: ndarray, axis: list[str] = ["x", "y", "z"]) -> ndarray:
        """
        Moves points on vector along provided axis.
        Returnes moved points.
        """

        translated_inner_points: ndarray = points.copy()

        if "x" in axis:
            translated_inner_points[:, 0] += vector[0]
        if "y" in axis:
            translated_inner_points[:, 1] += vector[1]
        if "z" in axis:
            translated_inner_points[:, 1] += vector[2]

        return translated_inner_points
