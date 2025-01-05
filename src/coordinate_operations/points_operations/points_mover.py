import numpy as np

from src.base_structure_classes import Points


class PointsMover:
    @staticmethod
    def move_on_vector(
            points: Points,
            vector: np.ndarray,
            axis: list[str] = ["x", "y", "z"],
    ) -> Points:
        """
        Moves points on vector along provided axis.
        Returnes moved points.
        """

        translated_inner_points: Points = points.copy()

        if "x" in axis:
            translated_inner_points.points[:, 0] += vector[0]
        if "y" in axis:
            translated_inner_points.points[:, 1] += vector[1]
        if "z" in axis:
            translated_inner_points.points[:, 1] += vector[2]

        return translated_inner_points
