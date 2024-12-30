import numpy as np

from src.utils import Logger
# from src.base_structure_classes import Points


logger = Logger(__name__)


class LinesOperations:
    @staticmethod
    def points_are_collinear(a: np.ndarray, b: np.ndarray, c: np.ndarray, eps: float = 1e-6) -> bool:
        """
        Returns True if points a, b, c lie on the same line in 3D,
        i.e. the norm of the cross product is near zero.
        """
        cross_prod: np.ndarray = np.cross(b - a, c - a)
        return bool(np.linalg.norm(cross_prod) < eps)
