import numpy as np
from numpy import ndarray


class PlanesBuilder:
    @staticmethod
    def build_plane_params(
            p1: ndarray | list[float],
            p2: ndarray | list[float],
            p3: ndarray | list[float],
    ) -> tuple[float, float, float, float]:
        """
        Define the plane like Ax + By + Cz + D = 0
        using the three provided points.

        Takes 3 points as a parameters as lists with 3 coordinates.
        Returns A, B, C, D parameters from the equation above.
        """

        p1_np: ndarray = np.array(p1)
        p2_np: ndarray = np.array(p2)
        p3_np: ndarray = np.array(p3)

        # Calculate the normal vector of the plane
        v1: ndarray = p2_np - p1_np
        v2: ndarray = p3_np - p1_np
        normal: ndarray = np.cross(v1, v2)
        A, B, C = normal
        D = -np.dot(normal, p1_np)

        return A, B, C, D
