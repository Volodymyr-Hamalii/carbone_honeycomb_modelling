import numpy as np


class PlanesBuilder:
    @staticmethod
    def build_plane_params(
            p1: np.ndarray | list[float],
            p2: np.ndarray | list[float],
            p3: np.ndarray | list[float],
    ) -> tuple[float, float, float, float]:
        """
        Define the plane like Ax + By + Cz + D = 0
        using the three provided points.

        Takes 3 points as a parameters as lists with 3 coordinates.
        Returns A, B, C, D parameters from the equation above.
        """

        p1_np: np.ndarray = np.array(p1)
        p2_np: np.ndarray = np.array(p2)
        p3_np: np.ndarray = np.array(p3)

        # Calculate the normal vector of the plane
        v1: np.ndarray = p2_np - p1_np
        v2: np.ndarray = p3_np - p1_np
        normal: np.ndarray = np.cross(v1, v2)

        A: float = normal[0]
        B: float = normal[1]
        C: float = normal[2]
        D: float = -np.dot(normal, p1_np)

        if B < 0:
            # Actually, I don't know why we need this check, but atoms filtering related planes
            # only works correctly under this condition.
            return -A, -B, -C, -D

        return A, B, C, D
