import numpy as np

from src.base_structure_classes import Points


class PointsRotator:
    @classmethod
    def rotate_on_angle_related_center(
            cls,
            points: Points,
            angle_x: float = 0,
            angle_y: float = 0,
            angle_z: float = 0,
    ) -> Points:
        if angle_x == 0 and angle_y == 0 and angle_z == 0:
            return points

        # Calculate the centroid (center) of the points
        centroid: np.ndarray = np.mean(points.points, axis=0)

        # Move the points to the origin (centroid becomes the origin)
        centered_points: np.ndarray = points.points - centroid

        # Define rotation matrices for X, Y, and Z axis
        rotation_matrix_x: np.ndarray = cls._get_rotation_matrix_x(angle_x)
        rotation_matrix_y: np.ndarray = cls._get_rotation_matrix_y(angle_y)
        rotation_matrix_z: np.ndarray = cls._get_rotation_matrix_z(angle_z)

        # Combine the rotation matrices
        rotation_matrix: np.ndarray = rotation_matrix_z @ rotation_matrix_y @ rotation_matrix_x

        # Apply the rotation to the points
        rotated_points: np.ndarray = centered_points @ rotation_matrix.T

        # Move the points back to the original position (reverse translation)
        final_points: np.ndarray = rotated_points + centroid

        return Points(
            points=final_points
        )

    @staticmethod
    def _get_rotation_matrix_x(angle_x: float) -> np.ndarray:
        return np.array([
            [1, 0, 0],
            [0, np.cos(angle_x), -np.sin(angle_x)],
            [0, np.sin(angle_x), np.cos(angle_x)]
        ])

    @staticmethod
    def _get_rotation_matrix_y(angle_y: float) -> np.ndarray:
        return np.array([
            [np.cos(angle_y), 0, np.sin(angle_y)],
            [0, 1, 0],
            [-np.sin(angle_y), 0, np.cos(angle_y)]
        ])

    @staticmethod
    def _get_rotation_matrix_z(angle_z: float) -> np.ndarray:
        return np.array([
            [np.cos(angle_z), -np.sin(angle_z), 0],
            [np.sin(angle_z), np.cos(angle_z), 0],
            [0, 0, 1]
        ])
