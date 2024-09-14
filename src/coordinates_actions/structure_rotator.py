import numpy as np
from numpy import ndarray

class StructureRotator:
    @classmethod
    def rotate_on_angle_related_center(
            cls, points: ndarray, angle_x: float = 0, angle_y: float = 0, angle_z: float = 0) -> ndarray:
        # Calculate the centroid (center) of the points
        centroid = np.mean(points, axis=0)

        # Move the points to the origin (centroid becomes the origin)
        centered_points = points - centroid

        # Define rotation matrices for X, Y, and Z axis
        rotation_matrix_x = np.array([
            [1, 0, 0],
            [0, np.cos(angle_x), -np.sin(angle_x)],
            [0, np.sin(angle_x), np.cos(angle_x)]
        ])

        rotation_matrix_y = np.array([
            [np.cos(angle_y), 0, np.sin(angle_y)],
            [0, 1, 0],
            [-np.sin(angle_y), 0, np.cos(angle_y)]
        ])

        rotation_matrix_z = np.array([
            [np.cos(angle_z), -np.sin(angle_z), 0],
            [np.sin(angle_z), np.cos(angle_z), 0],
            [0, 0, 1]
        ])

        # Combine the rotation matrices
        rotation_matrix = rotation_matrix_z @ rotation_matrix_y @ rotation_matrix_x

        # Apply the rotation to the points
        rotated_points = centered_points @ rotation_matrix.T

        # Move the points back to the original position (reverse translation)
        final_points = rotated_points + centroid

        return final_points
