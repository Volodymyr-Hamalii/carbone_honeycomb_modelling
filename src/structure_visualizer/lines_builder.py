import numpy as np
from numpy import ndarray, float16

from scipy.spatial.distance import pdist, squareform

class LinesBuilder:
    @classmethod
    def build_lines(cls, coordinates: np.ndarray, num_neighbors=2) -> list[list[ndarray]]:
        """
        Build lines between points (like, bonds between atoms).
         
        """
        lines: list[list[ndarray]] = []

        # Calculate the distance matrix for all atoms
        distances_matrix: np.ndarray = squareform(pdist(coordinates))

        # Find the nearest neighbors for all atoms
        nearest_neighbors = cls._find_nearest_neighbors(distances_matrix, num_neighbors)

        # Prepare the lines for the nearest neighbors
        for i, neighbors in enumerate(nearest_neighbors):
            for neighbor in neighbors:
                # Append the lines between atoms
                lines.append([coordinates[i], coordinates[neighbor]])

        return lines

    @staticmethod
    def _find_nearest_neighbors(distances_matrix: ndarray, num_neighbors: int) -> np.ndarray:
        # Add a large value to the diagonal to ignore self-distances
        np.fill_diagonal(distances_matrix, np.inf)

        # Find the indices of the nearest neighbors
        nearest_neighbors = np.argsort(distances_matrix, axis=1)[:, :num_neighbors]

        return nearest_neighbors
