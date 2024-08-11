import numpy as np
from numpy import ndarray
from matplotlib.pyplot import Axes  # type: ignore
from mpl_toolkits.mplot3d.art3d import Line3DCollection

from scipy.spatial.distance import pdist, squareform


class LinesBuilder:
    @classmethod
    def add_lines(
        cls,
        coordinates: ndarray,
        ax: Axes,
        color_bonds: str = "black",
        num_of_min_distances: int = 3,
    ):
        lines: list[list[ndarray]] = LinesBuilder.build_lines(
            coordinates, num_of_min_distances=num_of_min_distances)

        lc = Line3DCollection(lines, colors=color_bonds, linewidths=1)
        ax.add_collection3d(lc)  # type: ignore

    @classmethod
    def build_lines(
            cls,
            coordinates: ndarray,
            num_of_min_distances: int = 3,
    ) -> list[list[ndarray]]:
        """
        Build lines between points (like, bonds between atoms).

        Return lines: list[list[ndarray]]

        To add them to the structure add lines like
            lc = Line3DCollection(lines, colors='black', linewidths=1)
            ax.add_collection3d(lc)  # ax: Axes
        """
        lines: list[list[ndarray]] = []

        # Calculate the distance matrix for all atoms
        distances_matrix: ndarray = squareform(pdist(coordinates))

        # Round all values to 5th number of decimal place (to avoid duplicates like 1.44000006 and 1.44000053)
        distances_matrix = np.round(distances_matrix, decimals=2)

        # Add a large value to the diagonal to ignore self-distances
        np.fill_diagonal(distances_matrix, np.inf)

        min_distances: ndarray = cls._find_min_unique_values(arr=distances_matrix, num_of_values=num_of_min_distances)

        # Find the nearest neighbors for all atoms
        # nearest_neighbors = cls._find_nearest_neighbors(distances_matrix, num_neighbors)

        # Prepare the lines for the nearest neighbors
        # for i, neighbors in enumerate(nearest_neighbors):
        #    for neighbor in neighbors:
        #        # Append the lines between atoms
        #        lines.append([coordinates[i], coordinates[neighbor]])

        # Iterate over the distance matrix to find atom pairs with distances in min_distances
        for i in range(len(coordinates)):
            for j in range(i + 1, len(coordinates)):
                if distances_matrix[i, j] in min_distances:
                    # Append the line between atoms i and j
                    lines.append([coordinates[i], coordinates[j]])

        return lines

    @staticmethod
    def _find_nearest_neighbors(distances_matrix: ndarray, num_neighbors: int) -> ndarray:
        """ Find the indices of the nearest neighbors """
        return np.argsort(distances_matrix, axis=1)[:, :num_neighbors]

    @staticmethod
    def _find_min_unique_values(arr: ndarray, num_of_values: int) -> ndarray:
        # Flatten the array to 1D and extract unique values
        unique_values: ndarray = np.unique(arr)

        # Sort the unique values
        sorted_unique_values: ndarray = np.sort(unique_values)

        # Return the first num_of_values values
        return sorted_unique_values[:num_of_values]
