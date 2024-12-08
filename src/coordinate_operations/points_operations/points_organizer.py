import numpy as np
from itertools import combinations


class PointsOrganizer:
    @staticmethod
    def _group_by_unique_xy(
        coordinates: np.ndarray
    ) -> dict[tuple[np.float32, np.float32], np.ndarray]:
        """
        Returns dict like
        {(x, y): [points_with_these_x_and_y]}
        """

        # Extract x and y coordinates
        xy_coordinates: np.ndarray = coordinates[:, :2]  # Get only x and y columns

        # Find unique x, y pairs and their indices
        unique_xy, indices = np.unique(xy_coordinates, axis=0, return_inverse=True)

        # Group z-coordinates by unique x, y pairs
        groups = {tuple(xy): [] for xy in unique_xy}

        for i, xy_index in enumerate(indices):
            groups[tuple(unique_xy[xy_index])].append(coordinates[i])

        # Convert each group back to NumPy arrays
        return {key: np.array(value) for key, value in groups.items()}

    @staticmethod
    def _group_by_the_xy_lines(
        groups_by_xy: dict[tuple[np.float32, np.float32], np.ndarray],
        epsilon: float = 1e-4,
        min_points_in_line: int = 2,
    ) -> list[dict[tuple[np.float32, np.float32], np.ndarray]]:
        """
        Groups the coordinates into sets of points that lie on lines in the xOy plane.

        Takes 2 params:
        groups_by_xy - result of the PointsOrganizer._group_by_unique_xy(coordinates)
        epsilon - precision of the rounding to group points.
        """
        points = list(groups_by_xy.keys())
        grouped_lines: list[dict[tuple[np.float32, np.float32], np.ndarray]] = []

        # Check all combinations of 2 points to define candidate lines
        for p1, p2 in combinations(points, 2):
            line_group: dict[tuple[np.float32, np.float32], np.ndarray] = {
                p1: groups_by_xy[p1],
                p2: groups_by_xy[p2]
            }

            for p3 in points:
                if p3 == p1 or p3 == p2:
                    continue

                # Check if p3 is on the line formed by p1 and p2
                det: np.floating = (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p1[1] - p3[1]) * (p2[0] - p3[0])
                if abs(det) < epsilon:
                    line_group[p3] = groups_by_xy[p3]

            # Add to the grouped_lines if it contains at least 3 points
            if len(line_group) >= min_points_in_line:
                # Check if it's a new group (no complete overlap with existing groups)
                if not any(set(line_group.keys()).issubset(set(existing.keys())) for existing in grouped_lines):
                    grouped_lines.append(line_group)

        return grouped_lines
