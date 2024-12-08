import numpy as np
# from collections import defaultdict
# from scipy.spatial.distance import cdist


from src.coordinate_operations import PointsOrganizer, DistanceMeasure
from src.structure_visualizer import StructureVisualizer

from .channel import CarbonHoneycombChannel


class CarbonHoneycombActions:

    @staticmethod
    def _split_groups_by_max_distances(
        groups_by_the_xy_lines: list[dict[tuple[np.float32, np.float32], np.ndarray]],
        max_distance_between_xy_groups: np.floating | float,
    ) -> list[dict[tuple[np.float32, np.float32], np.ndarray]]:
        """
        Takes groups_by_the_xy_lines and split the group into separate
        if there are some groups with the min distances more than max_distance_between_xy_groups.
        """
        result: list[dict[tuple[np.float32, np.float32], np.ndarray]] = []

        for line_dict in groups_by_the_xy_lines:
            # Sort keys (usually by x then y)
            keys: list[tuple[np.floating, np.floating]] = sorted(
                line_dict.keys(), key=lambda p: (float(p[0]), float(p[1])))

            clusters = []
            current_cluster: list[tuple[np.floating, np.floating]] = [keys[0]]

            # Build clusters based on max distance
            for i in range(1, len(keys)):
                dist: float = DistanceMeasure.calculate_distance_between_2_points(keys[i - 1], keys[i])
                if dist <= max_distance_between_xy_groups:
                    # Same cluster
                    current_cluster.append(keys[i])
                else:
                    # New cluster starts here
                    clusters.append(current_cluster)
                    current_cluster = [keys[i]]

            # Append the last cluster if not empty
            if current_cluster:
                clusters.append(current_cluster)

            # Filter out clusters that have only one point
            filtered_clusters = [c for c in clusters if len(c) > 1]

            # Convert each cluster back to a dict
            for cluster in filtered_clusters:
                cluster_dict = {k: line_dict[k] for k in cluster}
                result.append(cluster_dict)

        return result

    @staticmethod
    def _filter_honeycomb_planes_groups(
        honeycomb_planes_groups: list[dict[tuple[np.float32, np.float32], np.ndarray]]
    ) -> list[dict[tuple[np.float32, np.float32], np.ndarray]]:

        lengths: list[int] = [len(i.values()) for i in honeycomb_planes_groups]

        # TODO: refactor (use list instead of the dict)
        lengths_counter_map: dict[int, int] = {}
        for i in lengths:
            if i not in lengths_counter_map:
                lengths_counter_map[i] = 0
            lengths_counter_map[i] += 1

        allowed_lengths: list[int] = []
        for i in sorted(lengths_counter_map.keys(), reverse=True):
            if (len(allowed_lengths) < 3) or (i > 4):
                allowed_lengths.append(i)

        return [
            i for i in honeycomb_planes_groups
            if len(i.values()) in allowed_lengths
        ]

    @classmethod
    def split_init_structure_into_separate_channels(
            cls, carbon_coordinates: np.ndarray
    ) -> list[CarbonHoneycombChannel]:

        # Create groups with the same X and Y coordinate
        # (like, to split all points into columns)
        groups_by_xy: dict[
            tuple[np.float32, np.float32], np.ndarray
        ] = PointsOrganizer.group_by_unique_xy(carbon_coordinates)

        # Define groups that lie on the same line
        groups_by_the_xy_lines: list[
            dict[tuple[np.float32, np.float32], np.ndarray]
        ] = PointsOrganizer.group_by_the_xy_lines(groups_by_xy, epsilon=1e-2, min_points_in_line=3)

        # Define 2D honeycomb structure
        x_y_points: np.ndarray = np.array([
            [i[0], i[1]] for i in groups_by_xy.keys()
        ])

        # Split by the max distance between groups (to define separate channel planes)
        distances_between_xy_groups: np.ndarray = DistanceMeasure.calculate_min_distances_between_points(x_y_points)
        max_distance_between_xy_groups: np.floating = np.max(distances_between_xy_groups) * 1.5

        honeycomb_planes_groups: list[
            dict[tuple[np.float32, np.float32], np.ndarray]
        ] = cls._split_groups_by_max_distances(groups_by_the_xy_lines, max_distance_between_xy_groups)

        honeycomb_planes_groups = cls._filter_honeycomb_planes_groups(honeycomb_planes_groups)

        StructureVisualizer.show_2d_graph(x_y_points, show_coordinates=True)
