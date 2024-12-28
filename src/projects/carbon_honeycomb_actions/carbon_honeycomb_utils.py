import numpy as np
# from collections import defaultdict
# from scipy.spatial.distance import cdist


from src.coordinate_operations import PointsOrganizer, DistanceMeasure
from src.structure_visualizer import StructureVisualizer


class CarbonHoneycombUtils:

    @staticmethod
    def find_end_points(
        points_on_line: np.ndarray
    ) -> tuple[tuple[np.float32, np.float32], tuple[np.float32, np.float32]]:
        first_point: np.ndarray = points_on_line[0]
        
        start: tuple[np.float32, np.float32] = first_point[0], first_point[1]
        end: tuple[np.float32, np.float32] = start

        for xy_points in points_on_line[1:]:
                x, y = xy_points
                sx, sy = start
                ex, ey = end

                if x == sx and x == ex:
                    if y < sy:
                        start = (x, y)
                    elif y > ey:
                        end = (x, y)
                elif x < sx:
                    start = (x, y)
                elif x > ex:
                    end = (x, y)

        return start, end

    @classmethod
    def find_end_points_of_honeycomb_planes_groups(
        cls, honeycomb_planes_groups: list[dict[tuple[np.float32, np.float32], np.ndarray]]
    ) -> list[tuple[tuple[np.float32, np.float32], tuple[np.float32, np.float32]]]:
        """
        Returns the coordinates of the ends of the segments
        on the xOy plane that form the channel planes.

        Keep the same index order as in honeycomb_planes_groups.
        """

        end_points_of_groups: list[tuple[tuple[np.float32, np.float32], tuple[np.float32, np.float32]]] = []

        for group in honeycomb_planes_groups:
            keys: np.ndarray = np.array(list(group.keys()))
            start, end = cls.find_end_points(points_on_line=keys)
            end_points_of_groups.append((start, end))

        # return sorted(end_points_of_groups, key=lambda p: p[0][0])
        return end_points_of_groups

    @staticmethod
    def found_hexagone_node_indexes(
        extreme_points_of_groups: list[tuple[tuple[np.float32, np.float32], tuple[np.float32, np.float32]]]
    ) -> list[list[int]]:
        """ Fount the lines that form hexagones and returns the list of indexes of such groups. """

        # Step 1: Extract all unique nodes (points) and map them to integer indices
        # Collect all points
        points_set = set()
        for i, (p1, p2) in enumerate(extreme_points_of_groups):
            points_set.add(p1)
            points_set.add(p2)
        points_list = list(points_set)
        point_to_index = {p: idx for idx, p in enumerate(points_list)}

        # Step 2: Build adjacency list from edges
        n = len(points_list)
        adjacency = [[] for _ in range(n)]
        for i, (p1, p2) in enumerate(extreme_points_of_groups):
            u = point_to_index[p1]
            v = point_to_index[p2]
            adjacency[u].append((v, i))  # store (node, edge_index)
            adjacency[v].append((u, i))  # undirected

        # Step 3: Find all 6-cycles
        # A brute force backtracking approach:
        hexagons = set()  # use a set to avoid duplicates, will store tuple of sorted edges

        def backtrack(start_node, current_node, visited_nodes, visited_edges):
            # If we have a path of length 6 (6 edges means 7 nodes), check if it forms a cycle back to start_node
            if len(visited_nodes) == 7:
                # Check if last node connects to start_node
                if visited_nodes[-1] == start_node:
                    # We found a 6-cycle
                    # visited_edges has the edge indices in the order they were visited
                    cycle_edges = tuple(sorted(visited_edges))
                    hexagons.add(cycle_edges)
                return

            if len(visited_nodes) > 7:
                return  # longer than needed

            # Explore neighbors
            for (nxt, eidx) in adjacency[current_node]:
                # We must not revisit nodes to ensure a simple cycle, except possibly to close the cycle
                if nxt == start_node and len(visited_nodes) == 6:
                    # Can we close the cycle with start_node?
                    # This adds one more edge, completing a 6-cycle
                    new_visited_edges = visited_edges + [eidx]
                    cycle_edges = tuple(sorted(new_visited_edges))
                    hexagons.add(cycle_edges)
                elif nxt not in visited_nodes:
                    # Continue path
                    # Add node and edge
                    backtrack(start_node, nxt, visited_nodes + [nxt], visited_edges + [eidx])

        # Try starting from each node and find 6-cycles
        # To avoid recounting the same cycle multiple times starting from different nodes,
        # we can impose an order: start from a specific node and only proceed to neighbors with greater index, etc.
        # However, using the set of sorted edges as a key will eliminate duplicates anyway.
        for start in range(n):
            backtrack(start, start, [start], [])

        # Step 4: Convert the set of cycles back into lists of edges
        # Each element in hexagons is a tuple of edge indices sorted.
        # If needed, we can return them as lists in sorted order.
        result: list[list[int]] = [list(cycle) for cycle in hexagons]

       # Remove duplicates by converting each list to a frozenset
        unique_cycles = {frozenset(cycle) for cycle in result}

        # Convert each frozenset back to a list
        result = [list(cycle) for cycle in unique_cycles]

        return result

    @staticmethod
    def split_xy_groups_by_max_distances(
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
    def split_xy_groups_by_max_distances(
        points_grouped_by_lines: list[np.ndarray],
        max_distance_between_xy_groups: np.floating | float,
    ) -> list[np.ndarray]:
        ...
