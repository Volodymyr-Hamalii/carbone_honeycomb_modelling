import numpy as np
# from collections import defaultdict
# from scipy.spatial.distance import cdist


from src.coordinate_operations import PointsOrganizer, DistanceMeasure
from src.structure_visualizer import StructureVisualizer
from .find_the_shortest_cycle import find_all_shortest_cycles

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

        # Filter not allowed lengths and duplicates
        filtered_honeycomb_planes_groups: list[dict[tuple[np.float32, np.float32], np.ndarray]] = []
        for i in honeycomb_planes_groups:
            if i not in filtered_honeycomb_planes_groups and len(i.values()) in allowed_lengths:
                filtered_honeycomb_planes_groups.append(i)

        return filtered_honeycomb_planes_groups

    # @staticmethod
    # def _build_common_points_matrix(
    #         honeycomb_planes_groups: list[dict[tuple[np.float32, np.float32], np.ndarray]],
    # ) -> np.ndarray:
    #     """
    #     Build a matrix of size N x N (where N is the number of groups),
    #     indicating whether two groups share at least one (x,y) coordinate.

    #     Returns:
    #         np.ndarray: A boolean matrix where:
    #                     matrix[i, j] = True if honeycomb_planes_groups[i] and honeycomb_planes_groups[j]
    #                     share at least one (x,y) point. Otherwise, False.
    #     """
    #     n = len(honeycomb_planes_groups)
    #     # Convert each group's keys to a set for efficient intersection
    #     sets_of_points = [set(group.keys()) for group in honeycomb_planes_groups]

    #     # Initialize the matrix with False
    #     common_matrix = np.zeros((n, n), dtype=bool)

    #     # Check each pair of groups
    #     for i in range(n):
    #         for j in range(n):
    #             if i != j:
    #                 # Check if there's any intersection in keys
    #                 if sets_of_points[i].intersection(sets_of_points[j]):
    #                     common_matrix[i, j] = True

    #     return common_matrix

    # @classmethod
    # def _get_channel_groups(
    #     cls,
    #     honeycomb_planes_groups: list[dict[tuple[np.float32, np.float32], np.ndarray]]
    # ) -> list[list[int]]:
    #     """Get the list of index groups that form channels (connected components that contain at least one cycle)."""
    #     common_points_matrix: np.ndarray = cls._build_common_points_matrix(honeycomb_planes_groups)
    #     return find_all_shortest_cycles(common_points_matrix)
        # n = len(honeycomb_planes_groups)

        # visited = [False] * n
        # channel_groups = []

        # def dfs(start):
        #     stack = [start]
        #     component_nodes = []
        #     edges_count = 0  # We'll count edges by summing adjacency in the component
        #     visited[start] = True

        #     while stack:
        #         node = stack.pop()
        #         component_nodes.append(node)
        #         # Count the edges: for each 'True' in that row, we have an edge
        #         # We'll count only edges from this node to others (avoid double-counting)
        #         edges_count += np.count_nonzero(common_points_matrix[node])
        #         for nxt, connected in enumerate(common_points_matrix[node]):
        #             if connected and not visited[nxt]:
        #                 visited[nxt] = True
        #                 stack.append(nxt)

        #     # Each edge counted twice (once from each endpoint), so divide by 2
        #     edges_count //= 2
        #     return component_nodes, edges_count

        # for i in range(n):
        #     if not visited[i]:
        #         component_nodes, edges_count = dfs(i)
        #         # Check if this component forms a cycle
        #         # A tree with V nodes has V-1 edges. If edges_count >= len(component_nodes),
        #         # then we have at least one cycle.
        #         if edges_count >= len(component_nodes):
        #             channel_groups.append(component_nodes)

        # return channel_groups

    @staticmethod
    def _find_end_points_of_groups(
        honeycomb_planes_groups: list[dict[tuple[np.float32, np.float32], np.ndarray]]
    ) -> list[tuple[tuple[np.float32, np.float32], tuple[np.float32, np.float32]]]:
        """
        Returns the coordinates of the ends of the segments
        on the xOy plane that form the channel planes.

        Keep the same index order as in honeycomb_planes_groups.
        """

        end_points_of_groups: list[tuple[tuple[np.float32, np.float32], tuple[np.float32, np.float32]]] = []

        for group in honeycomb_planes_groups:
            keys = list(group.keys())

            start: tuple[np.float32, np.float32] = keys[0]
            end: tuple[np.float32, np.float32] = keys[0]

            for xy_points in keys[1:]:
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

            end_points_of_groups.append((start, end))

        # return sorted(end_points_of_groups, key=lambda p: p[0][0])
        return end_points_of_groups

    @staticmethod
    def _found_hexagone_node_indexes(
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
        # channel_groups: list[list[int]] = cls._get_channel_groups(honeycomb_planes_groups)

        end_points_of_groups: list[
            tuple[tuple[np.float32, np.float32], tuple[np.float32, np.float32]]
        ] = cls._find_end_points_of_groups(honeycomb_planes_groups)

        plane_group_indexes: list[list[int]] = cls._found_hexagone_node_indexes(end_points_of_groups)

        StructureVisualizer.show_2d_graph(x_y_points, show_coordinates=True)
