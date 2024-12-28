import numpy as np

from src.coordinate_operations import PointsOrganizer, DistanceMeasure

from ...carbon_honeycomb_utils import CarbonHoneycombUtils
from .carbon_honeycomb_plane_hexagon import CarbonHoneycombHexagon


class CarbonHoneycombPlaneActions:
    @classmethod
    def define_plane_hexagons(cls, points: np.ndarray) -> list[CarbonHoneycombHexagon]:
        points_grouped_by_lines: list[np.ndarray] = PointsOrganizer.group_by_lines(points)

        distances_between_points: np.ndarray = DistanceMeasure.calculate_min_distances_between_points(points)
        max_distance_between_points: np.floating = np.max(distances_between_points) * 1.25

        points_grouped_by_lines = CarbonHoneycombUtils.split_groups_by_max_distances(
            points_grouped_by_lines, max_distance_between_points)

        end_points: list[tuple[tuple, tuple]] = cls._find_end_points(points_grouped_by_lines)

        plane_hexagons_indexes: list[list[int]] = CarbonHoneycombUtils.found_hexagone_node_indexes(end_points)

        plane_hexagons: list[CarbonHoneycombHexagon] = cls._build_plane_hexagons(
            points_grouped_by_lines, plane_hexagons_indexes)
        return plane_hexagons

    @staticmethod
    def _find_end_points(points_grouped_by_lines: list[np.ndarray]) -> list[tuple[tuple, tuple]]:

        end_points_of_groups: list[tuple[tuple, tuple]] = []

        for group in points_grouped_by_lines:
            coordinates_to_check: tuple[str, str] = ("x", "z")  # default

            # Chech if all x coordinates are equal
            x_coord: np.ndarray = group[:, 0]
            all_x_equal: np.bool_ = np.all(x_coord == x_coord[0])

            if all_x_equal:
                # Replace "x" with "y":
                coordinates_to_check: tuple[str, str] = ("y", "z")

            start, end = CarbonHoneycombUtils.find_end_points(
                points_on_line=group, coordinates_to_check=coordinates_to_check)
            end_points_of_groups.append((start, end))

        # return sorted(end_points_of_groups, key=lambda p: p[0][0])
        return end_points_of_groups

    @staticmethod
    def _build_plane_hexagons(
        points_grouped_by_lines: list[np.ndarray],
        plane_hexagons_indexes: list[list[int]],
    ) -> list[CarbonHoneycombHexagon]:
        plane_hexagons: list[CarbonHoneycombHexagon] = []

        for plane_group_indexes in plane_hexagons_indexes:
            unique_points_set = set()

            for i in plane_group_indexes:
                for point in points_grouped_by_lines[i]:
                    # Convert to tuple to store in a set
                    point_tuple = tuple(point)
                    unique_points_set.add(point_tuple)

            # Convert the set of tuples back to a numpy array
            honeycomb_points: np.ndarray = np.array(list(unique_points_set))
            honeycomb_channel = CarbonHoneycombHexagon(points=honeycomb_points)
            plane_hexagons.append(honeycomb_channel)

        return plane_hexagons
