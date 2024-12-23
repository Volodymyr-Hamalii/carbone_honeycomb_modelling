import numpy as np
from collections import defaultdict
from dataclasses import dataclass

from src.coordinate_operations import PointsOrganizer
from .planes import CarbonHoneycombPlane


@dataclass
class CarbonHoneycombChannel:
    points: np.ndarray

    @property
    def planes(self) -> list[CarbonHoneycombPlane]:

        # TODO: refactor this "property"

        groups_by_the_xy_lines: list[dict[tuple[np.float32, np.float32], np.ndarray]] = PointsOrganizer.group_by_the_xy_lines(
            coordinates_to_group=self.points, epsilon=1e-1, min_points_in_line=3)

        # Haneycomb plane groups on the on the xOy coordinate plane
        xy_sets: list[set] = [
            set(point_group.keys()) for point_group in groups_by_the_xy_lines
        ]

        point_group_neighbours: dict[int, list] = defaultdict(list)

        prev_plane_index: int = -1
        next_plane_index: int = -1

        for i, i_point_group in enumerate(xy_sets):
            for j, j_point_group in enumerate(xy_sets):
                if i != j and i_point_group & j_point_group:
                    point_group_neighbours[i].append(j)

            # Check first and second plane indexes
            if (0., 0.) in i_point_group:
                # Check if it's base plane
                if prev_plane_index == -1:
                    for point in i_point_group:
                        # Check at least one more point
                        if point[0] != 0. and point[1] == 0.:
                            prev_plane_index = i
                            break

                    # Check if it's second plane
                    if prev_plane_index == -1 and next_plane_index == -1:
                        next_plane_index = i
                else:
                    next_plane_index = i

        if prev_plane_index == -1:
            raise ValueError("Error during searching for the first plane index: index not found.")

        if next_plane_index == -1:
            raise ValueError("Error during searching for the second plane index: index not found.")

        planes: list[CarbonHoneycombPlane] = []
        for i in range(len(point_group_neighbours)):
            direction: bool = i % 2 == 0

            if len(planes) == 0:
                index: int = prev_plane_index
            else:
                index = next_plane_index

                for i in point_group_neighbours[next_plane_index]:
                    if i != prev_plane_index:
                        next_plane_index = i
                        break

                prev_plane_index = index

            plane_points: np.ndarray = np.array([
                point for point in groups_by_the_xy_lines[index].values()
            ])

            planes.append(
                CarbonHoneycombPlane(
                    points=np.concatenate(plane_points, axis=0),
                    direction=direction
                )
            )

        return planes

