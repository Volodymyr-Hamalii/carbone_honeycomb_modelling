import numpy as np
from collections import defaultdict
from dataclasses import dataclass

from src.utils import Logger
from src.coordinate_operations import PointsOrganizer
from .planes import CarbonHoneycombPlane


logger = Logger(__name__)

@dataclass
class CarbonHoneycombChannel:
    points: np.ndarray

    # def __post_init__(self) -> None:
    #     self.planes: list[CarbonHoneycombPlane] = self._build_planes()

    @property
    def planes(self) -> list[CarbonHoneycombPlane]:
        """
        Returns a list of CarbonHoneycombPlane objects
        representing planes in the honeycomb channel.
        """
        return self._build_planes()

    def _build_planes(self) -> list[CarbonHoneycombPlane]:
        # 1. Group points into lines in the xOy plane
        groups_by_xy_lines: list[dict[tuple[np.float32, np.float32], np.ndarray]] = PointsOrganizer.group_by_the_xy_lines(
            coordinates_to_group=self.points, epsilon=1e-1, min_points_in_line=3)

        # 2. Build sets of (x, y) for each group and neighbor relationships
        xy_sets: list[set[tuple[np.floating, np.floating]]] = [set(group.keys()) for group in groups_by_xy_lines]
        point_group_neighbors: dict[int, list[int]] = self._build_neighbors(xy_sets)

        # 3. Find the first (base) plane index and second plane index
        prev_plane_index, next_plane_index = self._find_first_and_second_plane(xy_sets)

        if prev_plane_index == -1 or next_plane_index == -1:
            logger.warning("Planes not build")
            return []

        # 4. Build the planes by traversing neighbor relationships
        planes: list[CarbonHoneycombPlane] = []
        for i in range(len(point_group_neighbors)):
            # For illustration: alternate direction every other plane
            direction: bool = (i % 2 == 0)

            # If no planes are built yet, pick the "base" plane; otherwise pick the "next"
            if not planes:
                index: int = prev_plane_index
            else:
                index: int = next_plane_index
                # Advance the next_plane_index to a neighbor (avoiding the previous plane)
                for neighbor_index in point_group_neighbors[next_plane_index]:
                    if neighbor_index != prev_plane_index:
                        next_plane_index = neighbor_index
                        break
                prev_plane_index = index

            # 5. Concatenate all arrays in the chosen group to build one plane's coordinates
            plane_points: np.ndarray = np.concatenate(list(groups_by_xy_lines[index].values()), axis=0)

            planes.append(CarbonHoneycombPlane(points=plane_points, direction=direction))

        return planes

    def _build_neighbors(
        self,
        xy_sets: list[set[tuple[np.float32, np.float32]]]
    ) -> dict[int, list[int]]:
        """ Build neighbor relationships between groups that share at least one (x,y). """
        neighbors = defaultdict(list)
        for i, i_points in enumerate(xy_sets):
            for j, j_points in enumerate(xy_sets):
                if i != j and i_points & j_points:
                    neighbors[i].append(j)
        return neighbors

    def _find_first_and_second_plane(
        self,
        xy_sets: list[set[tuple[np.float32, np.float32]]]
    ) -> tuple[int, int]:
        """
        Identify two special planes:
          - The 'base' plane (prev_plane_index): one containing (0,0) plus some other point with y=0.
          - The 'second' plane (next_plane_index): if there's another group with (0,0), pick that too.
        """
        prev_plane_index = -1
        next_plane_index = -1

        for i, group in enumerate(xy_sets):
            # If a group contains (0,0), it might be our base or second plane
            if (0., 0.) in group:
                # If we haven't found the base yet, verify the group has another point with y=0
                if prev_plane_index == -1:
                    for pt in group:
                        # e.g., at least one more point with x != 0 and y == 0
                        if pt[0] != 0. and pt[1] == 0.:
                            prev_plane_index = i
                            break
                    # If we still haven't assigned prev_plane_index, then assign next_plane_index
                    if prev_plane_index == -1 and next_plane_index == -1:
                        next_plane_index = i
                else:
                    # If base is found, the next group containing (0,0) is second plane
                    next_plane_index = i

        # if prev_plane_index == -1:
        #     raise ValueError("Base plane (first plane) not found.")

        # if next_plane_index == -1:
        #     raise ValueError("Second plane not found.")

        return prev_plane_index, next_plane_index
