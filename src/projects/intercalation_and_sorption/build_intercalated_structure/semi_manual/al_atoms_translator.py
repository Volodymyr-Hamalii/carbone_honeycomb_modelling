import numpy as np
from scipy.optimize import minimize

from src.projects.carbon_honeycomb_actions.channel.planes.carbon_honeycomb_plane import CarbonHoneycombPlane
from src.utils import Constants, Logger, FileReader
from src.base_structure_classes import Points
from src.coordinate_operations import PointsOrganizer, PointsMover, PointsRotator, DistanceMeasure
from src.projects.carbon_honeycomb_actions import (
    CarbonHoneycombChannel,
    CarbonHoneycombPlane,
)


class AlAtomsTranslator:
    @classmethod
    def translate_for_all_planes(
        cls,
        carbon_channel: CarbonHoneycombChannel,
        al_plane_coordinates: Points,
    ) -> Points:

        planes: list[CarbonHoneycombPlane] = carbon_channel.planes

        al_groups: list[np.ndarray] = cls._group_by_lines(al_plane_coordinates)
        plane_group_map: dict[int, np.ndarray] = cls._match_plane_with_group(
            al_groups, planes)

        all_al_points: Points = cls._copy_al_points_to_rest_planes(plane_group_map, carbon_channel)

        return all_al_points

    @staticmethod
    def _group_by_lines(al_plane_coordinates: Points) -> list[np.ndarray]:
        # Create groups with the same X and Y coordinate
        # (like, to split all points into columns)
        groups_by_xy: dict[
            tuple[np.float32, np.float32], np.ndarray
        ] = PointsOrganizer.group_by_unique_xy(al_plane_coordinates.points)

        # Define groups that lie on the same line
        groups_by_the_xy_lines: list[
            dict[tuple[np.float32, np.float32], np.ndarray]
        ] = PointsOrganizer.group_by_the_xy_lines(groups_by_xy, epsilon=1e-1, min_points_in_line=3)

        grouped_by_lines: list[np.ndarray] = []

        for i in groups_by_the_xy_lines:
            group: np.ndarray = np.concatenate(list(i.values()))
            grouped_by_lines.append(group)

        return grouped_by_lines

    @staticmethod
    def _match_plane_with_group(
        al_groups: list[np.ndarray],
        planes: list[CarbonHoneycombPlane],
    ) -> dict[int, np.ndarray]:
        """ Returns the dict with the plane index and corresponding Al points. """

        plane_group_map: dict[int, np.ndarray] = {}

        for group in al_groups:
            # Find the plane with the minimum distance sum
            closest_plane_index: int = min(
                range(len(planes)),
                key=lambda i: DistanceMeasure.calculate_min_distance_sum(group, planes[i].points)
            )
            plane_group_map[closest_plane_index] = group

        return plane_group_map

    @classmethod
    def _copy_al_points_to_rest_planes(
            cls,
            plane_group_map: dict[int, np.ndarray],
            carbon_channel: CarbonHoneycombChannel,
    ) -> Points:
        planes: list[CarbonHoneycombPlane] = carbon_channel.planes
        channel_center: np.ndarray = carbon_channel.center
        all_al_points: list[np.ndarray] = []

        for plane_i, plane in enumerate(planes):
            if plane_i in plane_group_map:
                # Just add the related points
                all_al_points.append(plane_group_map[plane_i])
                continue

            # Get the index of the related plane
            group_i: int = plane_i % len(plane_group_map)

            al_points = Points(plane_group_map[group_i])
            al_points_moved: Points = cls._move_al_to_other_plane(
                al_points, channel_center,
                al_points_plane=planes[group_i],
                target_plane=plane,
            )

            all_al_points.append(al_points_moved.points)
            continue

            angle: int = plane_i * 90  # TODO: check this logic
            al_points_rotated: Points = PointsRotator.rotate_on_angle_related_center(
                al_points_moved, angle_z=angle)

            al_points_adjusted: Points = cls._adjust_al_atoms(al_points_rotated, plane)

            all_al_points.append(al_points_adjusted.points)

        al_points_ndarray: np.ndarray = np.concatenate(all_al_points)
        return Points(al_points_ndarray)

    @staticmethod
    def _move_al_to_other_plane(
        al_points: Points,
        channel_center: np.ndarray,
        al_points_plane: CarbonHoneycombPlane,
        target_plane: CarbonHoneycombPlane,
    ) -> Points:
        """
        It does 3 steps:
        1) move al_points from al_points center to the center of the al_points_plane,
        2) move al_points from 1) step point to target_plane center,
        3) move al_points from 2) step point in the direction to the channel_center on the distance,
        that equals length of the vector from 1) step (to keep the same distance between Al and plane).
        """

        # Step 1: Move Al to `al_points_plane` center
        vector_1: np.ndarray = al_points_plane.center - al_points.center
        dist_from_plane: np.floating = np.linalg.norm(vector_1)

        # Step 2: Move from `al_points_plane.center` to `target_plane.center`
        vector_2: np.ndarray = target_plane.center - al_points_plane.center

        # Step 3: Move in the direction of `channel_center` to maintain distance
        direction_vector = channel_center - target_plane.center

        # Compute length and normalize direction vector
        direction_length: np.floating = np.linalg.norm(direction_vector)
        if direction_length == 0:
            return Points(points=al_points.points + vector_1 + vector_2)  # No need to move further

        unit_direction = direction_vector / direction_length
        scaled_vector = unit_direction * dist_from_plane  # Keep the same distance

        # Compute final transformed points
        transformed_points: np.ndarray = al_points.points + vector_1 + vector_2 + scaled_vector

        return Points(points=transformed_points)

    @classmethod
    def _adjust_al_atoms(
        cls,
        al_points: Points,
        plane: CarbonHoneycombPlane,
    ) -> Points:
        """
        Rotate a bit to find the min distance variation 
        between Al points and the plane. 
        """
        init_angle: float = 0

        result = minimize(
            cls._func_to_minimize,
            init_angle,
            args=(al_points, plane),
            method="BFGS",
            options={"disp": True},
        )

        result_angle: float = result.x

        rotated_points: Points = PointsRotator.rotate_on_angle_related_center(
            al_points, angle_z=result_angle)

        return rotated_points

    @staticmethod
    def _func_to_minimize(
        angle_z: float,
        al_points: Points,
        plane: CarbonHoneycombPlane,
    ) -> np.floating:
        rotated_points: Points = PointsRotator.rotate_on_angle_related_center(
            al_points, angle_z=angle_z)

        min_dists: np.ndarray = DistanceMeasure.calculate_min_distances(
            rotated_points.points, plane.points)

        var_result: np.floating = np.var(min_dists)

        if np.min(min_dists) < 1:
            # Avoid situations when Al atoms are too close to the plane
            return var_result * 10
        
        return var_result
