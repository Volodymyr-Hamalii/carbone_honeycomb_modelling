import numpy as np
from scipy.optimize import minimize


from src.base_structure_classes import Points
from src.coordinate_operations import (
    PointsOrganizer,
    PointsRotator,
    PointsMover,
    DistanceMeasure,
)
from src.projects.carbon_honeycomb_actions import (
    CarbonHoneycombChannel,
    CarbonHoneycombPlane,
    CarbonHoneycombUtils,
)
# from src.structure_visualizer import StructureVisualizer

from ..based_on_planes_configs import AtomsFilter


class AlAtomsTranslator:
    @classmethod
    def translate_for_all_channels(
        cls,
        coordinates_carbon: Points,
        carbon_channels: list[CarbonHoneycombChannel],
        al_channel_coordinates: Points,
    ) -> Points:

        al_final_coordinates: np.ndarray = np.empty((0, 3))
        al_final_coordinates = np.vstack((al_final_coordinates, al_channel_coordinates.points))

        al_center: np.ndarray = al_channel_coordinates.center

        for carbon_channel in carbon_channels:
            # Check that it's not the channel for which we have already translated atoms
            channel_center: np.ndarray = carbon_channel.center
            # if np.linalg.norm(channel_center - al_center) < 1:
            #     continue

            # Translate on the vector from the channel center to the al center
            vector: np.ndarray = channel_center - al_center
            al_final_coordinates = np.vstack((al_final_coordinates, al_final_coordinates + vector))

        edge_channel_centers: list[np.float32] = cls._get_centers_of_edge_carbon_channels(coordinates_carbon)

        for center in edge_channel_centers:
            vector: np.ndarray = center - al_center
            al_final_coordinates = np.vstack((al_final_coordinates, al_channel_coordinates.points + vector))

        return Points(al_final_coordinates)

    @classmethod
    def translate_for_all_planes(
        cls,
        carbon_channel: CarbonHoneycombChannel,
        al_plane_coordinates: Points,
    ) -> Points:

        planes: list[CarbonHoneycombPlane] = carbon_channel.planes

        # al_groups: list[np.ndarray] = cls._group_by_lines(al_plane_coordinates)
        plane_group_map: dict[int, np.ndarray] = cls._match_plane_with_group(
            al_plane_coordinates, planes)

        all_al_points: Points = cls._copy_al_points_to_rest_planes(plane_group_map, carbon_channel)

        al_points_upd: Points = AtomsFilter.replace_nearby_atoms_with_one_atom(all_al_points)

        return al_points_upd

    # @staticmethod
    # def _group_by_lines(al_plane_coordinates: Points) -> list[np.ndarray]:
    #     # Create groups with the same X and Y coordinate
    #     # (like, to split all points into columns)
    #     groups_by_xy: dict[
    #         tuple[np.float32, np.float32], np.ndarray
    #     ] = PointsOrganizer.group_by_unique_xy(al_plane_coordinates.points)

    #     # Define groups that lie on the same line
    #     groups_by_the_xy_lines: list[
    #         dict[tuple[np.float32, np.float32], np.ndarray]
    #     ] = PointsOrganizer.group_by_the_xy_lines(groups_by_xy, epsilon=1e-3, min_points_in_line=3)

    #     grouped_by_lines: list[np.ndarray] = []

    #     for i in groups_by_the_xy_lines:
    #         group: np.ndarray = np.concatenate(list(i.values()))
    #         grouped_by_lines.append(group)

    #     return grouped_by_lines

    @staticmethod
    def _match_plane_with_group(
        al_plane_coordinates: Points,
        planes: list[CarbonHoneycombPlane],
    ) -> dict[int, np.ndarray]:
        """ Returns the dict with the plane index and corresponding Al points. """

        plane_group_map: dict[int, list[np.ndarray]] = {
            plane_i: [] for plane_i in range(len(planes))
        }

        for atom in al_plane_coordinates.points:
            min_dist: float = np.inf
            plane_indexes: list[int] = []

            for plane_i, plane in enumerate(planes):
                min_distances: np.ndarray = DistanceMeasure.calculate_min_distances(np.array([atom]), plane.points)
                min_dist_plane: float = round(min(min_distances), 2)

                if min_dist_plane < min_dist:
                    min_dist: float = min_dist_plane
                    plane_indexes = [plane_i]

                elif min_dist == min_dist_plane:
                    plane_indexes.append(plane_i)

            for plane_i in plane_indexes:
                plane_group_map[plane_i].append(atom)
            # Find the plane with the minimum distance sum
            # closest_plane_index: float = min(
            #     range(len(planes)),
            #     key=lambda i: DistanceMeasure.calculate_min_distances(np.array([atom]), planes[i].points)
            # )
            # plane_group_map[closest_plane_index] = atom

        max_len: int = max(len(value) for value in plane_group_map.values())

        # Filter plane groups with the max length
        plane_group_map = {
            key: value
            for key, value in plane_group_map.items()
            if len(value) == max_len
        }

        return {
            key: np.array(value)
            for key, value
            in plane_group_map.items()
        }

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

            # all_al_points.append(al_points_moved.points)
            # continue

            angle: float = (group_i - plane_i) * np.pi / 3

            al_points_rotated: Points = PointsRotator.rotate_on_angle_related_center(
                al_points_moved, angle_z=angle)

            al_points_adjusted: Points = cls._adjust_al_atoms(al_points_rotated, plane)

            # Check if we need to reflect points
            al_points_reflected: Points = PointsMover.reflect_through_vertical_axis(al_points_adjusted)

            al_points_reflected_min_dists: np.floating = DistanceMeasure.calculate_min_distance_sum(
                al_points_reflected.points, plane.points)
            al_points_adjusted_min_dists: np.floating = DistanceMeasure.calculate_min_distance_sum(
                al_points_adjusted.points, plane.points)

            if al_points_adjusted_min_dists > al_points_reflected_min_dists:
                all_al_points.append(al_points_adjusted.points)
            else:
                all_al_points.append(al_points_reflected.points)

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
        init_angle: np.ndarray = np.array([0.0])

        result = minimize(
            cls._func_to_minimize,
            init_angle,
            args=(al_points, plane),
            method="Powell",
            options={"disp": True},
        )

        result_angle: float = result.x.item()

        rotated_points: Points = PointsRotator.rotate_on_angle_related_center(
            al_points, angle_z=result_angle)

        return rotated_points

    @staticmethod
    def _func_to_minimize(
        angle_z: np.ndarray,
        al_points: Points,
        plane: CarbonHoneycombPlane,
    ) -> np.floating:
        rotated_points: Points = PointsRotator.rotate_on_angle_related_center(
            al_points, angle_z=angle_z[0])

        min_dists: np.ndarray = DistanceMeasure.calculate_min_distances(
            rotated_points.points, plane.points)

        var_result: np.floating = np.var(min_dists)

        if np.any(min_dists < 1.0):
            # Avoid situations when Al atoms are too close to the plane
            return var_result * 10

        return var_result

    @staticmethod
    def _get_centers_of_edge_carbon_channels(
        coordinates_carbon: Points,
    ) -> list[np.float32]:
        """ 
        Takes the lines with the pointns with the same Y coordinate, 
        takes left and right (with the min and max X coordinate) parallel segments from the lines
        and returns the centers of the segment pairs.
        """

        # Get the lines with the points with the same Y coordinate
        groups_by_the_xy_lines: list[dict[tuple[np.float32, np.float32], np.ndarray]] = PointsOrganizer.group_by_the_xy_lines(
            PointsOrganizer.group_by_unique_xy(coordinates_carbon.points), epsilon=1e-1, min_points_in_line=3)

        honeycomb_planes_groups: list[
            dict[tuple[np.float32, np.float32], np.ndarray]
        ] = CarbonHoneycombUtils.split_xy_groups_by_max_distances(
            groups_by_the_xy_lines, max_distance_between_xy_groups=3)

        # Filter honeycomb_planes_groups groups with the same Y in the group
        honeycomb_planes_groups = [
            group for group in honeycomb_planes_groups
            if list(group.keys())[0][1] == list(group.keys())[1][1]
        ]

        # Keep only planes, that contains max or min X coordinate along the all groups
        max_x: np.floating = max(
            subgroup[0]
            for group in honeycomb_planes_groups
            for subgroup in group.keys()
        )
        min_x: np.floating = min(
            subgroup[0]
            for group in honeycomb_planes_groups
            for subgroup in group.keys()
        )

        honeycomb_planes_groups_left = [
            group for group in honeycomb_planes_groups
            if any(key[0] in [min_x] for key in group.keys())
        ]

        honeycomb_planes_groups_right = [
            group for group in honeycomb_planes_groups
            if any(key[0] in [max_x] for key in group.keys())
        ]

        honeycomb_planes_groups_left_coordinates = [
            point for group in honeycomb_planes_groups_left
            for point in group.values()
        ]

        honeycomb_planes_groups_right_coordinates = [
            point for group in honeycomb_planes_groups_right
            for point in group.values()
        ]

        left_centers = np.mean(np.concatenate(honeycomb_planes_groups_left_coordinates), axis=0)
        right_centers = np.mean(np.concatenate(honeycomb_planes_groups_right_coordinates), axis=0)

        return [left_centers, right_centers]
