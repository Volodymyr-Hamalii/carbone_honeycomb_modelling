import numpy as np

from src.utils import Logger, Constants
from src.coordinate_operations import DistanceMeasure
from src.base_structure_classes import Points
from src.projects.carbon_honeycomb_actions import CarbonHoneycombChannel


logger = Logger("AtomConfigurator")


class AtomConfigurator:
    @classmethod
    def space_al_atoms_equidistantly(
        cls,
        coordinates_al: Points,
        carbon_channel: CarbonHoneycombChannel,
    ) -> Points:
        """ 
        To space Al atoms equally apart from each other
        to keep the minimal distance between atoms equals dist_between_al_atoms.

        Algorithm: iterate all Al points until the cycle is not stopped:
        1. Calculate the distance matrix between the Al points.
        2. Get the Al points that have a dist to the closest neighbor < min_dist_between_al_atoms. 
            If there are no such points -- stop the cycle.
        3. Find the point that have the minimal average distance to the neighbors.
        4. Move the point in the direction to the channel center to have the distance
            to the closest neighbor == dist_between_al_atoms.
            If there is no movement -- stop the cycle.
        """
        # TODO: check other possible other algorithm:
        # 1. Split the Al points by groups with the closest honeycomb plane.
        # 2. Iterate the plane groups:
        #     2.1. Build dist matrix between Al points.
        #     2.2. Get the Al point that have a dist to the closest neighbor < min_dist_between_al_atoms. If there are no such points -- stop the cycle.
        #     2.3. Find the point that have the minimal average distance to the neighbors
        #     2.4. Move the point in the direction from a plane. If there is no movement -- stop the cycle.

        al_points: np.ndarray = coordinates_al.points

        while True:
            dist_matrix: np.ndarray = DistanceMeasure.calculate_dist_matrix(al_points)
            point_index: np.intp | None = cls._find_the_point_with_min_ave_dist(dist_matrix, min_dist_between_al_atoms)

            if point_index is None:
                break

            al_points_upd: np.ndarray | None = cls._upd_points_with_moved_point(
                al_points=al_points,
                point_index=point_index,
                distances=dist_matrix[point_index],
                carbon_channel_center=carbon_channel.channel_center,
                min_dist_between_al_atoms=min_dist_between_al_atoms,
            )

            if al_points_upd is None:
                break

            al_points = al_points_upd

        return Points(points=al_points)

    @staticmethod
    def _find_the_point_with_min_ave_dist(
        dist_matrix: np.ndarray,
        min_dist_between_al_atoms: float
    ) -> np.intp | None:
        """
        Returns the index of the point with the minimal average distance to the neighbors.
        If there is no such index -- returns None.
        """
        try:
            # Create a mask for values less than `min_dist_between_al_atoms`
            mask: np.ndarray = dist_matrix < min_dist_between_al_atoms

            # Replace values not satisfying the condition with `NaN`
            filtered_matrix: np.ndarray = np.where(mask, dist_matrix, np.nan)

            # Compute the average for each row, ignoring NaN values
            min_ave_dists: np.ndarray = np.nanmean(filtered_matrix, axis=1)

            # Find the lowest average distance
            return np.nanargmin(min_ave_dists)

        except ValueError:
            return

    @staticmethod
    def _upd_points_with_moved_point(
        al_points: np.ndarray,
        point_index: np.intp,
        distances: np.ndarray,
        min_dist_between_al_atoms: float,
        carbon_channel_center: np.ndarray,
    ) -> np.ndarray | None:
        """ 
        Move the point from al_point with point_index in the direction to the carbon channel main axis
        to set the minimal distance to the closest neighbor from the al_point equals min_dist_between_al_atoms.

        If such move is not possible, returns None.
        """

        point_to_move: np.ndarray = al_points[point_index]

        # Find the point on the channel axis to which we can move the point
        point_on_main_axis: np.ndarray = np.array([
            carbon_channel_center[0],
            carbon_channel_center[1],
            point_to_move[2],
        ])

        # Vector from the point_to_move to the channel axis
        direction_vector = point_on_main_axis - point_to_move
        direction_vector_length: np.floating = np.linalg.norm(direction_vector)

        # If the direction vector length is zero, the point is already on the axis
        if direction_vector_length == 0:
            return None

        # Normalize the direction vector
        direction_unit_vector = direction_vector / direction_vector_length

        # Find the closest neighbor distance
        closest_neighbor_index: np.intp = np.argmin(distances)
        closest_distance = distances[closest_neighbor_index]

        # Calculate the required distance to move
        distance_to_move = min_dist_between_al_atoms - closest_distance

        # If the closest neighbor already satisfies the minimum distance, no move is needed
        if distance_to_move <= 0:
            return None

        # Move the point along the direction vector by the required distance
        moved_point: np.ndarray = point_to_move + distance_to_move * direction_unit_vector

        # Validate that the new point satisfies the required distance
        new_distances: np.ndarray = np.linalg.norm(al_points - moved_point, axis=1)
        if np.any(new_distances < min_dist_between_al_atoms):
            return None

        # Replace the old point with the moved point
        updated_points: np.ndarray = al_points.copy()
        updated_points[point_index] = moved_point

        return updated_points
