from dataclasses import dataclass
from typing import Any, Generator
import numpy as np

from src.utils import Logger, execution_time_logger, Constants
from src.coordinate_operations import DistanceMeasure
from src.base_structure_classes import Points
from src.projects.carbon_honeycomb_actions import CarbonHoneycombChannel

from .intercalation_exceptions import NotOptimalStructureError
from .atoms_filter import AtomsFilter

logger = Logger("AtomConfigurator")


@dataclass(frozen=True)
class ConfigParams:
    max_points_to_move_before_reset_coef: float
    num_of_points_to_skip: int
    percent_to_remove: float

    def __repr__(self) -> str:
        params_dict: dict = asdict(self)
        return ', '.join(f"{key}: {value}" for key, value in params_dict.items())


class AtomConfigurator:
    @staticmethod
    def _iterate_params_for_reorganizing_al_atoms() -> Generator[ConfigParams, Any, None]:
        # max_points_to_move_before_reset_coef: list[float] = [1.1, 1.3, 1.5, 1.7, 2]
        max_points_to_move_before_reset_coef_list: list[float] = [1, 0.9, 0.75, 0.5]
        num_of_points_to_skip_list: list[int] = [0, 1, 2, 3]
        percent_to_remove_params: list[float] = [1, 0.9, 0.8, 0.7]

        for max_points_to_move_before_reset_coef in max_points_to_move_before_reset_coef_list:
            for num_of_points_to_skip in num_of_points_to_skip_list:
                for percent_to_remove in percent_to_remove_params:
                    yield ConfigParams(
                        max_points_to_move_before_reset_coef=max_points_to_move_before_reset_coef,
                        num_of_points_to_skip=num_of_points_to_skip,
                        percent_to_remove=percent_to_remove,
                    )

    @classmethod
    @execution_time_logger
    def reorganize_al_atoms(cls, coordinates_al: Points, carbon_channel: CarbonHoneycombChannel) -> Points:
        """
        Move Al atoms to have the correct minimal distance between them
        (or remove them if we can't achive optimal configuration just by moving).
        """

        logger.info(f"Number of atoms before reorganizing: {len(coordinates_al.points)}")

        # # Take a step as 1% of the distance between atoms
        # step: float = Constants.phys.al.DIST_BETWEEN_ATOMS * 0.05

        # start: float = Constants.phys.al.DIST_BETWEEN_ATOMS * 1.2
        # end: float = Constants.phys.al.DIST_BETWEEN_ATOMS * 1.21

        # best_result: Points = coordinates_al
        # best_result_len: int = 0

        # for min_dist in np.arange(start, end, step):
        #     logger.info(f"Min dist: {min_dist}")

        #     coordinates_al_upd: Points = AtomsFilter.remove_some_close_atoms(coordinates_al, min_dist)

        #     if len(coordinates_al_upd.points) > best_result_len:
        #         best_result = coordinates_al_upd
        #         best_result_len = len(coordinates_al_upd.points)

        #     elif len(coordinates_al_upd.points) == len(coordinates_al.points):
        #         break

        # Filter and move Al atoms to have the correct minimal distance between them
        dist_between_al_atoms: float = Constants.phys.al.DIST_BETWEEN_ATOMS
        min_dist: float = dist_between_al_atoms * 1.2
        step: float = dist_between_al_atoms * 0.025

        max_neighbours_start: int = -1
        max_neighbours: int = -1

        # Bool to rotate between increasing min_dist and decreasing max_neighbours
        # to_change_step: bool = False

        best_result: Points = Points(points=np.array([]))
        min_distances: np.ndarray = np.array([])

        for config_params in cls._iterate_params_for_reorganizing_al_atoms():
            logger.info(f"Max points to move before reset coef: {config_params.max_points_to_move_before_reset_coef}")

            while True:
                try:
                    coordinates_al_upd, max_neighbours = cls.space_al_atoms_equidistantly(
                        coordinates_al,
                        carbon_channel,
                        min_dist,
                        max_neighbours,
                        config_params,
                    )

                    if max_neighbours > max_neighbours_start:
                        max_neighbours_start = max_neighbours

                    # logger.info(f"Number of atoms after reorganizing: {len(coordinates_al.points)}")
                    if len(coordinates_al_upd.points) > len(best_result.points):
                        best_result = coordinates_al_upd
                        min_distances = DistanceMeasure.calculate_min_distances(
                            carbon_channel.points, coordinates_al_upd.points)

                    elif len(coordinates_al.points) == len(coordinates_al_upd.points):
                        # Get the set with the minimal sum of distances to carbon atoms
                        min_distances_upd: np.ndarray = DistanceMeasure.calculate_min_distances(
                            carbon_channel.points, coordinates_al_upd.points)

                        if np.sum(min_distances_upd) < np.sum(min_distances):
                            best_result = coordinates_al_upd
                            min_distances = min_distances_upd

                except NotOptimalStructureError as e:
                    logger.warning("Not optimal structure:", e)

                if min_dist < dist_between_al_atoms:
                    # msg: str = "Min dist is less than the distance between atoms, Al building failed."
                    # raise AlBuildError(msg)
                    break

                if max_neighbours <= max_neighbours_start // 2:
                    # Reset max_neighbours
                    max_neighbours = -1
                    max_neighbours_start = -1
                    min_dist -= step
                else:
                    max_neighbours -= 1

        return best_result

    @classmethod
    # @execution_time_logger
    def space_al_atoms_equidistantly(
        cls,
        coordinates_al: Points,
        carbon_channel: CarbonHoneycombChannel,
        min_dist: float,
        max_neighbours: int,
        config_params: ConfigParams,
    ) -> tuple[Points, int]:
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

        max_points_to_move_before_reset_coef: float = config_params.max_points_to_move_before_reset_coef
        num_of_points_to_skip: float = config_params.num_of_points_to_skip
        percent_to_remove: float = config_params.percent_to_remove

        al_points: np.ndarray = coordinates_al.points

        dist_between_al_atoms: float = Constants.phys.al.DIST_BETWEEN_ATOMS
        min_dist_between_al_atoms: float = Constants.phys.al.MIN_RECOMENDED_DIST_BETWEEN_ATOMS
        moved_point_indexes: set[np.intp] = set()

        counter = 0
        max_counter: int = len(coordinates_al) * 15

        logger.info(f"Min dist: {min_dist}; max neighbours: {max_neighbours}")

        result: tuple[Points, int] = AtomsFilter.remove_some_close_atoms(
            coordinates_al, min_dist, max_neighbours, num_of_points_to_skip=0)

        coordinates_al_upd: Points = result[0]
        max_neighbours = result[1]
        logger.info(f"Number of atoms after filtering to reorganize: {len(coordinates_al_upd.points)}")

        # if max_neighbours > max_neighbours_start:
        #     max_neighbours_start = max_neighbours

        # For debugging purposes
        # if len(coordinates_al_upd) <= 30:
        #     pass
        max_points_to_move_before_reset: int = round(
            len(coordinates_al_upd) * max_points_to_move_before_reset_coef)

        while True:
            dist_matrix: np.ndarray = DistanceMeasure.calculate_dist_matrix(al_points)
            point_index: np.intp | None = cls._find_the_point_with_min_ave_dist(
                dist_matrix,
                min_dist_between_al_atoms,
                moved_point_indexes,
                num_of_points_to_skip,
            )

            if point_index is None:
                # Check if all points are already spaced correctly

                # Reset moved_point_indexes
                moved_point_indexes = set()

                point_index = cls._find_the_point_with_min_ave_dist(
                    dist_matrix,
                    min_dist_between_al_atoms,
                    moved_point_indexes,
                    num_of_points_to_skip,
                )

                if point_index is None:
                    # if len(DistanceMeasure.
                    #        calculate_min_distances_between_points(al_points)
                    #        [DistanceMeasure.calculate_min_distances_between_points(al_points) <
                    #         min_dist_between_al_atoms]) > 0:
                    #     pass
                    break

            al_points_upd: np.ndarray | None = cls._upd_points_with_moved_point(
                al_points=al_points,
                point_index=point_index,
                distances=dist_matrix[point_index],
                carbon_channel_center=carbon_channel.channel_center,
                dist_between_al_atoms=dist_between_al_atoms,
                min_dist_between_al_atoms=min_dist_between_al_atoms,
            )

            if al_points_upd is None:
                # if len(
                #         DistanceMeasure.
                #         calculate_min_distances_between_points(al_points)
                #         [DistanceMeasure.calculate_min_distances_between_points(al_points) < min_dist_between_al_atoms]) > 0:
                #     pass
                break

            dists_between_al_and_carbon: np.ndarray = DistanceMeasure.calculate_min_distances(
                al_points_upd, carbon_channel.points)

            # Filter case when Al atoms are too far from carbon atoms
            k = 1.15  # Max dist from Al to carbon atoms coef related dist_between_al_atoms
            too_far_indices = np.where(dists_between_al_and_carbon > dist_between_al_atoms * k)[0]

            if len(too_far_indices) > 0:
                al_points_upd = np.delete(al_points_upd, too_far_indices, axis=0)
                moved_point_indexes = set()
            else:
                moved_point_indexes.add(point_index)

            al_points = al_points_upd

            if len(moved_point_indexes) >= max_points_to_move_before_reset:
                # Reset moved_point_indexes
                moved_point_indexes = set()

            counter += 1
            if counter == max_counter:
                raise NotOptimalStructureError(f"Too many iterations to space Al atoms equidistantly ({max_counter}).")

        logger.info(f"Number of iterations in space_al_atoms_equidistantly: {counter}/{max_counter}")

        result_al_points: Points = AtomsFilter.remove_too_close_atoms(
            Points(points=al_points), min_allowed_dist=min_dist_between_al_atoms)

        return result_al_points, max_neighbours

    @classmethod
    def _find_the_point_with_min_ave_dist(
        cls,
        dist_matrix: np.ndarray,
        min_dist_between_al_atoms: float,
        moved_point_indexes: set[np.intp],
        num_of_points_to_skip: int = 0,
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

            num_of_too_close_neighbours: np.ndarray = np.sum(~np.isnan(filtered_matrix), axis=1)

            # Find the maximum number of too-close neighbors
            max_too_close_neighbours = np.max(num_of_too_close_neighbours)

            if max_too_close_neighbours == 0:
                return

            # Filter rows where the number of too-close neighbors matches the maximum
            candidate_indexes = np.where(num_of_too_close_neighbours == max_too_close_neighbours)[0]

            # Compute the average distance for the filtered rows
            min_ave_dists: np.ndarray = np.nanmean(filtered_matrix[candidate_indexes], axis=1)

            if num_of_points_to_skip == 0:
                # Find the index of the candidate with the lowest average distance
                candidate_index: np.intp = np.nanargmin(min_ave_dists)
                global_index: np.intp = candidate_indexes[candidate_index]

            else:
                # Skip first `num_of_points_to_skip` points
                # Sort candidates by their average distance
                sorted_candidates: np.ndarray = np.argsort(min_ave_dists)

                # Ensure we don't try to skip more points than available
                if num_of_points_to_skip >= len(sorted_candidates):
                    return None  # No valid points remain after skipping

                # Get the index after skipping
                candidate_index: np.intp = sorted_candidates[num_of_points_to_skip]
                global_index: np.intp = candidate_indexes[candidate_index]

            # Check if the selected point has already been moved
            if global_index in moved_point_indexes:
                # Create a copy of the distance matrix to exclude the current point
                dist_matrix_copy: np.ndarray = dist_matrix.copy()
                dist_matrix_copy[global_index] = np.inf  # Effectively exclude this point
                return cls._find_the_point_with_min_ave_dist(
                    dist_matrix_copy, min_dist_between_al_atoms, moved_point_indexes,
                    num_of_points_to_skip=0,
                )

            return global_index

        except ValueError:
            return

    @staticmethod
    def _upd_points_with_moved_point(
        al_points: np.ndarray,
        point_index: np.intp,
        distances: np.ndarray,
        carbon_channel_center: np.ndarray,
        dist_between_al_atoms: float,
        min_dist_between_al_atoms: float,
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
        direction_vector_length: np.float64 = np.linalg.norm(direction_vector)

        # If the direction vector length is zero, the point is already on the axis
        if direction_vector_length == 0:
            return None

        # Normalize the direction vector
        direction_unit_vector = direction_vector / direction_vector_length

        # Find the closest neighbor distance
        closest_neighbor_index: np.intp = np.argmin(distances)
        closest_distance: np.float64 = distances[closest_neighbor_index]

        # If the closest neighbor already satisfies the minimum distance, no move is needed
        if round(closest_distance, 2) >= round(min_dist_between_al_atoms, 2):
            return None

        # Project the closest neighbor onto the line formed by point_to_move and point_on_main_axis
        closest_neighbor: np.ndarray = al_points[closest_neighbor_index]
        vector_to_closest: np.ndarray = closest_neighbor - point_to_move
        projection_length: np.float64 = np.dot(vector_to_closest, direction_unit_vector)
        projection_point: np.ndarray = point_to_move + projection_length * direction_unit_vector

        # Adjust the position along the line to achieve the desired distance
        total_distance_to_move: np.float64 | float = np.sqrt(dist_between_al_atoms**2 - closest_distance**2)

        max_dist_to_move: float = DistanceMeasure.calculate_distance_between_2_points(
            point_to_move, point_on_main_axis
        )

        # Move no more than to the center of the channel
        if total_distance_to_move > max_dist_to_move:
            total_distance_to_move = max_dist_to_move

        moved_point: np.ndarray = projection_point + total_distance_to_move * direction_unit_vector

        # if total_distance_to_move < 1:
        #     logger.info(f"Total distance to move: {total_distance_to_move}")

        # Validate that the new point satisfies the required distance
        # new_distances: np.ndarray = np.linalg.norm(al_points - moved_point, axis=1)
        # if np.any(new_distances < min_dist_between_al_atoms):
        #     return None

        # Replace the old point with the moved point
        updated_points: np.ndarray = al_points.copy()

        # dists_before = DistanceMeasure.calculate_min_distances_between_points(updated_points)
        updated_points[point_index] = moved_point
        # dists_after = DistanceMeasure.calculate_min_distances_between_points(updated_points)

        return updated_points
