import math
import numpy as np
from numpy import ndarray, floating

from src.utils import Logger
from src.base_structure_classes import Points, CoordinateLimits
from src.data_preparation import StructureSettings
from src.coordinate_operations import (
    DistanceMeasure,
    PointsFilter,
)

from src.projects.carbon_honeycomb_actions import CarbonHoneycombChannel
from .variance_calculator import VarianceCalculator


logger = Logger("AlAtomsFilter")


class AlAtomsFilter:
    @classmethod
    def get_filtered_al_atoms(
            cls,
            carbon_channel: CarbonHoneycombChannel,
            coordinates_al: Points,
            coordinates_al_prev: Points,
            structure_settings: StructureSettings,
            max_atoms: int,
            min_dist_between_al_sum_prev: float,
            dist_and_rotation_variance_prev: float | floating,
    ) -> tuple[Points, float, float | floating, int]:

        # Prev values by default
        coordinates_al_result: Points = coordinates_al_prev
        min_dist_between_al_sum_result: float = min_dist_between_al_sum_prev
        dist_and_rotation_variance_result: float | floating = dist_and_rotation_variance_prev
        max_atoms_result: int = max_atoms

        coordinates_al_filtered: Points = cls.filter_al_atoms_related_carbon(
            coordinates_al, carbon_channel, structure_settings)

        num_of_atoms: int = len(coordinates_al_filtered)

        if num_of_atoms > max_atoms:
            coordinates_al_result = coordinates_al_filtered
            max_atoms_result = len(coordinates_al_result)

            # Update min distances between Al atoms
            min_dist_between_al_sum_result = cls._calculate_min_dist_between_al_sum(
                coordinates_al_filtered)

            # Reset dist_and_rotation_variance
            dist_and_rotation_variance_result = 0

            if num_of_atoms > 1:
                logger.info(f"Number of Al atoms: {num_of_atoms}")

                # Print average min distance between Al atoms
                ave_min_dist: float = round(min_dist_between_al_sum_result / max_atoms_result, 4)
                logger.info(f"Average min distance between Al atoms: {ave_min_dist}")

        elif num_of_atoms > 0 and num_of_atoms == max_atoms:

            # Calculate min distances between Al atoms
            current_min_dist_between_al_sum: float = cls._calculate_min_dist_between_al_sum(
                coordinates_al_filtered)

            # The nearest atoms have priority
            if current_min_dist_between_al_sum < min_dist_between_al_sum_prev:
                coordinates_al_result = coordinates_al_filtered
                max_atoms_result = len(coordinates_al_result)
                min_dist_between_al_sum_result = current_min_dist_between_al_sum

                # Print average min distance between Al atoms
                ave_min_dist: float = round(min_dist_between_al_sum_result / max_atoms_result, 4)
                logger.info(f"Average min distance between Al atoms: {ave_min_dist}")

                # Reset dist_and_rotation_variance
                dist_and_rotation_variance_result = 0

            else:
                # Check variance
                current_min_dist_sum: floating = DistanceMeasure.calculate_min_distance_sum(
                    coordinates_al_filtered.points, carbon_channel.points)
                variance_related_channel: floating = VarianceCalculator.calculate_variance_related_channel(
                    coordinates_al_filtered, carbon_channel)

                current_dist_and_rotation_variance: floating = current_min_dist_sum - variance_related_channel

                if current_dist_and_rotation_variance > dist_and_rotation_variance_prev:
                    coordinates_al_result = coordinates_al_filtered
                    max_atoms_result = len(coordinates_al_result)
                    dist_and_rotation_variance_result = current_dist_and_rotation_variance

        return (
            coordinates_al_result,
            min_dist_between_al_sum_result,
            dist_and_rotation_variance_result,
            max_atoms_result,
        )

    @staticmethod
    def _calculate_min_dist_between_al_sum(coordinates_al: Points) -> float:
        min_dist_between_al: ndarray = DistanceMeasure.calculate_min_distances_between_points(
            coordinates_al.points)
        return round(np.sum(min_dist_between_al), 4)

    # @classmethod
    # def find_max_filtered_atoms_by_minimize(
    #         cls,
    #         coordinates_carbon: ndarray,
    #         coordinates_al: ndarray,
    #         structure_settings: StructureSettings,
    # ) -> ndarray:
    #     """
    #     Parallel move and filter Al atoms related carbon atoms
    #     to find the option with the maximum Al atoms after filtering.
    #     """

    #     al_lattice_parameter: float = structure_settings.al_lattice_parameter

    #     # Objective function to optimize (to maximize number of atoms)
    #     def objective_function(params):
    #         # Unpack translation and rotation parameters
    #         step_x, step_y, angle_x, angle_y, angle_z = params

    #         # Apply translations
    #         moved_coordinates_al = coordinates_al.copy()
    #         moved_coordinates_al[:, 0] += step_x
    #         moved_coordinates_al[:, 1] += step_y

    #         # Apply rotations
    #         rotated_coordinates_al: ndarray = PointsRotator.rotate_on_angle_related_center(
    #             moved_coordinates_al, angle_x=angle_x, angle_y=angle_y, angle_z=angle_z)

    #         # Filter coordinates based on the given criteria
    #         coordinates_al_filtered: ndarray = cls.filter_al_atoms_related_carbon(
    #             rotated_coordinates_al, coordinates_carbon, structure_settings)

    #         # Calculate the negative number of atoms (since we want to maximize it)
    #         num_of_atoms: int = len(coordinates_al_filtered)

    #         # Return negative of number of atoms for minimization
    #         return -num_of_atoms

    #     # Another objective function that also considers xy variance
    #     def variance_function(params) -> float:
    #         step_x, step_y, angle_x, angle_y, angle_z = params

    #         # Apply translations
    #         moved_coordinates_al: ndarray = coordinates_al.copy()
    #         moved_coordinates_al[:, 0] += step_x
    #         moved_coordinates_al[:, 1] += step_y

    #         # Apply rotations
    #         rotated_coordinates_al: ndarray = PointsRotator.rotate_on_angle_related_center(
    #             moved_coordinates_al, angle_x=angle_x, angle_y=angle_y, angle_z=angle_z)

    #         # Filter coordinates
    #         coordinates_al_filtered: ndarray = cls.filter_al_atoms_related_carbon(
    #             rotated_coordinates_al, coordinates_carbon, structure_settings)

    #         num_of_atoms: int = len(coordinates_al_filtered)

    #         # If number of atoms is less than the maximum, discard this configuration
    #         if num_of_atoms < max_atoms:
    #             return np.inf

    #         # Otherwise, calculate xy variance
    #         xy_variance = float(
    #             np.var(coordinates_al_filtered[:, 0]) + np.var(coordinates_al_filtered[:, 1]))

    #         return -len(coordinates_al_filtered) + xy_variance

    #     # Optimization settings
    #     initial_guess: list[float] = [0.0, 0.0, 0.0, 0.0, 0.0]
    #     bounds: list[tuple[float, float]] = [
    #         (0, al_lattice_parameter),  # bounds for x translation
    #         (0, al_lattice_parameter),  # bounds for y translation
    #         (-math.pi / 4, math.pi / 2),  # bounds for x rotation
    #         (-math.pi / 4, math.pi / 2),  # bounds for y rotation
    #         (-math.pi / 4, math.pi / 2),  # bounds for z rotation
    #     ]

    #     res = minimize(objective_function, initial_guess, method='L-BFGS-B', bounds=bounds)

    #     # Extract optimal translation and rotation parameters that maximize atoms
    #     optimal_translation_and_rotation = res.x

    #     # Calculate the maximum number of atoms found (the negative of the optimized result)
    #     max_atoms = -res.fun

    #     # Now, minimize the variance using the best translation/rotation found so far
    #     res_variance = minimize(variance_function, optimal_translation_and_rotation, method='L-BFGS-B', bounds=bounds)

    #     # Extract optimized results
    #     optimized_params = res_variance.x
    #     step_x, step_y, angle_x, angle_y, angle_z = optimized_params

    #     # Apply final optimized translation and rotation
    #     final_coordinates_al = coordinates_al.copy()
    #     final_coordinates_al[:, 0] += step_x
    #     final_coordinates_al[:, 1] += step_y
    #     final_rotated_coordinates_al = PointsRotator.rotate_on_angle_related_center(
    #         final_coordinates_al, angle_x=angle_x, angle_y=angle_y, angle_z=angle_z)

    #     # Filter final coordinates
    #     coordinates_al_filtered: ndarray = cls.filter_al_atoms_related_carbon(
    #         final_rotated_coordinates_al, coordinates_carbon, structure_settings)

    #     return coordinates_al_filtered

    @classmethod
    def filter_al_atoms_related_carbon(
            cls,
            coordinates_al: Points,
            carbon_channel: CarbonHoneycombChannel,
            structure_settings: StructureSettings,
    ) -> Points:
        """Filter Al atoms related planes and then related carbon atoms."""

        coordinates_al_filtered: Points = cls._filter_atoms_related_clannel_planes(
            coordinates_al=coordinates_al,
            carbon_channel=carbon_channel,
            structure_settings=structure_settings,
        )

        limits: CoordinateLimits = carbon_channel.coordinate_limits

        coordinates_al_filtered = PointsFilter.filter_by_min_max_z(
            points_to_filter=coordinates_al_filtered,
            z_min=limits.z_min,
            z_max=limits.z_max,
            move_align_z=True)

        return cls._filter_atoms_relates_carbon_atoms(
            coordinates_al=coordinates_al_filtered,
            coordinates_carbon=carbon_channel,
            max_distance_to_carbon_atoms=structure_settings.max_distance_to_carbon_atoms)

    @staticmethod
    def _filter_atoms_related_clannel_planes(
            coordinates_al: Points,
            carbon_channel: CarbonHoneycombChannel,
            structure_settings: StructureSettings,
    ) -> Points:
        """
        Filter points from coordinates array by planes
        (remove atoms that are outside channel and atoms inside that are closer than distance_from_plane param).
        """

        filtered_coordinates: Points = coordinates_al.copy()
        carbon_channel_center: np.ndarray = carbon_channel.channel_center

        for plane in carbon_channel.planes:
            # Build plane parameters
            A, B, C, D = plane.plane_params

            direction = bool(carbon_channel_center[1] > plane.center[1])
            filtered_coordinates = PointsFilter.filter_coordinates_related_to_plane(
                filtered_coordinates,
                A, B, C, D,
                direction=direction,
                min_distance=structure_settings.distance_from_plane)

            # from src.structure_visualizer import StructureVisualizer
            # StructureVisualizer.show_two_structures(carbon_channel.points, filtered_coordinates.points)

            if len(filtered_coordinates) == 0:
                return filtered_coordinates

        return filtered_coordinates

    @staticmethod
    def _filter_atoms_relates_carbon_atoms(
        coordinates_al: Points,
        coordinates_carbon: Points,
        max_distance_to_carbon_atoms: float,
    ) -> Points:
        """
        Filter points from the coordinates_al array
        by the max distance (max_distance_to_carbon_atoms param)
        to the points in the coordinates_carbon array.
        """

        # Find the minimum distance for each atom in coordinates_al to any atom in coordinates_carbon
        min_distances: np.ndarray = DistanceMeasure.calculate_min_distances(
            coordinates_al.points, coordinates_carbon.points
        )

        filtered_al_coordinates: Points = Points(
            points=coordinates_al.points[min_distances >= max_distance_to_carbon_atoms]
        )

        return filtered_al_coordinates
