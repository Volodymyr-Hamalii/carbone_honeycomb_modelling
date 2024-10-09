import math
import numpy as np
from numpy import ndarray, floating
from scipy.spatial.distance import cdist
from scipy.optimize import minimize

from ..utils import Logger, execution_time_logger
from ..data_preparation import StructureSettings, ChannelPoints
from ..coordinates_actions import (
    PlanesBuilder,
    CoordinatesFilter,
    StructureRotator,
)
from ..calculators import VarianceCalculator, DistanceCalculator

logger = Logger("AlAtomsFilter")


class AlAtomsFilter:

    @classmethod
    @execution_time_logger
    def find_max_filtered_atoms(
            cls,
            coordinates_carbone: ndarray,
            coordinates_al: ndarray,
            structure_settings: StructureSettings,
    ) -> ndarray:
        """
        Parallel move and filter Al atoms related carbone atoms
        to find the option with the maximum Al atoms after filtering.
        """

        max_atoms: int = 0
        min_dist_between_al_sum: float = np.inf
        dist_and_rotation_variance: float | floating = 0

        coordinates_al_result: ndarray = coordinates_al.copy()

        al_lattice_parameter: float = structure_settings.al_lattice_parameter

        range_to_move_step: float = al_lattice_parameter / 20
        range_to_move: ndarray = np.arange(0, al_lattice_parameter, range_to_move_step)

        angle_range_to_rotate_step: float = math.pi / 30
        angle_range_to_rotate: ndarray = np.arange(0, math.pi / 3, angle_range_to_rotate_step)

        for step_x in range_to_move:
            moved_x_coordinates_al: ndarray = coordinates_al.copy()
            moved_x_coordinates_al[:, 0] += step_x

            for step_y in range_to_move:
                moved_xy_coordinates_al: ndarray = moved_x_coordinates_al.copy()
                moved_xy_coordinates_al[:, 1] += step_y

                for angle_x in angle_range_to_rotate:
                    x_rotaded_coordinates_al: ndarray = StructureRotator.rotate_on_angle_related_center(
                        moved_xy_coordinates_al.copy(), angle_x=angle_x)

                    for angle_y in angle_range_to_rotate:
                        xy_rotaded_coordinates_al: ndarray = StructureRotator.rotate_on_angle_related_center(
                            x_rotaded_coordinates_al.copy(), angle_y=angle_y)

                        for angle_z in angle_range_to_rotate:
                            xyz_rotaded_coordinates_al: ndarray = StructureRotator.rotate_on_angle_related_center(
                                xy_rotaded_coordinates_al.copy(), angle_z=angle_z)
                            
                            result: tuple = cls._get_filtered_al_atoms(
                                coordinates_carbone=coordinates_carbone,
                                coordinates_al=xyz_rotaded_coordinates_al,
                                structure_settings=structure_settings,
                                coordinates_al_prev=coordinates_al_result,
                                max_atoms=max_atoms,
                                min_dist_between_al_sum_prev=min_dist_between_al_sum,
                                dist_and_rotation_variance_prev=dist_and_rotation_variance,
                            )

                            coordinates_al_result = result[0]
                            min_dist_between_al_sum = result[1]
                            dist_and_rotation_variance = result[2]

                            max_atoms = len(coordinates_al_result)

        return coordinates_al_result

    @classmethod
    def _get_filtered_al_atoms(
            cls,
            coordinates_carbone: ndarray,
            coordinates_al: ndarray,
            structure_settings: StructureSettings,
            coordinates_al_prev: ndarray,
            max_atoms: int,
            min_dist_between_al_sum_prev: float,
            dist_and_rotation_variance_prev: float | floating,
    ) -> tuple[ndarray, float, float | floating]:

        # Prev values by default
        coordinates_al_result: ndarray = coordinates_al_prev
        min_dist_between_al_sum_result: float = min_dist_between_al_sum_prev
        dist_and_rotation_variance_result: float | floating = dist_and_rotation_variance_prev

        coordinates_al_filtered: ndarray = cls.filter_al_atoms_related_carbone(
            coordinates_al, coordinates_carbone, structure_settings)

        num_of_atoms: int = len(coordinates_al_filtered)

        if num_of_atoms > max_atoms:
            logger.info("max_atoms", max_atoms)

            coordinates_al_result = coordinates_al_filtered

            # Update min distances between Al atoms
            min_dist_between_al_sum_result = cls._calculate_min_dist_between_al_sum(
                coordinates_al_filtered)
            logger.info("min_dist_between_al_sum", min_dist_between_al_sum_result)

            # Reset dist_and_rotation_variance
            dist_and_rotation_variance_result = 0

        elif num_of_atoms > 0 and num_of_atoms == max_atoms:
            # The advantage is the maximum equidistant points from the channel
            # current_min_distances_sum: floating = PointsOrganizer.calculate_min_distance_sum(
            #     coordinates_al_filtered, coordinates_carbone)

            # Calculate min distances between Al atoms
            current_min_dist_between_al_sum: float = cls._calculate_min_dist_between_al_sum(
                coordinates_al_filtered)

            # The nearest atoms have priority
            if current_min_dist_between_al_sum < min_dist_between_al_sum_prev:
                coordinates_al_result = coordinates_al_filtered
                min_dist_between_al_sum_result = current_min_dist_between_al_sum

                # Reset dist_and_rotation_variance
                dist_and_rotation_variance_result = 0

                logger.info("min_dist_between_al_sum", min_dist_between_al_sum_result)
            else:
                # Check variance
                current_min_dist_sum: floating = DistanceCalculator.calculate_min_distance_sum(
                    coordinates_al_filtered, coordinates_carbone)
                variance_related_channel: floating = VarianceCalculator.calculate_variance_related_channel(
                    coordinates_al_filtered, coordinates_carbone)

                current_dist_and_rotation_variance: floating = current_min_dist_sum - variance_related_channel

                if current_dist_and_rotation_variance > dist_and_rotation_variance_prev:
                    coordinates_al_result = coordinates_al_filtered
                    dist_and_rotation_variance_result = current_dist_and_rotation_variance

        return (
            coordinates_al_result,
            min_dist_between_al_sum_result,
            dist_and_rotation_variance_result,
        )

    @staticmethod
    def _calculate_min_dist_between_al_sum(coordinates_al: ndarray) -> float:
        min_dist_between_al: ndarray = DistanceCalculator.calculate_min_distances_between_points(
            coordinates_al)
        return round(np.sum(min_dist_between_al), 4)

    @classmethod
    def find_max_filtered_atoms_by_minimize(
            cls,
            coordinates_carbone: ndarray,
            coordinates_al: ndarray,
            structure_settings: StructureSettings,
    ) -> ndarray:
        """
        Parallel move and filter Al atoms related carbone atoms
        to find the option with the maximum Al atoms after filtering.
        """

        al_lattice_parameter: float = structure_settings.al_lattice_parameter

        # Objective function to optimize (to maximize number of atoms)
        def objective_function(params):
            # Unpack translation and rotation parameters
            step_x, step_y, angle_x, angle_y, angle_z = params

            # Apply translations
            moved_coordinates_al = coordinates_al.copy()
            moved_coordinates_al[:, 0] += step_x
            moved_coordinates_al[:, 1] += step_y

            # Apply rotations
            rotated_coordinates_al = StructureRotator.rotate_on_angle_related_center(
                moved_coordinates_al, angle_x=angle_x, angle_y=angle_y, angle_z=angle_z)

            # Filter coordinates based on the given criteria
            coordinates_al_filtered: ndarray = cls.filter_al_atoms_related_carbone(
                rotated_coordinates_al, coordinates_carbone, structure_settings)

            # Calculate the negative number of atoms (since we want to maximize it)
            num_of_atoms = len(coordinates_al_filtered)

            # Return negative of number of atoms for minimization
            return -num_of_atoms

        # Another objective function that also considers xy variance
        def variance_function(params):
            # Unpack translation and rotation parameters
            step_x, step_y, angle_x, angle_y, angle_z = params

            # Apply translations
            moved_coordinates_al = coordinates_al.copy()
            moved_coordinates_al[:, 0] += step_x
            moved_coordinates_al[:, 1] += step_y

            # Apply rotations
            rotated_coordinates_al = StructureRotator.rotate_on_angle_related_center(
                moved_coordinates_al, angle_x=angle_x, angle_y=angle_y, angle_z=angle_z)

            # Filter coordinates based on the given criteria
            coordinates_al_filtered: ndarray = cls.filter_al_atoms_related_carbone(
                rotated_coordinates_al, coordinates_carbone, structure_settings)

            # Count the number of atoms
            num_of_atoms = len(coordinates_al_filtered)

            # If number of atoms is less than the maximum, discard this configuration
            if num_of_atoms < max_atoms:
                return np.inf  # Higher value to discourage this solution

            # Otherwise, calculate xy variance
            xy_variance = float(
                np.var(coordinates_al_filtered[:, 0]) + np.var(coordinates_al_filtered[:, 1]))

            return -len(coordinates_al_filtered) + xy_variance

        # Optimization settings
        initial_guess = [0.0, 0.0, 0.0, 0.0, 0.0]
        bounds = [
            (0, al_lattice_parameter),  # bounds for x translation
            (0, al_lattice_parameter),  # bounds for y translation
            (-math.pi / 4, math.pi / 2),  # bounds for x rotation
            (-math.pi / 4, math.pi / 2),  # bounds for y rotation
            (-math.pi / 4, math.pi / 2),  # bounds for z rotation
        ]

        # First, maximize the number of atoms
        res = minimize(objective_function, initial_guess, method='L-BFGS-B', bounds=bounds)

        # Extract optimal translation and rotation parameters that maximize atoms
        optimal_translation_and_rotation = res.x

        # Calculate the maximum number of atoms found (the negative of the optimized result)
        max_atoms = -res.fun

        # Now, minimize the variance using the best translation/rotation found so far
        res_variance = minimize(variance_function, optimal_translation_and_rotation, method='L-BFGS-B', bounds=bounds)

        # Extract optimized results
        optimized_params = res_variance.x
        step_x, step_y, angle_x, angle_y, angle_z = optimized_params

        # Apply final optimized translation and rotation
        final_coordinates_al = coordinates_al.copy()
        final_coordinates_al[:, 0] += step_x
        final_coordinates_al[:, 1] += step_y
        final_rotated_coordinates_al = StructureRotator.rotate_on_angle_related_center(
            final_coordinates_al, angle_x=angle_x, angle_y=angle_y, angle_z=angle_z)

        # Filter final coordinates
        coordinates_al_filtered = cls.filter_al_atoms_related_carbone(
            final_rotated_coordinates_al, coordinates_carbone, structure_settings)

        return coordinates_al_filtered

    @classmethod
    def filter_al_atoms_related_carbone(
            cls,
            coordinates_al: ndarray,
            coordinates_carbone: ndarray,
            structure_settings: StructureSettings
    ) -> ndarray:
        """Filter Al atoms related planes and then related carbone atoms."""

        coordinates_al_filtered: ndarray = cls._filter_atoms_related_clannel_planes(
            coordinates=coordinates_al,
            points_to_set_channel_planes=structure_settings.points_to_set_channel_planes)

        coordinates_al_filtered = CoordinatesFilter.filter_by_min_max_z(
            coordinates_to_filter=coordinates_al_filtered,
            coordinates_with_min_max_z=coordinates_carbone,
            move_align_z=True)

        return cls._filter_atoms_relates_carbone_atoms(
            coordinates_al=coordinates_al_filtered,
            coordinates_carbone=coordinates_carbone,
            max_distance_to_carbone_atoms=structure_settings.al_lattice_parameter / 2)

    @staticmethod
    def _filter_atoms_related_clannel_planes(
            coordinates: ndarray,
            points_to_set_channel_planes: list[ChannelPoints],
            distance_from_plane: float = 0,
    ) -> ndarray:
        """
        Filter points from coordinates array by planes
        (remove atoms that are outside channel and atoms inside that are closer than distance_from_plane param).
        """

        filtered_coordinates: ndarray = coordinates

        for plane_data in points_to_set_channel_planes:
            # Build plane parameters
            A, B, C, D = PlanesBuilder.build_plane_parameters(
                p1=plane_data.points[0],
                p2=plane_data.points[1],
                p3=plane_data.points[2])

            filtered_coordinates = CoordinatesFilter.filter_coordinates_related_to_plane(
                filtered_coordinates,
                A, B, C, D,
                min_distance=distance_from_plane,
                direction=plane_data.direction)

        return filtered_coordinates

    @staticmethod
    def _filter_atoms_relates_carbone_atoms(
        coordinates_al: ndarray,
        coordinates_carbone: ndarray,
        max_distance_to_carbone_atoms: float,
    ) -> ndarray:
        """
        Filter points from the coordinates_al array
        by the max distance (max_distance_to_carbone_atoms param)
        to the points in the coordinates_carbone array.
        """

        # Calculate the distance matrix
        distances_matrix: ndarray = cdist(coordinates_al, coordinates_carbone)

        # Find the minimum distance for each atom in coordinates_al to any atom in coordinates_carbone
        min_distances: ndarray = np.min(distances_matrix, axis=1)

        # Filter the atoms in coordinates_al based on the max distance to carbone atoms
        filtered_coordinates: ndarray = coordinates_al[min_distances >= max_distance_to_carbone_atoms]

        return filtered_coordinates
