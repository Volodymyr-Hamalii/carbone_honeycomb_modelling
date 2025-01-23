import numpy as np

from src.utils import Logger, execution_time_logger, Constants
from src.base_structure_classes import Points
from src.coordinate_operations import DistanceMeasure


logger = Logger("AtomsFilter")


class AtomsFilter:
    @classmethod
    # @execution_time_logger
    def replace_nearby_atoms_with_one_atom(cls, coordinates_al: Points) -> Points:
        """
        If there are 2 or more atoms have the distance between them less than dist_between_al / 3
        replace these atoms with one that places on the center between the replaced group.

        Returnes the updated coordinates_al (with replaced atoms).
        """

        dist_matrix: np.ndarray = DistanceMeasure.calculate_dist_matrix(coordinates_al.points)

        min_allowed_dist: float = Constants.phys.al.DIST_BETWEEN_ATOMS / 3

        # Keep track of groups of atoms to merge
        merged_indices = set()
        points_upd: list[np.ndarray] = []

        # Iterate through each atom
        for i, point in enumerate(coordinates_al.points):
            if i in merged_indices:
                continue  # Skip atoms that are already merged

            # Find atoms within the min_allowed_dist
            close_atoms: np.ndarray = np.where(dist_matrix[i] < min_allowed_dist)[0]

            # Merge these atoms if some close atoms are found
            if len(close_atoms) > 0:
                close_atoms = np.concatenate((close_atoms, [i]))
                merged_indices.update(close_atoms)

                # Calculate the center of the group
                group_center: np.ndarray = coordinates_al.points[close_atoms].mean(axis=0)
                points_upd.append(group_center)
            else:
                # Keep this atom as is if no nearby atoms are found
                points_upd.append(point)

        coordinates_al_upd = Points(points=np.array(points_upd))
        dist_matrix_upd: np.ndarray = DistanceMeasure.calculate_dist_matrix(coordinates_al_upd.points)

        if np.any(dist_matrix_upd < min_allowed_dist):
            # Run replacement one more time
            return cls.replace_nearby_atoms_with_one_atom(coordinates_al_upd)

        return coordinates_al_upd

    @staticmethod
    @execution_time_logger
    def remove_too_close_atoms(
            coordinates_al: Points,
            min_allowed_dist: float = Constants.phys.al.MIN_ALLOWED_DIST_BETWEEN_ATOMS,
    ) -> Points:
        """
        If there are 2 or more atoms have the distance between them less allowed
        keep only one atom and remove the rest.

        Returnes the updated coordinates_al (with removed atoms).
        """

        dist_matrix: np.ndarray = DistanceMeasure.calculate_dist_matrix(coordinates_al.points)
        # min_dists: np.ndarray = np.min(dist_matrix, axis=1)

        # Keep track of groups of atoms to merge
        merged_indices = set()
        points_upd: list[np.ndarray] = []

        # Iterate through each atom
        for i, point in enumerate(coordinates_al.points):
            if i in merged_indices:
                continue  # Skip atoms that are already merged

            points_upd.append(point)

            # Find atoms within the min_allowed_dist
            close_atoms: np.ndarray = np.where(dist_matrix[i] < min_allowed_dist)[0]

            # Skip too close atoms some close atoms are found
            if len(close_atoms) > 0:
                close_atoms = np.concatenate((close_atoms, [i]))
                merged_indices.update(close_atoms)

        return Points(points=np.array(points_upd))

    @classmethod
    # @execution_time_logger
    def remove_some_close_atoms(
            cls,
            coordinates_al: Points,
            min_dist: float,
            max_neighbours: int = -1,
            num_of_points_to_skip: int = 0,

            # inner params
            percent_to_remove: float = 1,
            init_amount: int = 0,
            removed_amount: int = 0,
    ) -> tuple[Points, int]:
        """
        Works recursively:
        1. Sort atoms by the biggest num of close neighbours and the min ave dist (if there are no such atoms - return the initial coordinates_al)
        2. Take the indexes of the atoms with the biggest num of close neighbours
        3. Iterate through each atom index
            3.1. Delete atom with provided index (=0 by default)
            3.2. Run the recursive function
            3.3. If the returned number of atoms bigger than the previous one - save the result
        4. Return the best result (with the biggest number of atoms)

        Returnes the updated coordinates_al (with removed atoms).
        """

        if percent_to_remove < 1 and init_amount == 0:
            init_amount = len(coordinates_al)

        dist_matrix: np.ndarray = DistanceMeasure.calculate_dist_matrix(coordinates_al.points)

        # Create a mask for values less than `min_dist`
        mask: np.ndarray = dist_matrix < min_dist

        # Replace values not satisfying the condition with `NaN`
        filtered_matrix: np.ndarray = np.where(mask, dist_matrix, np.nan)

        num_of_too_close_neighbours: np.ndarray = np.sum(~np.isnan(filtered_matrix), axis=1)

        if np.all(num_of_too_close_neighbours == 0):
            return coordinates_al, max_neighbours

        # Get the first max_neighbours
        if max_neighbours == -1:
            max_neighbours = np.max(num_of_too_close_neighbours)

        # Sort atoms by the biggest num of close neighbours and the min ave dist
        max_neighbour_indices = np.where(num_of_too_close_neighbours >= max_neighbours)[0]

        if len(max_neighbour_indices) == 0:
            return coordinates_al, max_neighbours

        if num_of_points_to_skip >= len(max_neighbour_indices):
            num_of_points_to_skip = len(max_neighbour_indices) - 1

        # try:
        new_points: np.ndarray = np.delete(coordinates_al.points, max_neighbour_indices[num_of_points_to_skip], axis=0)
        # except IndexError:
        #     pass

        new_coordinates_al = Points(points=new_points)

        if percent_to_remove < 1:
            removed_amount += 1
            if removed_amount > init_amount * percent_to_remove:
                return new_coordinates_al, max_neighbours

        result, max_neighbours = cls.remove_some_close_atoms(
            new_coordinates_al, min_dist, max_neighbours,

            percent_to_remove=percent_to_remove,
            init_amount=init_amount,
            removed_amount=removed_amount,
        )

        return result, max_neighbours

        # best_result: Points = coordinates_al
        # max_atoms_count: int = 0

        # for idx in max_neighbour_indices:
        #     new_points: np.ndarray = np.delete(coordinates_al.points, idx, axis=0)
        #     new_coordinates_al = Points(points=new_points)
        #     result: Points = cls.remove_some_close_atoms(new_coordinates_al, min_dist, max_neighbours)

        #     if len(result.points) > max_atoms_count:
        #         best_result = result
        #         max_atoms_count = len(result.points)

        # return best_result
