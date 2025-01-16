import numpy as np

from src.utils import Logger, Constants
from src.base_structure_classes import Points
from src.coordinate_operations import DistanceMeasure


logger = Logger("AtomsFilter")


class AtomsFilter:
    @staticmethod
    def replace_nearby_atoms_with_one_atom(coordinates_al: Points, k: float = 0.3) -> Points:
        """
        If there are 2 or more atoms have the distance between them less than dist_between_al / 3
        replace these atoms with one that places on the center between the replaced group.

        Returnes the updated coordinates_al (with replaced atoms).
        """

        dist_matrix: np.ndarray = DistanceMeasure.calculate_dist_matrix(coordinates_al.points)
        # min_dists: np.ndarray = np.min(dist_matrix, axis=1)

        dist_between_al: float = Constants.phys.AL_DIST_BETWEEN_ATOMS
        min_allowed_dist: float = k * dist_between_al

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
            return cls.replace_nearby_atoms_with_one_atom(coordinates_al_upd, k)

        return coordinates_al_upd

    @staticmethod
    def remove_too_close_atoms(coordinates_al: Points, k: float = 0.7) -> Points:
        """
        If there are 2 or more atoms have the distance between them more than k*dist_between_al
        keep only one atom and remove the rest.

        Returnes the updated coordinates_al (with replaced atoms).
        """

        dist_matrix: np.ndarray = DistanceMeasure.calculate_dist_matrix(coordinates_al.points)
        # min_dists: np.ndarray = np.min(dist_matrix, axis=1)

        dist_between_al: float = Constants.phys.AL_DIST_BETWEEN_ATOMS
        min_allowed_dist: float = k * dist_between_al

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
