import numpy as np
from scipy.optimize import minimize
from scipy.spatial.distance import cdist

from src.utils import Constants, ConstantsAtomParams, Logger
from src.base_structure_classes import Points
from src.coordinate_operations import PointsMover, DistanceMeasure
from src.projects.carbon_honeycomb_actions import (
    CarbonHoneycombChannel,
)
from ..by_variance import (
    AtomsFilter,
)
from .al_atoms_translator import AlAtomsTranslator


logger = Logger("FullChannelBuilder")


class FullChannelBuilder:
    @classmethod
    def build_full_channel(
            cls,
            carbon_channel: CarbonHoneycombChannel,
            al_channel_planes_coordinates: Points,
            al_bulk_coordinates: Points,
            atom_params: ConstantsAtomParams,
    ) -> Points:

        al_bulk_filtered_related_channel_planes: Points = cls._filter_related_channel_planes_al_atoms(
            carbon_channel=carbon_channel,
            al_bulk=al_bulk_coordinates,
            atom_params=atom_params,
        )

        al_bulk_optimized_positions: Points = cls._find_optimal_positions_for_al_atoms(
            al_bulk_coordinates=al_bulk_filtered_related_channel_planes,
            al_channel_planes_coordinates=al_channel_planes_coordinates,
            atom_params=atom_params,
        )

        al_bulk_adjusted: Points = cls._adjust_the_closest_al_atoms(
            al_bulk_coordinates=al_bulk_optimized_positions,
            al_channel_planes_coordinates=al_channel_planes_coordinates,
            atom_params=atom_params,
        )

        logger.info(f"Al bulk atoms added: {len(al_bulk_adjusted.points)}")

        return Points(
            points=np.vstack([al_channel_planes_coordinates.points, al_bulk_adjusted.points])
        )

    @staticmethod
    def _filter_related_channel_planes_al_atoms(
        carbon_channel: CarbonHoneycombChannel,
        al_bulk: Points,
        atom_params: ConstantsAtomParams,
    ) -> Points:

        al_bulk_coordinates: Points = AtomsFilter.filter_atoms_related_clannel_planes(
            inter_points=al_bulk,
            carbon_channel=carbon_channel,
            distance_from_plane=atom_params.MIN_ALLOWED_DIST_BETWEEN_ATOMS,
        )

        return al_bulk_coordinates

    @staticmethod
    def _filter_related_plane_al_atoms(
        al_bulk_coordinates: Points,
        al_channel_planes_coordinates: Points,
        atom_params: ConstantsAtomParams,
    ) -> Points:
        # Find the minimum distance for each atom in coordinates_al to any atom in coordinates_carbon
        min_distances: np.ndarray = DistanceMeasure.calculate_min_distances(
            al_bulk_coordinates.points, al_channel_planes_coordinates.points
        )

        filtered_al_coordinates: Points = Points(
            points=al_bulk_coordinates.points[min_distances >= atom_params.MIN_RECOMENDED_DIST_BETWEEN_ATOMS]
        )

        return filtered_al_coordinates

    @classmethod
    def _find_optimal_positions_for_al_atoms(
        cls,
        al_bulk_coordinates: Points,
        al_channel_planes_coordinates: Points,
        atom_params: ConstantsAtomParams,
    ) -> Points:
        init_vector: np.ndarray = np.array([0.0, 0.0])

        result = minimize(
            cls._objective_function,
            init_vector,
            args=(al_bulk_coordinates, al_channel_planes_coordinates, atom_params),
            method="BFGS",
            options={"disp": True}
        )

        vector_to_move: np.ndarray = result.x

        if len(vector_to_move) == 2:
            vector_to_move = np.append(vector_to_move, 0.0)

        moved_al_bulk_coordinates: Points = PointsMover.move_on_vector(
            points=al_bulk_coordinates,
            vector=vector_to_move,
        )

        filtered_al_bulk_coordinates: Points = cls._filter_related_plane_al_atoms(
            al_bulk_coordinates=moved_al_bulk_coordinates,
            al_channel_planes_coordinates=al_channel_planes_coordinates,
            atom_params=atom_params,
        )

        return filtered_al_bulk_coordinates

    @classmethod
    def _objective_function(
        cls,
        vector_to_move: np.ndarray,
        al_bulk_coordinates: Points,
        al_channel_planes_coordinates: Points,
        atom_params: ConstantsAtomParams,
    ) -> float | np.floating:

        if len(vector_to_move) == 2:
            vector_to_move = np.append(vector_to_move, 0.0)

        moved_al_points: Points = PointsMover.move_on_vector(
            points=al_bulk_coordinates,
            vector=vector_to_move,
        )

        filtered_al_points: Points = cls._filter_related_plane_al_atoms(
            al_bulk_coordinates=moved_al_points,
            al_channel_planes_coordinates=al_channel_planes_coordinates,
            atom_params=atom_params,
        )

        if len(filtered_al_points) == 0:
            return np.inf

        # Calculate the maximum possible number of atoms (before filtering)
        max_possible_atoms: int = len(al_bulk_coordinates)

        # Calculate atom count penalty (increases as we lose more atoms)
        atoms_penalty: float = ((max_possible_atoms - len(filtered_al_points)) / max_possible_atoms) ** 2

        # Calculate distance variance as before
        min_dists_var: np.floating = cls._calculate_min_dists_var(
            filtered_al_points, al_channel_planes_coordinates
        )

        # Combine penalties with appropriate weights
        # The atom count penalty is dominant (multiplied by 1000)
        # The distance variance is secondary (multiplied by 1)
        return 1000 * atoms_penalty + min_dists_var

    @staticmethod
    def _calculate_min_dists_var(
        points_1: Points,
        points_2: Points,
    ) -> np.floating:
        min_distances: np.ndarray = DistanceMeasure.calculate_min_distances(
            points_1.points, points_2.points)

        return np.var(min_distances)

    @staticmethod
    def _adjust_the_closest_al_atoms(
        al_bulk_coordinates: Points,
        al_channel_planes_coordinates: Points,
        atom_params: ConstantsAtomParams,
    ) -> Points:
        """
        Check distances between atoms in al_bulk_coordinates and al_channel_planes_coordinates.
        If the distance is less than Constants.phys.al.DIST_BETWEEN_ATOMS,
        get 3 the closest atoms in al_bulk_coordinates and all atoms in al_channel_planes_coordinates
        that are closer than Constants.phys.al.DIST_BETWEEN_ATOMS and move the atom to the center of these atoms.
        """
        min_dists: np.ndarray = DistanceMeasure.calculate_min_distances(
            al_bulk_coordinates.points, al_channel_planes_coordinates.points
        )

        counter: int = 0

        for i in range(len(min_dists)):
            if min_dists[i] < atom_params.DIST_BETWEEN_ATOMS:
                dists_plane_bulk: np.ndarray = cdist(al_bulk_coordinates.points, al_channel_planes_coordinates.points)
                dists_bulk: np.ndarray = DistanceMeasure.calculate_dist_matrix(al_bulk_coordinates.points)

                # Get the closest atoms from the plane
                closest_atoms_from_plane: np.ndarray = al_channel_planes_coordinates.points[
                    dists_plane_bulk[i] < atom_params.DIST_BETWEEN_ATOMS * 1.01]

                # Get the closest atoms from the bulk
                closest_atoms_from_bulk: np.ndarray = al_bulk_coordinates.points[
                    dists_bulk[i] < atom_params.DIST_BETWEEN_ATOMS * 1.01]

                if len(closest_atoms_from_plane) == 0 or len(closest_atoms_from_bulk) == 0:
                    continue

                adjusted_atom: np.ndarray = np.mean(
                    np.vstack([closest_atoms_from_plane, closest_atoms_from_bulk]), axis=0)

                al_bulk_coordinates.points[i] = adjusted_atom
                counter += 1

        logger.info(f"Adjusted {counter} atoms.")

        return al_bulk_coordinates
