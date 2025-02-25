import numpy as np
from scipy.optimize import minimize

from src.utils import Constants, Logger
from src.base_structure_classes import Points
from src.coordinate_operations import PointsMover, DistanceMeasure
from src.projects.carbon_honeycomb_actions import (
    CarbonHoneycombChannel,
)
from ..by_variance import (
    AlAtomsFilter,
)
from .al_atoms_translator import AlAtomsTranslator


logger = Logger(__name__)


class FullChannelBuilder:
    @classmethod
    def build_full_channel(
            cls,
            carbon_channel: CarbonHoneycombChannel,
            al_channel_planes_coordinates: Points,
            al_bulk_coordinates: Points
    ) -> Points:

        al_bulk_filtered_related_channel_planes: Points = cls._filter_related_channel_planes_al_atoms(
            carbon_channel=carbon_channel,
            al_bulk=al_bulk_coordinates,
        )

        al_bulk_optimized_positions: Points = cls._find_optimal_positions_for_al_atoms(
            al_bulk_coordinates=al_bulk_filtered_related_channel_planes,
            al_channel_planes_coordinates=al_channel_planes_coordinates,
        )

        return Points(
            points=np.vstack([al_channel_planes_coordinates.points, al_bulk_optimized_positions.points])
        )

    @staticmethod
    def _filter_related_channel_planes_al_atoms(
        carbon_channel: CarbonHoneycombChannel,
        al_bulk: Points,
    ) -> Points:

        al_bulk_coordinates: Points = AlAtomsFilter.filter_atoms_related_clannel_planes(
            coordinates_al=al_bulk,
            carbon_channel=carbon_channel,
            distance_from_plane=Constants.phys.al.MIN_ALLOWED_DIST_BETWEEN_ATOMS,
        )

        return al_bulk_coordinates

    @staticmethod
    def _filter_related_plane_al_atoms(
        al_bulk_coordinates: Points,
        al_channel_planes_coordinates: Points,
    ) -> Points:
        # Find the minimum distance for each atom in coordinates_al to any atom in coordinates_carbon
        min_distances: np.ndarray = DistanceMeasure.calculate_min_distances(
            al_bulk_coordinates.points, al_channel_planes_coordinates.points
        )

        filtered_al_coordinates: Points = Points(
            points=al_bulk_coordinates.points[min_distances >= Constants.phys.al.MIN_RECOMENDED_DIST_BETWEEN_ATOMS]
        )

        return filtered_al_coordinates

    @classmethod
    def _find_optimal_positions_for_al_atoms(
        cls,
        al_bulk_coordinates: Points,
        al_channel_planes_coordinates: Points,
    ) -> Points:
        init_vector: np.ndarray = np.array([0.0, 0.0])

        result = minimize(
            cls.objective_function,
            init_vector,
            args=(al_bulk_coordinates, al_channel_planes_coordinates),
            method="BFGS",
            options={"disp": True}
        )

        vector_to_move: np.ndarray = result.x

        if len(vector_to_move) == 2:
            vector_to_move = np.append(vector_to_move, 0.0)

        moved_al_points: Points = PointsMover.move_on_vector(
            points=al_bulk_coordinates,
            vector=vector_to_move,
        )

        return moved_al_points

    @classmethod
    def objective_function(
        cls,
        vector_to_move: np.ndarray,
        al_bulk_coordinates: Points,
        al_channel_planes_coordinates: Points,
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
