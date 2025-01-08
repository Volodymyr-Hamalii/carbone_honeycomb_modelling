import math
from pathlib import Path

import numpy as np
# from scipy.spatial.distance import cdist
# from scipy.optimize import minimize

from src.projects.carbon_honeycomb_actions.channel.planes.plane_polygons.carbon_honeycomb_plane_polygons import CarbonHoneycombHexagon
from src.utils import Constants, PathBuilder, Logger, execution_time_logger
from src.structure_visualizer import AtomsUniverseBuilder
from src.data_preparation import StructureSettings
from src.base_structure_classes import Points, FlatFigure, AlLatticeType, CoordinateLimits
from src.coordinate_operations import PointsRotator
from src.projects.carbon_honeycomb_actions import CarbonHoneycombChannel, CarbonHoneycombPlane

from ...structure_operations import StructureTranslator


logger = Logger("IntercalatedChannelBuilder")


class IntercalatedChannelBuilderBasedOnPlaneConfigs:
    @classmethod
    def build_al_in_carbon(
            cls,
            carbon_channel: CarbonHoneycombChannel,
            structure_settings: StructureSettings,
            equidistant_al_points: bool = True,
    ) -> Points:
        """ Returns coordinates_al """

        coordinates_al: Points = cls.calculate_coordinates_al(
            carbon_channel=carbon_channel,
            structure_settings=structure_settings)

        # if equidistant_al_points:
        #     # Set Al atoms maximally equidistant from the channel atoms
        #     coordinates_al: Points = AlAtomsSetter.equidistant_points_sets_in_channel(
        #         carbon_channel, coordinates_al, structure_settings)

        return coordinates_al

    @classmethod
    def calculate_coordinates_al(
        cls,
        carbon_channel: CarbonHoneycombChannel,
        structure_settings: StructureSettings,
    ) -> Points:

        coordinates_al: Points = cls._build_al_atoms_near_planes(carbon_channel)

        return coordinates_al

    @classmethod
    def _build_al_atoms_near_planes(cls, carbon_channel: CarbonHoneycombChannel) -> Points:
        coordinates_al: list[np.ndarray] = []
        carbon_channel_center: np.ndarray = carbon_channel.channel_center

        for plane in carbon_channel.planes:
            plane_coordinates_al: list[np.ndarray] = cls._build_al_atoms_near_polygons(
                polygons=plane.hexagons,  # type: ignore
                carbon_channel_center=carbon_channel_center,
            )

            coordinates_al.extend(plane_coordinates_al)

            plane_coordinates_al: list[np.ndarray] = cls._build_al_atoms_near_polygons(
                polygons=plane.pentagons,  # type: ignore
                carbon_channel_center=carbon_channel_center,
            )

            coordinates_al.extend(plane_coordinates_al)

        return Points(points=np.array(coordinates_al))

    @classmethod
    def _build_al_atoms_near_polygons(
            cls,
            polygons: list[FlatFigure],
            carbon_channel_center: np.ndarray,
    ) -> list[np.ndarray]:
        plane_coordinates_al: list[np.ndarray] = []
        for polygon in polygons:
            al_atom: np.ndarray = cls._build_al_atom_near_polygon(polygon, carbon_channel_center)
            plane_coordinates_al.append(al_atom)
        return plane_coordinates_al

    @staticmethod
    def _build_al_atom_near_polygon(
            polygon: FlatFigure,
            carbon_channel_center: np.ndarray,
            max_iter: int = 50,
            tol: float = 1e-5,
    ) -> np.ndarray:
        distance_from_carbon_atoms: float = 2.

        polygon_center: np.ndarray = polygon.center
        point_for_line: np.ndarray = np.array([
            carbon_channel_center[0],
            carbon_channel_center[1],
            polygon_center[2],
        ], dtype=float)

        # If polygon_center and point_for_line are identical, avoid degenerate case
        if np.allclose(polygon_center, point_for_line):
            # fallback: just return polygon_center or raise an error
            return polygon_center

        # 2) Parameterize the line: P(t) = start + t*(end - start)
        start = polygon_center.astype(float)
        end = point_for_line.astype(float)

        # 3) Define the function for average distance
        all_points = polygon.points.astype(float)  # shape (N,3)

        def average_distance(t: float) -> float:
            pt = start + t * (end - start)  # current candidate
            dists = np.linalg.norm(all_points - pt, axis=1)
            return dists.mean()

        # 4) Binary search for t in [0, large_upper_bound]
        left, right = 0.0, 5.0  # pick an upper bound you expect is "far enough"

        # If even at t=0 or t=right, we can't get close to desired distance, no guaranteed solution
        # but we'll proceed with best effort or expand the search if needed
        for _ in range(max_iter):
            mid: float = 0.5 * (left + right)
            mid_val: float = average_distance(mid)

            # compare to the desired distance
            if abs(mid_val - distance_from_carbon_atoms) < tol:
                # Good enough
                break

            if mid_val < distance_from_carbon_atoms:
                # We need to increase distance => move further along the line
                left: float = mid
            else:
                # We need to decrease distance => move closer
                right: float = mid

        # mid is our approximate solution
        t_best: float = 0.5 * (left + right)
        al_atom = start + t_best * (end - start)
        return al_atom
