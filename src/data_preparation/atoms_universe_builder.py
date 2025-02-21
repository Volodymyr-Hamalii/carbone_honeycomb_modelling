from pathlib import Path

import MDAnalysis as mda
import numpy as np
from numpy import ndarray

from src.utils import Logger
from src.base_structure_classes import Points, CoordinateLimits


import warnings
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')


logger = Logger("AtomsUniverseBuilder")


class AtomsUniverseBuilder:
    @staticmethod
    def builds_atoms_coordinates(
        path_to_pdb_file: Path,
    ) -> Points:
        """ 
        Builds atom's Universe based on the {path_to_pdb_file}.pdb file.
        If coordinate_limits is not None -- also filters atoms by provided limits.
        """

        # Load a structure from a file
        u = mda.Universe(path_to_pdb_file)

        # Access atoms and coordinates
        atoms = u.atoms
        points: ndarray = atoms.positions  # type: ignore

        return Points(points=points)

    @classmethod
    def build_hcp_lattice_type(
        cls,
        dist_between_atoms: float,
        coordinate_limits: CoordinateLimits,
    ) -> Points:
        """ Generate coordinates for planes with a 'ABAB' close-packed stacking sequence (Hexagonal close-packed) """

        # Extract limits
        # x_min, x_max, y_min, y_max, z_min, z_max = cls._get_extended_limits(
        #     coordinate_limits)
        x_min: float = coordinate_limits.x_min
        x_max: float = coordinate_limits.x_max
        y_min: float = coordinate_limits.y_min
        y_max: float = coordinate_limits.y_max
        z_min: float = coordinate_limits.z_min
        z_max: float = coordinate_limits.z_max

        a: float = dist_between_atoms

        # Regular tetrahedron altitude
        reg_tet_alt: float = np.sqrt(6) * a / 3

        # Equilateral triangle altitude
        eq_tr_alt: float = np.sqrt(3) * a / 2

        # Define the base HCP layer (A and B layers are different)
        base_layer_a: ndarray = np.array([
            [0, 0, 0],
            [a, 0, 0],
            [a / 2, eq_tr_alt, 0],
            [3 * a / 2, eq_tr_alt, 0],
        ])

        base_layer_b: ndarray = np.array([
            [a / 2, eq_tr_alt / 3, 0],
            [3 * a / 2, eq_tr_alt / 3, 0],
            [a, 4 * eq_tr_alt / 3, 0],
            [2 * a, 4 * eq_tr_alt / 3, 0],
        ])

        # Extend the unit cell to fill the volume
        x_step: float = 2 * a
        x_range: ndarray = np.arange(x_min - x_step, x_max + x_step, x_step)

        y_step: float = 2 * eq_tr_alt
        y_range: ndarray = np.arange(y_min - y_step, y_max + y_step, y_step)

        z_step: float = reg_tet_alt
        z_range: ndarray = np.arange(z_min - z_step, z_max + z_step, z_step)

        all_points = []
        for i, z in enumerate(z_range):
            base_layer = base_layer_a if i % 2 == 0 else base_layer_b
            for x in x_range:
                for y in y_range:
                    translated_layer = base_layer + np.array([x, y, z])
                    all_points.append(translated_layer)

        # Flatten the list of layers into a single array
        full_structure: ndarray = np.concatenate(all_points)

        # Filter the structure to be within the channel limits
        within_limits = (
            (full_structure[:, 0] >= x_min) & (full_structure[:, 0] <= x_max) &
            (full_structure[:, 1] >= y_min) & (full_structure[:, 1] <= y_max) &
            (full_structure[:, 2] >= z_min) & (full_structure[:, 2] <= z_max)
        )

        return Points(points=full_structure[within_limits])

    @classmethod
    def build_fcc_lattice_type(
        cls,
        dist_between_atoms: float,
        coordinate_limits: CoordinateLimits,
    ) -> Points:
        """
        Generate coordinates for planes with a 'ABCABC' close-packed stacking sequence (Face-centered cubic).
        """

        # Extract limits
        # x_min, x_max, y_min, y_max, z_min, z_max = cls._get_extended_limits(
        #     coordinate_limits)
        x_min: float = coordinate_limits.x_min
        x_max: float = coordinate_limits.x_max
        y_min: float = coordinate_limits.y_min
        y_max: float = coordinate_limits.y_max
        z_min: float = coordinate_limits.z_min
        z_max: float = coordinate_limits.z_max

        a: float = dist_between_atoms

        # Regular tetrahedron altitude
        reg_tet_alt: float = np.sqrt(6) * a / 3

        # Equilateral triangle altitude
        eq_tr_alt: float = np.sqrt(3) * a / 2

        # FCC base planes: A, B, and C layers are different
        base_layer_a: ndarray = np.array([
            [0, 0, 0],
            [a, 0, 0],
            [a / 2, eq_tr_alt, 0],
            [3 * a / 2, eq_tr_alt, 0],
        ])

        base_layer_b: ndarray = np.array([
            [a / 2, eq_tr_alt / 3, 0],
            [3 * a / 2, eq_tr_alt / 3, 0],
            [a, 4 * eq_tr_alt / 3, 0],
            [2 * a, 4 * eq_tr_alt / 3, 0],
        ])

        base_layer_c: ndarray = np.array([
            [0, 2 * eq_tr_alt / 3, 0],
            [a, 2 * eq_tr_alt / 3, 0],
            [a / 2, 5 * eq_tr_alt / 3, 0],
            [3 * a / 2, 5 * eq_tr_alt / 3, 0],
        ])

        # return np.concatenate([base_layer_a, base_layer_b, base_layer_c])

        layers: list[ndarray] = [base_layer_a, base_layer_b, base_layer_c]

        # Extend the unit cell to fill the volume
        x_step: float = 2 * a
        x_range: ndarray = np.arange(x_min - x_step, x_max + x_step, x_step)

        y_step: float = 2 * eq_tr_alt
        y_range: ndarray = np.arange(y_min - y_step, y_max + y_step, y_step)

        z_step: float = reg_tet_alt
        z_range: ndarray = np.arange(z_min - z_step, z_max + z_step, z_step)

        all_points: list[ndarray] = []
        for i, z in enumerate(z_range):
            base_layer: ndarray = layers[i % 3]  # Cycle through A, B, C layers
            for x in x_range:
                for y in y_range:
                    translated_layer = base_layer + np.array([x, y, z])
                    all_points.append(translated_layer)

        # Flatten the list of layers into a single array
        full_structure: ndarray = np.concatenate(all_points)

        # Filter the structure to be within the channel limits
        within_limits = (
            (full_structure[:, 0] >= x_min) & (full_structure[:, 0] <= x_max) &
            (full_structure[:, 1] >= y_min) & (full_structure[:, 1] <= y_max) &
            (full_structure[:, 2] >= z_min) & (full_structure[:, 2] <= z_max)
        )

        return Points(points=full_structure[within_limits])

    @classmethod
    def _get_extended_limits(
        cls, coordinate_limits: CoordinateLimits
    ) -> tuple[float, float, float, float, float, float]:

        x_min, x_max = cls._get_extended_limits_for_coordinate(
            min_coord=coordinate_limits.x_min,
            max_coord=coordinate_limits.x_max,
        )

        y_min, y_max = cls._get_extended_limits_for_coordinate(
            min_coord=coordinate_limits.y_min,
            max_coord=coordinate_limits.y_max,
        )

        z_min, z_max = cls._get_extended_limits_for_coordinate(
            min_coord=coordinate_limits.z_min,
            max_coord=coordinate_limits.z_max,
        )

        return x_min, x_max, y_min, y_max, z_min, z_max

    @staticmethod
    def _get_extended_limits_for_coordinate(min_coord: float, max_coord: float) -> tuple[float, float]:
        delta: float = abs((max_coord - min_coord)) / 2
        min_coord_extended: float = min_coord - delta
        max_coord_extended: float = max_coord + delta
        return min_coord_extended, max_coord_extended
