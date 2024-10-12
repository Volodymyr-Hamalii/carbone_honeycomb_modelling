from pathlib import Path

import MDAnalysis as mda
import numpy as np
from numpy import ndarray

from ..utils import Logger, Constants
from ..data_preparation import ChannelLimits
from ..base_structure_classes import LatticeType

from .structure_utils import StructureUtils

import warnings
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')


logger = Logger("AtomsUniverseBuilder")


class AtomsUniverseBuilder:
    @staticmethod
    def builds_atoms_coordinates(
        path_to_pdb_file: Path,
        channel_coordinate_limits: ChannelLimits | None = None,
    ) -> ndarray:
        """ 
        Builds atom's Universe based on the {path_to_pdb_file}.pdb file.
        If channel_coordinate_limits is not None -- also filters atoms by provided limits.
        """

        # Load a structure from a file
        u = mda.Universe(path_to_pdb_file)

        # Access atoms and coordinates
        atoms = u.atoms
        coordinates: ndarray = atoms.positions  # type: ignore

        if channel_coordinate_limits is not None:
            x_min: float = channel_coordinate_limits.x_min
            x_max: float = channel_coordinate_limits.x_max
            y_min: float = channel_coordinate_limits.y_min
            y_max: float = channel_coordinate_limits.y_max

            # Filter coordinates based on the criteria
            filtered_indices = np.where(
                (coordinates[:, 0] >= x_min) & (coordinates[:, 0] <= x_max) &
                (coordinates[:, 1] >= y_min) & (coordinates[:, 1] <= y_max)
            )[0]

            # Select atoms based on the filtered indices
            single_channel_atoms = atoms[filtered_indices]  # type: ignore

            # Build a path to result file
            one_channel_pdb_file_path: Path = Path(
                str(path_to_pdb_file).replace(
                    Constants.filenames.INIT_PDB_FILE, Constants.filenames.PDB_FILE_ONE_CHANNEL))

            # Create result-data/{structure_folder}/ljout-result-one-channel.pdb file if it didn't exist
            if not one_channel_pdb_file_path.exists():
                logger.info("Created", one_channel_pdb_file_path)
                StructureUtils.write_pdb_from_mda(one_channel_pdb_file_path, single_channel_atoms)

            return single_channel_atoms.positions

        return coordinates

    @classmethod
    def build_hcp_lattice_type(
        cls,
        lattice_parameter: float,
        channel_coordinate_limits: ChannelLimits
    ) -> ndarray:
        """ Generate coordinates for planes with a 'ABAB' close-packed stacking sequence (Hexagonal close-packed) """

        # Extract limits
        x_min, x_max, y_min, y_max, z_min, z_max = cls._get_extended_limits(
            channel_coordinate_limits)

        a: float = lattice_parameter

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

        return full_structure[within_limits]

    @classmethod
    def build_fcc_lattice_type(
        cls,
        lattice_parameter: float,
        channel_coordinate_limits: ChannelLimits,
    ) -> ndarray:
        """
        Generate coordinates for planes with a 'ABCABC' close-packed stacking sequence (Face-centered cubic).
        """

        # Extract limits
        x_min, x_max, y_min, y_max, z_min, z_max = cls._get_extended_limits(
            channel_coordinate_limits)

        a: float = lattice_parameter

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

        return full_structure[within_limits]

    @staticmethod
    def _get_extended_limits(
        channel_coordinate_limits: ChannelLimits
    ) -> tuple[float, float, float, float, float, float]:

        x_extend: float = channel_coordinate_limits.x_max / 2
        x_min: float = channel_coordinate_limits.x_min - x_extend
        x_max: float = channel_coordinate_limits.x_max + x_extend

        y_extend: float = channel_coordinate_limits.y_max / 2
        y_min: float = channel_coordinate_limits.y_min - y_extend
        y_max: float = channel_coordinate_limits.y_max + y_extend

        # z_extend: float = channel_coordinate_limits.z_max or z_max_default / 2
        z_extend: float = 0
        z_min: float = channel_coordinate_limits.z_min
        z_max: float = channel_coordinate_limits.z_max

        return x_min, x_max, y_min, y_max, z_min, z_max
