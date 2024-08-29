import os
import MDAnalysis as mda
import numpy as np
from numpy import ndarray

from ..utils import Logger
from ..data_preparation import ChannelLimits
from ..base_structure_classes import LatticeType

from .structure_utils import StructureUtils

import warnings
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')


logger = Logger(__name__)


class AtomsUniverseBuilder:
    @staticmethod
    def builds_atoms_coordinates(
        path_to_pdb_file: str,
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

            # Write a result
            one_channel_pdb_file_name: str = path_to_pdb_file.replace(".pdb", "-one-channel.pdb")

            # Create result-data/{structure_folder}/ljout-result-one-channel.pdb file if it didn't exist
            if not os.path.exists(one_channel_pdb_file_name):
                logger.info("Created", one_channel_pdb_file_name)
                StructureUtils.write_pdb_from_mda(one_channel_pdb_file_name, single_channel_atoms)

            return single_channel_atoms.positions

        return coordinates

    @staticmethod
    def build_hcp_lattice_type(
        lattice_parameter: float,
        channel_coordinate_limits: ChannelLimits,
        z_min_default: float = 0,
        z_max_default: float = 12,
    ) -> ndarray:
        """ Generate coordinates for planes with a 'ABAB' close-packed stacking sequence (Hexagonal close-packed) """

        # Extract limits
        x_min, x_max = channel_coordinate_limits.x_min, channel_coordinate_limits.x_max
        y_min, y_max = channel_coordinate_limits.y_min, channel_coordinate_limits.y_max
        z_min = channel_coordinate_limits.z_min or z_min_default
        z_max = channel_coordinate_limits.z_max or z_max_default

        c_to_a_ratio: float = np.sqrt(8 / 3)  # Ideal c/a ratio for HCP
        a: float = lattice_parameter
        c: float = c_to_a_ratio * a

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

        z_step: float = c / 2
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

    @staticmethod
    def build_fcc_lattice_type(
        lattice_parameter: float,
        channel_coordinate_limits: ChannelLimits,
    ) -> ndarray:
        """ Generate coordinates for planes with a 'ABCABC' close-packed stacking sequence (Face-centered cubic) """

        # Extract limits
        x_min, x_max = channel_coordinate_limits.x_min, channel_coordinate_limits.x_max
        y_min, y_max = channel_coordinate_limits.y_min, channel_coordinate_limits.y_max
        z_min, z_max = channel_coordinate_limits.z_min, channel_coordinate_limits.z_max

        a = lattice_parameter

        # Define the base FCC layers (A, B, and C layers)
        base_layer_a = np.array([
            [0, 0, 0],
            [a / 2, a / 2, 0],
            [a, a, 0]
        ])

        base_layer_b = base_layer_a + np.array([0.5 * a, 0.5 * a, 0.5 * a])
        base_layer_c = base_layer_a + np.array([a, a, a])

        layers = [base_layer_a, base_layer_b, base_layer_c]

        # Extend the unit cell to fill the volume
        x_range = np.arange(x_min, x_max, a)
        y_range = np.arange(y_min, y_max, a)
        z_range = np.arange(z_min, z_max, a / np.sqrt(2))

        all_points = []
        for i, z in enumerate(z_range):
            layer_index = i % 3  # Cycle through A, B, C layers
            base_layer = layers[layer_index]
            for x in x_range:
                for y in y_range:
                    translated_layer = base_layer + np.array([x, y, z])
                    all_points.append(translated_layer)

        # Flatten the list of layers into a single array
        full_structure = np.concatenate(all_points)

        # Filter the structure to be within the channel limits
        within_limits = (
            (full_structure[:, 0] >= x_min) & (full_structure[:, 0] <= x_max) &
            (full_structure[:, 1] >= y_min) & (full_structure[:, 1] <= y_max) &
            (full_structure[:, 2] >= z_min) & (full_structure[:, 2] <= z_max)
        )

        return full_structure[within_limits]
