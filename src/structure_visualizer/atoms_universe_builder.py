import os
import MDAnalysis as mda
import numpy as np
from numpy import ndarray

from ..utils import Logger
from ..data_preparation import ChannelLimits
from ..intercalation import AlLatticeType

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
    def build_close_packed_structure(
            lattice_parameter: float,
            lattice_type: AlLatticeType,
            channel_coordinate_limits: ChannelLimits,
            to_translate_al: bool = True,  # TODO
    ) -> ndarray:
        """
        Build coordinates for HCP (hexagonal close-packed) or FCC (face-centered cubic) close-packed lattices
        with lattice_parameter distance between atoms
        """

        x_min: float = channel_coordinate_limits.x_min
        x_max: float = channel_coordinate_limits.x_max
        y_min: float = channel_coordinate_limits.y_min
        y_max: float = channel_coordinate_limits.y_max
        z_min: float = channel_coordinate_limits.z_min
        z_max: float = channel_coordinate_limits.z_max

        # For HCP:
        if lattice_type.is_hcp:
            c_to_a_ratio: float = np.sqrt(8 / 3)  # Ideal c/a ratio for HCP
            a: float = lattice_parameter
            c: float = c_to_a_ratio * a

            # Base HCP coordinates within a unit cell
            base_points: ndarray = np.array([
                [0, 0, 0],
                [2 / 3 * a, 1 / 3 * a * np.sqrt(3), 0],
                [1 / 3 * a, 2 / 3 * a * np.sqrt(3), 0],
                [0, 0, c / 2],
                [2 / 3 * a, 1 / 3 * a * np.sqrt(3), c / 2],
                [1 / 3 * a, 2 / 3 * a * np.sqrt(3), c / 2],
            ])

            # Extend the unit cell to fill the volume
            x_range: ndarray = np.arange(x_min, x_max, a)
            y_range: ndarray = np.arange(y_min, y_max, a * np.sqrt(3))
            z_range: ndarray = np.arange(z_min, z_max, c)

        # For FCC:
        elif lattice_type.is_fcc:
            a: float = lattice_parameter
            # Base FCC coordinates within a unit cell
            base_points = np.array([
                [0, 0, 0],
                [0, a / 2, a / 2],
                [a / 2, 0, a / 2],
                [a / 2, a / 2, 0],
            ])

            # Extend the unit cell to fill the volume
            x_range: ndarray = np.arange(x_min, x_max, a)
            y_range: ndarray = np.arange(y_min, y_max, a)
            z_range: ndarray = np.arange(z_min, z_max, a)
        else:
            raise ValueError("Set correct packing: 'HCP' or 'FCC'")

        # Generate the full structure by replicating the base points in the x, y, and z ranges
        all_points = []
        for x in x_range:
            for y in y_range:
                for z in z_range:
                    new_points = base_points + np.array([x, y, z])
                    all_points.append(new_points)

        # Concatenate all generated points
        full_structure: ndarray = np.concatenate(all_points)

        # Filter the structure to be within the channel limits
        within_limits = (
            (full_structure[:, 0] >= x_min) & (full_structure[:, 0] <= x_max) &
            (full_structure[:, 1] >= y_min) & (full_structure[:, 1] <= y_max) &
            (full_structure[:, 2] >= z_min) & (full_structure[:, 2] <= z_max)
        )
        return full_structure[within_limits]
