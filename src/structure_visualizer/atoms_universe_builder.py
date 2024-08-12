import MDAnalysis as mda
import numpy as np
from numpy import ndarray

from ..data_preparation import ChannelLimits
from .structure_utils import StructureUtils

import warnings
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')


class AtomsUniverseBuilder:
    @staticmethod
    def builds_atoms_coordinates(
        path_to_pdb_file: str,
        channel_coordinate_limits: ChannelLimits | None = None,
    ) -> ndarray:

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
            StructureUtils.write_pdb_from_mda(one_channel_pdb_file_name, single_channel_atoms)

            return single_channel_atoms.positions

        return coordinates
