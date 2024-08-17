import numpy as np
from numpy import ndarray

from ..utils import Logger
from ..data_preparation import ChannelLimits
from .coordinates_filter import CoordinatesFilter


logger = Logger(__name__)


class StructureTranslator:
    @classmethod
    def translate_cell(
            cls,
            cell_coordinates: ndarray,
            translation_limits: ChannelLimits | None,
            translation_step_x: float = 0,
            translation_step_y: float = 0,
            translation_step_z: float = 0,
    ) -> ndarray:
        """
        Translate elemental cell along x, y and z axis.
        If translation_steps are not provided -- they will be calculated
        as a distances along specific axis in the elemental cell.
        """

        if translation_limits is None:
            logger.error("Translation cannot be done without translation_limits.")
            return cell_coordinates

        # Determine the translation steps based on the provided coordinates
        translation_step_x = cls._find_translation_step("x", translation_step_x, cell_coordinates[:, 0])
        translation_step_y = cls._find_translation_step("y", translation_step_y, cell_coordinates[:, 1])
        translation_step_z = cls._find_translation_step("z", translation_step_z, cell_coordinates[:, 2])

        translated_coordinates_list: list[ndarray] = []

        x_min: float = translation_limits.x_min
        x_max: float = translation_limits.x_max
        y_min: float = translation_limits.y_min
        y_max: float = translation_limits.y_max
        z_min: float = translation_limits.z_min
        z_max: float = translation_limits.z_max

        # Translate the structure in x, y, and z directions
        for x_translation in np.arange(x_min, x_max, translation_step_x):
            for y_translation in np.arange(y_min, y_max, translation_step_y):
                for z_translation in np.arange(z_min, z_max, translation_step_z):
                    translation_vector: ndarray = np.array([x_translation, y_translation, z_translation])
                    translated_structure: ndarray = cell_coordinates + translation_vector
                    translated_coordinates_list.append(translated_structure)

        # Concatenate all translated coordinates
        translated_coordinates: ndarray = np.vstack(translated_coordinates_list)

        # Remove duplicate coordinates
        translated_coordinates: ndarray = np.unique(translated_coordinates, axis=0)

        return CoordinatesFilter.filter_by_max_min_z(translated_coordinates, z_min, z_max)

    @staticmethod
    def _find_translation_step(axis: str, default_step: float, axis_coordinates: ndarray) -> float:
        translation_step: float = default_step or axis_coordinates.ptp()

        if translation_step == 0:
            raise ValueError(f"Translation step for {axis} axis = 0.")
        
        return translation_step
