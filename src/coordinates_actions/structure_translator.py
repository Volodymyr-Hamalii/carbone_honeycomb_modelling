import numpy as np
from numpy import ndarray

from .coordinate_limits import CoordinateLimits


class StructureTranslator:
    @staticmethod
    def translate(
            coordinates: ndarray,
            translation_limits: CoordinateLimits,
    ) -> ndarray:
        x_min: float = translation_limits.x_min
        x_max: float = translation_limits.x_max
        y_min: float = translation_limits.y_min
        y_max: float = translation_limits.y_max
        z_min: float = translation_limits.z_min
        z_max: float = translation_limits.z_max

        # Determine the translation steps based on the provided coordinates
        translation_step_x: float = coordinates[:, 0].ptp() if coordinates[:, 0].ptp() != 0 else 1.0
        translation_step_y: float = coordinates[:, 1].ptp() if coordinates[:, 1].ptp() != 0 else 1.0
        translation_step_z: float = coordinates[:, 2].ptp() if coordinates[:, 2].ptp() != 0 else 1.0

        translated_coordinates_list: list[ndarray] = []

        # Translate the structure in x, y, and z directions
        for x_translation in np.arange(x_min, x_max + translation_step_x, translation_step_x):
            for y_translation in np.arange(y_min, y_max + translation_step_y, translation_step_y):
                for z_translation in np.arange(z_min, z_max + translation_step_z, translation_step_z):
                    translation_vector: ndarray = np.array([x_translation, y_translation, z_translation])
                    translated_structure: ndarray = coordinates + translation_vector
                    translated_coordinates_list.append(translated_structure)

        # Concatenate all translated coordinates
        translated_coordinates: ndarray = np.vstack(translated_coordinates_list)

        # Remove duplicate coordinates
        translated_coordinates: ndarray = np.unique(translated_coordinates, axis=0)

        return translated_coordinates
