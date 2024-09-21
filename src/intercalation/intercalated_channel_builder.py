from pathlib import Path
from numpy import ndarray

from ..utils import PathBuilder, Logger
from ..structure_visualizer import AtomsUniverseBuilder
from ..coordinates_actions import StructureTranslator, PointsOrganizer
from ..data_preparation import StructureSettings

from .al_lattice_type import AlLatticeType
from .al_atoms_filter import AlAtomsFilter


logger = Logger(__name__)


class IntercalatedChannelBuilder:
    @staticmethod
    def build_carbone_coordinates(structure_folder: str, structure_settings: None | StructureSettings) -> ndarray:
        path_to_pdb_file: Path = PathBuilder.build_path_to_result_data_file(structure_folder)

        return AtomsUniverseBuilder.builds_atoms_coordinates(
            path_to_pdb_file,
            channel_coordinate_limits=structure_settings.channel_limits if structure_settings else None)

    @staticmethod
    def build_al_coordinates_for_cell(
            to_translate_al: bool,
            al_file: str,
            structure_settings: None | StructureSettings = None,
    ) -> ndarray:
        path_to_al_pdb_file: Path = PathBuilder.build_path_to_init_data_file(file=al_file)
        coordinates_al: ndarray = AtomsUniverseBuilder.builds_atoms_coordinates(path_to_al_pdb_file)

        if to_translate_al:
            if structure_settings is None:
                raise ValueError(
                    "To translate Al please provide structure_settings.json file with channel_limits.")

            return StructureTranslator.translate_cell(
                cell_coordinates=coordinates_al, translation_limits=structure_settings.channel_limits)

        return coordinates_al

    @staticmethod
    def build_al_coordinates_for_close_packed(
            al_lattice_type: AlLatticeType,
            structure_settings: None | StructureSettings,
            to_translate_al: bool,  # TODO
    ) -> ndarray:

        if structure_settings is None or structure_settings.al_lattice_parameter == 0:
            raise ValueError(
                "To translate Al please provide structure_settings.json file with channel_limits and al_lattice_parameter.")

        if al_lattice_type.is_fcc:
            return AtomsUniverseBuilder.build_fcc_lattice_type(
                lattice_parameter=structure_settings.al_lattice_parameter,
                channel_coordinate_limits=structure_settings.channel_limits,
            )

        if al_lattice_type.is_hcp:
            return AtomsUniverseBuilder.build_hcp_lattice_type(
                lattice_parameter=structure_settings.al_lattice_parameter,
                channel_coordinate_limits=structure_settings.channel_limits,
            )

        logger.warning("No Al atoms built for closed packed structure.")
        return ndarray([])

    @classmethod
    def build_al_in_carbone(
            cls,
            coordinates_carbone: ndarray,
            coordinates_al: ndarray,
            structure_settings: None | StructureSettings,
            to_filter_al_atoms: bool = True,
            equidistant_al_points: bool = True
    ) -> tuple[ndarray, ndarray]:
        """ Return coordinates_carbone, coordinates_al """

        if to_filter_al_atoms is True:
            if structure_settings is None:
                logger.error("Not able to filter AL atoms without structure_settings.json")
            else:
                coordinates_al_filtered: ndarray = AlAtomsFilter.find_max_filtered_atoms(
                    coordinates_carbone=coordinates_carbone,
                    coordinates_al=coordinates_al,
                    structure_settings=structure_settings)

                if equidistant_al_points:
                    # Set Al atoms maximally equidistant from the channel atoms
                    coordinates_al_filtered: ndarray = PointsOrganizer.equidistant_points_sets_in_channel(
                        coordinates_carbone, coordinates_al_filtered)

                return coordinates_carbone, coordinates_al_filtered

        return coordinates_carbone, coordinates_al
