from pathlib import Path
from numpy import ndarray

from src.utils import Constants, PathBuilder, Logger, Inputs
from src.structure_visualizer import StructureVisualizer, AtomsUniverseBuilder, VisualizationParameters
from src.data_preparation import StructureSettings, StructureSettingsManager, ChannelLimits
from src.projects.intercalation_and_sorption import IntercalatedChannelBuilder
from src.base_structure_classes import AlLatticeType

logger = Logger("Actions")


class AppActionsShowInitData:

    @staticmethod
    def show_init_structure(structure_folder: str, to_set: bool) -> None:
        """ Show 3D model of result_data/{structure_folder}/ljout-from-init-dat.pdb """

        path_to_init_pdb_file: Path = PathBuilder.build_path_to_result_data_file(structure_folder)
        coordinates: ndarray = AtomsUniverseBuilder.builds_atoms_coordinates(path_to_init_pdb_file)

        to_build_bonds: bool = Inputs.bool_input(to_set, default_value=True, text="To build bonds between atoms")
        StructureVisualizer.show_structure(
            coordinates, to_build_bonds=to_build_bonds, set_equal_scale=False, title=structure_folder)

    @staticmethod
    def show_init_al_structure(structure_folder: str, to_set: bool) -> None:
        """ Show 3D model of init_data/al.pdb """

        to_translate_al: bool = Inputs.bool_input(
            to_set, default_value=True, text="To translate AL atomes to fill full volume")

        al_lattice_type_str: str = Inputs.text_input(
            to_set, default_value="FCC",
            # to_set, default_value="HCP",
            text=AlLatticeType.get_info(),
            available_values=AlLatticeType.get_available_types())
        al_lattice_type = AlLatticeType(al_lattice_type_str)

        structure_settings: None | StructureSettings = StructureSettingsManager.get_structure_settings(
            structure_folder=structure_folder)

        if al_lattice_type.is_cell:
            al_file: str = Inputs.text_input(to_set, default_value=Constants.filenames.AL_FILE, text="Init AL file")

            coordinates_al: ndarray = IntercalatedChannelBuilder.build_al_coordinates_for_cell(
                to_translate_al=to_translate_al,
                al_file=al_file,
                structure_settings=structure_settings)

            num_of_min_distances = 1
            skip_first_distances = 1
        else:
            # Fill the volume with aluminium for close-packed lattice
            coordinates_al: ndarray = IntercalatedChannelBuilder.build_al_coordinates_for_close_packed(
                al_lattice_type=al_lattice_type,
                structure_settings=structure_settings,
                to_translate_al=to_translate_al)

            num_of_min_distances = 1
            skip_first_distances = 0

        to_build_bonds: bool = Inputs.bool_input(to_set, default_value=True, text="To build bonds between atoms")
        StructureVisualizer.show_structure(
            coordinates=coordinates_al,
            to_build_bonds=to_build_bonds,
            visual_parameters=VisualizationParameters.al,
            num_of_min_distances=num_of_min_distances,
            skip_first_distances=skip_first_distances,
            title="Aluminium")

    @staticmethod
    def show_one_channel_structure(structure_folder: str, to_set: bool) -> None:
        """
        Build one channel model from result_data/{structure_folder}/ljout-from-init-dat.pdb atoms
        based on result_data/{structure_folder}/structure_settings.json channel limits.

        Write result to result_data/{structure_folder}/ljout-result-one-channel.pdb if it didn't exist.
        """

        path_to_init_pdb_file: Path = PathBuilder.build_path_to_result_data_file(structure_folder)

        structure_settings: None | StructureSettings = StructureSettingsManager.get_structure_settings(
            structure_folder=structure_folder)

        channel_limits: ChannelLimits | None = structure_settings.channel_limits if structure_settings else None

        coordinates: ndarray = AtomsUniverseBuilder.builds_atoms_coordinates(
            path_to_init_pdb_file, channel_limits)

        to_build_bonds: bool = Inputs.bool_input(to_set, default_value=True, text="To build bonds between atoms")
        StructureVisualizer.show_structure(coordinates, to_build_bonds=to_build_bonds, title=structure_folder)
