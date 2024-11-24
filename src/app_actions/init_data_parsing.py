from src.utils import FileConverter, Logger
from src.data_preparation import StructureSettingsManager

logger = Logger("Actions")


class AppActionsInitDataParsing:
    @staticmethod
    def convert_init_dat_to_pdb(structure_folder: str, to_set: bool) -> None:
        """
        Convert init_data/{structure_folder}/ljout.dat into result_data/{structure_folder}/ljout-from-init-dat.pdb
        Also create result_data/{structure_folder}/structure_settings.json template if it didn't exist.
        """

        FileConverter.dat_to_pdb(structure_folder=structure_folder)

        # Create template for coordinates if it doesn't exists
        StructureSettingsManager.create_structure_settings_template(structure_folder=structure_folder)
