from src.utils import FileConverter, Logger

logger = Logger("Actions")


class AppActionsInitDataParsing:
    @staticmethod
    def convert_init_dat_to_pdb(structure_folder: str, to_set: bool) -> None:
        """
        Convert init_data/{structure_folder}/ljout.dat into result_data/{structure_folder}/ljout-from-init-dat.pdb
        Also create result_data/{structure_folder}/structure_settings.json template if it didn't exist.
        """

        FileConverter.dat_to_pdb(structure_folder=structure_folder)
