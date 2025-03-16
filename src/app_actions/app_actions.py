from src.utils import Logger

from .init_data_parsing import AppActionsInitDataParsing
from .show_init_data import AppActionsShowInitData
from .intercalation_and_sorption import AppActionsIntercalationAndSorption


logger = Logger("Actions")


class AppActions:
    @classmethod
    def help(cls, structure_folder: str) -> None:
        """ Print all available methods and general info. """

        # Get all attributes of the class
        attributes: list[str] = dir(cls)

        # Filter out only methods
        methods: list[str] = [attr for attr in attributes if callable(getattr(cls, attr)) and not attr.startswith("__")]

        message: str = "Available actions:\n"
        for method in methods:
            message = f"{message}* {method}\n"

        message = message + "\nmain.py takes the following parameters: action, structure_folder, set.\n" + \
            "   1) action -- any action from the list above,\n" + \
            "   2) structure_folder -- structure to process from 'init_data' or 'result_data' folders,\n" + \
            "   3) set -- parameters custumization (to build bonds, to filter atoms etc.); just put 'set' to customize building,\n" + \
            "   4) optional arguments -- specific parameters if some method needs them.\n" + \
            f"For example:\npython3 main.py show_init_structure {structure_folder} set\n" + \
            "If you don't want to specify some argument -- just provide '_' on this place."

        logger.info(message)

    @staticmethod
    def convert_init_dat_to_pdb(structure_folder: str, to_set: bool) -> None:
        return AppActionsInitDataParsing.convert_init_dat_to_pdb(structure_folder, to_set)

    @staticmethod
    def convert_excel_to_dat(structure_folder: str, to_set: bool) -> None:
        return AppActionsInitDataParsing.convert_excel_to_dat(structure_folder, to_set)

    @staticmethod
    def show_init_structure(structure_folder: str, to_set: bool) -> None:
        return AppActionsShowInitData.show_init_structure(structure_folder, to_set)

    @staticmethod
    def get_channel_details(structure_folder: str, to_set: bool) -> None:
        return AppActionsShowInitData.get_channel_details(structure_folder, to_set)

    @staticmethod
    def show_init_al_structure(structure_folder: str, to_set: bool) -> None:
        return AppActionsShowInitData.show_init_al_structure(structure_folder, to_set)

    @staticmethod
    def show_one_channel_structure(structure_folder: str, to_set: bool) -> None:
        return AppActionsShowInitData.show_one_channel_structure(structure_folder, to_set)

    @staticmethod
    def fill_all_channels(structure_folder: str, to_set: bool) -> None:
        return AppActionsIntercalationAndSorption.fill_all_channels(structure_folder, to_set)

    # @staticmethod
    # def show_filtered_al_one_channel_structure(structure_folder: str, to_set: bool) -> None:
    #     return AppActionsIntercalationAndSorption.show_filtered_al_one_channel_structure(structure_folder, to_set)

    @staticmethod
    def update_al_coordinates_tbl(structure_folder: str, to_set: bool) -> None:
        AppActionsIntercalationAndSorption.update_al_coordinates_tbl(structure_folder, to_set)

    @staticmethod
    def full_flow(structure_folder: str, to_set: bool) -> None:
        return AppActionsIntercalationAndSorption.full_flow(structure_folder, to_set)

    @staticmethod
    def translate_al_to_other_planes(structure_folder: str, to_set: bool) -> None:
        return AppActionsIntercalationAndSorption.translate_al_to_other_planes(structure_folder, to_set)

    @staticmethod
    def get_al_in_channel_details(structure_folder: str, to_set: bool) -> None:
        return AppActionsIntercalationAndSorption.get_al_in_channel_details(structure_folder, to_set)

    @staticmethod
    def translate_al_to_all_channels(structure_folder: str, to_set: bool) -> None:
        return AppActionsIntercalationAndSorption.translate_al_to_all_channels(structure_folder, to_set)

    @staticmethod
    def update_al_full_channel_coordinates_tbl(structure_folder: str, to_set: bool) -> None:
        return AppActionsIntercalationAndSorption.update_al_full_channel_coordinates_tbl(structure_folder, to_set)
