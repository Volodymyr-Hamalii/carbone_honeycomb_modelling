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
    def show_init_structure(structure_folder: str, to_set: bool) -> None:
        return AppActionsShowInitData.show_init_structure(structure_folder, to_set)

    @staticmethod
    def show_init_al_structure(structure_folder: str, to_set: bool) -> None:
        return AppActionsShowInitData.show_init_al_structure(structure_folder, to_set)

    @staticmethod
    def show_one_channel_structure(structure_folder: str, to_set: bool) -> None:
        return AppActionsShowInitData.show_one_channel_structure(structure_folder, to_set)

    @staticmethod
    def show_al_in_one_channel_structure(structure_folder: str, to_set: bool) -> None:
        return AppActionsIntercalationAndSorption.show_al_in_one_channel_structure(structure_folder, to_set)

    @staticmethod
    def show_filtered_al_one_channel_structure(structure_folder: str, to_set: bool) -> None:
        return AppActionsIntercalationAndSorption.show_filtered_al_one_channel_structure(structure_folder, to_set)

    @staticmethod
    def full_flow(structure_folder: str, to_set: bool) -> None:
        return AppActionsIntercalationAndSorption.full_flow(structure_folder, to_set)
