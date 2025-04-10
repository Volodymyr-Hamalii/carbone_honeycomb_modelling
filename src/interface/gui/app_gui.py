from pathlib import Path
from tkinter import messagebox
import customtkinter as ctk

from src.utils import Constants, FileReader

from .components import *
from .viewmodels import *
from .windows import *


class AppGui(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Carbon Honeycomb Manager")

        # self.geometry("400x700")
        self.pack_propagate(True)
        self.grid_propagate(True)

        folder_path: Path = Constants.path.INIT_DATA_PATH
        structure_folders: list[str] = FileReader.read_list_of_dirs()

        if not structure_folders:
            messagebox.showerror(
                "Error",
                "Data folders not found. Please, put 'init_data' folder with CH channels data in the root directory:\n"
                f"{folder_path}."
            )
            self.structure_folder: str = "None"
        else:
            self.structure_folder: str = structure_folders[0]

        # Initialize ViewModels
        self.view_model_show_init_data = VMShowInitData()
        self.view_model_data_operations = VMDataConverter()
        self.view_model_intercalation_and_sorption = VMIntercalationAndSorption()

        # Create GUI components
        self.create_widgets(structure_folders)

    def create_widgets(self, structure_folders: list[str]) -> None:
        # Create a dropdown for structure folders
        self.structure_folder_dropdown = DropdownList(
            self,
            options=structure_folders,
            command=self.set_structure_folder,
        )
        self.structure_folder_dropdown.pack(pady=10, padx=10)

        self._create_init_data_info_frame()
        self._create_data_operations_frame()
        self._create_intercalation_and_sorption_frame()

    def set_structure_folder(self, value: str) -> None:
        self.structure_folder: str = value

    ######### Init data info #########

    def _create_init_data_info_frame(self) -> None:
        """ Create a frame for "Init data info" section """
        init_data_info_frame = ctk.CTkFrame(self)
        init_data_info_frame.pack(pady=10, padx=10, fill="x")

        # Add a label for the section
        init_data_info_label = ctk.CTkLabel(init_data_info_frame, text="CH channel general info")
        init_data_info_label.pack(pady=5)

        # Button to show init structure
        self.show_init_structure_btn = Button(
            init_data_info_frame, text="Show initial CH structure", command=self.open_show_init_structure_window
        )
        self.show_init_structure_btn.pack(pady=10, padx=10)

        # Button to show one channel structure
        self.show_one_channel_structure_btn = Button(
            init_data_info_frame, text="Show one channel structure", command=self.open_show_one_channel_structure_window
        )
        self.show_one_channel_structure_btn.pack(pady=10, padx=10)

        # Button to show channel parameters
        self.show_channel_parameters_btn = Button(
            init_data_info_frame, text="Show channel parameters", command=self.open_get_channel_details_window
        )
        self.show_channel_parameters_btn.pack(pady=10, padx=10)

    def open_show_init_structure_window(self) -> None:
        InitDataWindow(self.view_model_show_init_data, self.structure_folder)

    def open_show_one_channel_structure_window(self) -> None:
        InitDataWindow(self.view_model_show_init_data, self.structure_folder, is_one_channel=True)

    def open_get_channel_details_window(self) -> None:
        ChannelDetailsWindow(self.view_model_show_init_data, self.structure_folder)

    ######### Data operations #########

    def _create_data_operations_frame(self) -> None:
        """ Create a frame for "Data operations" section """
        data_operations_frame = ctk.CTkFrame(self)
        data_operations_frame.pack(pady=10, padx=10, fill="x")

        # Add a label for the section
        data_operations_label = ctk.CTkLabel(data_operations_frame, text="Data operations")
        data_operations_label.pack(pady=5)

        # Button to show data converter
        self.show_data_converter_btn = Button(
            data_operations_frame, text="Data converter", command=self.open_data_converter_window
        )
        self.show_data_converter_btn.pack(pady=10, padx=10)

    def open_data_converter_window(self) -> None:
        DataConverterWindow(self.view_model_data_operations, self.structure_folder)

    ######### Intercalation and sorption #########

    def _create_intercalation_and_sorption_frame(self) -> None:
        """ Create a frame for "Intercalation and sorption" section """
        intercalation_and_sorption_frame = ctk.CTkFrame(self)
        intercalation_and_sorption_frame.pack(pady=10, padx=10, fill="x")

        # Add a label for the section
        intercalation_and_sorption_label = ctk.CTkLabel(
            intercalation_and_sorption_frame, text="Al intercalation and sorption")
        intercalation_and_sorption_label.pack(pady=5)

        # Button to show intercalation and sorption
        self.show_intercalation_and_sorption_btn = Button(
            intercalation_and_sorption_frame, text="Build Al atoms in the CH channel",
            command=self.open_update_al_coordinates_table_window)
        self.show_intercalation_and_sorption_btn.pack(pady=10, padx=10)

        # Button to show intercalation and sorption
        self.show_intercalation_and_sorption_btn = Button(
            intercalation_and_sorption_frame, text="Translate Al to other planes",
            command=self.open_translate_al_to_other_planes_window)
        self.show_intercalation_and_sorption_btn.pack(pady=10, padx=10)

        # Button to show intercalation and sorption
        self.show_intercalation_and_sorption_btn = Button(
            intercalation_and_sorption_frame, text="Get Al in channel details",
            command=self.open_get_al_in_channel_details_window)
        self.show_intercalation_and_sorption_btn.pack(pady=10, padx=10)

        # Button to show intercalation and sorption
        self.show_intercalation_and_sorption_btn = Button(
            intercalation_and_sorption_frame, text="Translate Al to all channels",
            command=self.open_translate_al_to_all_channels_window)
        self.show_intercalation_and_sorption_btn.pack(pady=10, padx=10)

    def open_update_al_coordinates_table_window(self) -> None:
        UpdateAlCoordinatesTableWindow(self.view_model_intercalation_and_sorption, self.structure_folder)

    def open_translate_al_to_other_planes_window(self) -> None:
        TranslateAlToOtherPlanesWindow(self.view_model_intercalation_and_sorption, self.structure_folder)

    def open_translate_al_to_all_channels_window(self) -> None:
        TranslateAlToAllChannelsWindow(self.view_model_intercalation_and_sorption, self.structure_folder)

    def open_get_al_in_channel_details_window(self) -> None:
        GetAlInChannelDetailsWindow(self.view_model_intercalation_and_sorption, self.structure_folder)
