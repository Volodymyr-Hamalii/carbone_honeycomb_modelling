import customtkinter as ctk

from src.utils import FileReader

from .viewmodels import VMShowInitData, VMDataConverter
from .components import Button, DropdownList
from .windows import *


class AppGui(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Carbon Honeycomb Manager")
        self.geometry("600x400")

        structure_folders: list[str] = FileReader.read_list_of_dirs()
        self.structure_folder: str = structure_folders[0]

        # Initialize ViewModels
        self.view_model_show_init_data = VMShowInitData()
        self.view_model_data_operations = VMDataConverter()

        # Create GUI components
        self.create_widgets(structure_folders)

    def create_widgets(self, structure_folders: list[str]) -> None:
        # Create a dropdown for structure folders
        self.structure_folder_dropdown = DropdownList(
            self,
            options=structure_folders,
            command=self.set_structure_folder,
        )
        self.structure_folder_dropdown.pack(pady=10)

        self._create_init_data_info_frame()
        self._create_data_operations_frame()

    def _create_init_data_info_frame(self) -> None:
        """ Create a frame for "Init data info" section """
        init_data_info_frame = ctk.CTkFrame(self)
        init_data_info_frame.pack(pady=10, fill="x")

        # Add a label for the section
        init_data_info_label = ctk.CTkLabel(init_data_info_frame, text="CH channel general info")
        init_data_info_label.pack(pady=5)

        # Button to show init structure
        self.show_init_structure_btn = Button(
            init_data_info_frame, text="Show init CH structure", command=self.open_show_init_structure_window
        )
        self.show_init_structure_btn.pack(pady=10)

        # Button to show one channel structure
        self.show_one_channel_structure_btn = Button(
            init_data_info_frame, text="Show one channel structure", command=self.open_show_one_channel_structure_window
        )
        self.show_one_channel_structure_btn.pack(pady=10)

        # Button to show channel parameters
        self.show_channel_parameters_btn = Button(
            init_data_info_frame, text="Show channel parameters", command=self.open_get_channel_details_window
        )
        self.show_channel_parameters_btn.pack(pady=10)

    def _create_data_operations_frame(self) -> None:
        """ Create a frame for "Data operations" section """
        data_operations_frame = ctk.CTkFrame(self)
        data_operations_frame.pack(pady=10, fill="x")

        # Add a label for the section
        data_operations_label = ctk.CTkLabel(data_operations_frame, text="Data operations")
        data_operations_label.pack(pady=5)

        # Button to show data converter
        self.show_data_converter_btn = Button(
            data_operations_frame, text="Data converter", command=self.open_data_converter_window
        )
        self.show_data_converter_btn.pack(pady=10)

    def set_structure_folder(self, value: str) -> None:
        self.structure_folder: str = value

    def open_show_init_structure_window(self) -> None:
        StructureWindow(self.view_model_show_init_data, self.structure_folder)

    def open_show_one_channel_structure_window(self) -> None:
        StructureWindow(self.view_model_show_init_data, self.structure_folder, one_channel=True)

    def open_get_channel_details_window(self) -> None:
        ChannelDetailsWindow(self.view_model_show_init_data, self.structure_folder)

    def open_data_converter_window(self) -> None:
        DataConverterWindow(self.view_model_data_operations, self.structure_folder)
