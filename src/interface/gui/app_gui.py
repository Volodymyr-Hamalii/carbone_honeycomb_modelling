import customtkinter as ctk

from src.utils import FileReader

from .viewmodels.show_init_data import ViewModelShowInitData
from .components import Button, DropdownList
from .windows import StructureWindow, ChannelDetailsWindow

class AppGui(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Carbon Honeycomb Manager")
        self.geometry("600x400")

        structure_folders: list[str] = FileReader.read_list_of_dirs()
        self.structure_folder: str = structure_folders[0]

        # Initialize ViewModel
        self.view_model = ViewModelShowInitData(self.structure_folder)

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

        # Button to show init structure
        self.show_init_structure_btn = Button(
            self, text="Show init CH structure", command=self.open_show_init_structure_window
        )
        self.show_init_structure_btn.pack(pady=10)

        # Button to show one channel structure
        self.show_one_channel_structure_btn = Button(
            self, text="Show one channel structure", command=self.open_show_one_channel_structure_window
        )
        self.show_one_channel_structure_btn.pack(pady=10)

        # Button to show channel parameters
        self.show_channel_parameters_btn = Button(
            self, text="Show channel parameters", command=self.open_get_channel_details_window
        )
        self.show_channel_parameters_btn.pack(pady=10)

    def set_structure_folder(self, value: str) -> None:
        self.structure_folder: str = value
        self.view_model.set_structure_folder(value)

    def open_show_init_structure_window(self) -> None:
        StructureWindow(self.view_model, self.structure_folder)

    def open_show_one_channel_structure_window(self) -> None:
        StructureWindow(self.view_model, self.structure_folder, one_channel=True)

    def open_get_channel_details_window(self) -> None:
        ChannelDetailsWindow(self.view_model)
