import customtkinter as ctk
from src.interface.gui.viewmodels.show_init_data import ViewModelShowInitData
from src.interface.gui.components import Button, CheckBox, PlotWindow, DropdownList

from src.utils import FileReader


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

    def create_widgets(self, structure_folders: list[str]):
        # Create a dropdown for structure folders
        self.structure_folder_dropdown = DropdownList(
            self,
            text="Select structure folder",
            options=structure_folders,
            command=self.set_structure_folder,
        )
        self.structure_folder_dropdown.pack(pady=10)

        # Button to show init structure
        self.show_init_structure_btn = Button(
            self, text="Show init CH structure", command=self.open_show_init_structure_window
        )
        self.show_init_structure_btn.pack(pady=10)

        # Button to show init structure
        self.show_one_channel_structure_btn = Button(
            self, text="Show one channel structure", command=self.open_show_one_channel_structure_window
        )
        self.show_one_channel_structure_btn.pack(pady=10)

        # Button to show init structure
        self.show_one_channel_structure_btn = Button(
            self, text="Show one channel structure", command=self.open_show_one_channel_structure_window
        )
        self.show_one_channel_structure_btn.pack(pady=10)

    def set_structure_folder(self, value: str) -> None:
        self.structure_folder: str = value
        self.view_model.set_structure_folder(value)

    def open_show_init_structure_window(self) -> None:
        # Create a new window for input parameters
        self.input_window = ctk.CTkToplevel(self)
        self.input_window.title("Input Parameters")
        self.input_window.geometry("300x200")

        # Checkbox for to_build_bonds
        self.to_build_bonds_checkbox = CheckBox(
            self.input_window, text="Build Bonds", command=self.update_to_build_bonds)
        self.to_build_bonds_checkbox.pack(pady=10)

        # Checkbox for to_show_coordinates
        self.to_show_coordinates_checkbox = CheckBox(
            self.input_window, text="Show coordinates", command=self.update_to_show_coordinates)
        self.to_show_coordinates_checkbox.pack(pady=10)

        # Checkbox for to_show_indexes
        self.to_show_indexes_checkbox = CheckBox(
            self.input_window, text="Show indexes", command=self.update_to_show_indexes)
        self.to_show_indexes_checkbox.pack(pady=10)

        # Button to proceed to the next step
        self.next_btn = Button(
            self.input_window, text="Show structure", command=self.show_init_structure)
        self.next_btn.pack(pady=10)

    def show_init_structure(self) -> None:
        # Close the input window
        self.input_window.destroy()

        self.view_model.show_init_structure()
        self.show_plot_window()

    def open_show_one_channel_structure_window(self) -> None:
        # Create a new window for input parameters
        self.input_window = ctk.CTkToplevel(self)
        self.input_window.title("Input Parameters")
        self.input_window.geometry("300x200")

        # Checkbox for to_build_bonds
        self.to_build_bonds_checkbox = CheckBox(
            self.input_window, text="Build Bonds", command=self.update_to_build_bonds)
        self.to_build_bonds_checkbox.pack(pady=10)

        # Checkbox for to_show_coordinates
        self.to_show_coordinates_checkbox = CheckBox(
            self.input_window, text="Show coordinates", command=self.update_to_show_coordinates)
        self.to_show_coordinates_checkbox.pack(pady=10)

        # Checkbox for to_show_indexes
        self.to_show_indexes_checkbox = CheckBox(
            self.input_window, text="Show indexes", command=self.update_to_show_indexes)
        self.to_show_indexes_checkbox.pack(pady=10)

        # Button to proceed to the next step
        self.next_btn = Button(
            self.input_window, text="Show structure", command=self.show_one_channel_structure)
        self.next_btn.pack(pady=10)

    def show_one_channel_structure(self) -> None:
        # Close the input window
        self.input_window.destroy()

        self.view_model.show_one_channel_structure()
        self.show_plot_window()

    def open_get_channel_details_window(self) -> None:
        # Create a new window for input parameters
        self.input_window = ctk.CTkToplevel(self)
        self.input_window.title("Input Parameters")
        self.input_window.geometry("300x200")

        # Checkbox for to_build_bonds
        self.to_build_bonds_checkbox = CheckBox(
            self.input_window, text="Build Bonds", command=self.update_to_build_bonds)
        self.to_build_bonds_checkbox.pack(pady=10)

        # Checkbox for to_show_coordinates
        self.to_show_coordinates_checkbox = CheckBox(
            self.input_window, text="Show coordinates", command=self.update_to_show_coordinates)
        self.to_show_coordinates_checkbox.pack(pady=10)

        # Button to proceed to the next step
        self.next_btn = Button(
            self.input_window, text="Show channel parameters", command=self.get_channel_details)
        self.next_btn.pack(pady=10)

    def get_channel_details(self) -> None:
        # Close the input window
        self.input_window.destroy()

        self.view_model.get_channel_details()
        self.show_plot_window()

    def update_to_build_bonds(self) -> None:
        value = bool(self.to_build_bonds_checkbox.get())
        self.view_model.set_to_build_bonds(value)

    def update_to_show_coordinates(self) -> None:
        value = bool(self.to_show_coordinates_checkbox.get())
        self.view_model.set_to_show_coordinates(value)

    def update_to_show_indexes(self) -> None:
        value = bool(self.to_show_indexes_checkbox.get())
        self.view_model.set_to_show_indexes(value)

    def show_plot_window(self) -> None:
        # Create a new window to display the plot
        plot_window = PlotWindow(self, title=self.structure_folder)

        # Close the plot window
        plot_window.destroy()


app = AppGui()
app.mainloop()
