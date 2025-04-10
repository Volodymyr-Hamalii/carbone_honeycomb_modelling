from pathlib import Path
import customtkinter as ctk

from src.utils import Constants, FileReader
from ..viewmodels import VMShowInitData
from ..components import (
    Button,
    CheckBox,
    DropdownList,
    InputField,
    InputFieldCoordLimits,
)


class InitDataWindow:
    def __init__(self, view_model: VMShowInitData, structure_folder: str, is_one_channel: bool = False) -> None:
        self.view_model: VMShowInitData = view_model
        self.view_model.set_data_dir(Constants.filenames.INIT_DATA_DIR)

        self.structure_folder: str = structure_folder
        self.is_one_channel: bool = is_one_channel

        self.file_names: list[str] = []
        self._refresh_file_name_lists()

        self.create_window()

    def _refresh_file_name_lists(self) -> None:
        path: Path = self.view_model.data_dir / self.structure_folder
        self.file_names: list[str] = FileReader.read_list_of_files(path) or ["None"]
        if not self.view_model.file_name or self.view_model.file_name == "None":
            self.view_model.set_file_name(self.file_names[0])

    def create_window(self) -> None:
        self.input_window = ctk.CTkToplevel()
        self.input_window.pack_propagate(True)
        self.input_window.grid_propagate(True)

        title: str = (
            f"Show one channel structure ({self.structure_folder})"
            if self.is_one_channel
            else f"Show init full CH structure ({self.structure_folder}) "
        )
        self.input_window.title(title)

        self.file_names_dropdown: DropdownList = DropdownList(
            self.input_window,
            options=self.file_names,
            command=self.view_model.set_file_name,
        )
        self.file_names_dropdown.pack(pady=10, padx=10)

        # Checkbox for to_build_bonds
        self.to_build_bonds_checkbox = CheckBox(
            self.input_window, text="Build Bonds",
            command=self.update_to_build_bonds,
            default=self.view_model.to_build_bonds,
        )
        self.to_build_bonds_checkbox.pack(pady=10, padx=10)

        # Checkbox for to_show_coordinates
        self.to_show_coordinates_checkbox = CheckBox(
            self.input_window, text="Show coordinates",
            command=self.update_to_show_coordinates,
            default=self.view_model.to_show_coordinates,
        )
        self.to_show_coordinates_checkbox.pack(pady=10, padx=10)

        # Checkbox for to_show_indexes
        self.to_show_c_indexes_checkbox = CheckBox(
            self.input_window, text="Show C atoms indexes",
            command=self.update_to_show_c_indexes,
            default=self.view_model.to_show_c_indexes,
        )
        self.to_show_c_indexes_checkbox.pack(pady=10, padx=10)

        # Input field for bonds_num_of_min_distances
        self.bonds_num_of_min_distances_input_field = InputField(
            self.input_window, text="Number of min distances for bonds",
            command=self.update_bonds_num_of_min_distances,
            default_value=self.view_model.bonds_num_of_min_distances,
        )
        self.bonds_num_of_min_distances_input_field.pack(pady=10, padx=10)

        # Input field for bonds_skip_first_distances
        self.bonds_skip_first_distances_input_field = InputField(
            self.input_window, text="Skip first distances for bonds",
            command=self.update_bonds_skip_first_distances,
            default_value=self.view_model.bonds_skip_first_distances,
        )
        self.bonds_skip_first_distances_input_field.pack(pady=10, padx=10)

        # Input field for coord_limits
        self.coord_x_limits_input_field = InputFieldCoordLimits(
            self.input_window, text="X plot limits",
            command=self.update_x_coord_limits,
        )
        self.coord_x_limits_input_field.pack(pady=10, padx=10)

        self.coord_y_limits_input_field = InputFieldCoordLimits(
            self.input_window, text="Y plot limits",
            command=self.update_y_coord_limits,
        )
        self.coord_y_limits_input_field.pack(pady=10, padx=10)

        self.coord_z_limits_input_field = InputFieldCoordLimits(
            self.input_window, text="Z plot limits",
            command=self.update_z_coord_limits,
        )
        self.coord_z_limits_input_field.pack(pady=10, padx=10)

        # Button to proceed to the next step
        self.next_btn = Button(
            self.input_window, text="Show structure",
            command=self.show_structure,
        )
        self.next_btn.pack(pady=(10, 25), padx=10)

    def show_structure(self) -> None:
        if self.is_one_channel:
            self.view_model.show_one_channel_structure(self.structure_folder)
        else:
            self.view_model.show_init_structure(self.structure_folder)
        # self.show_plot_window()

    def update_to_build_bonds(self) -> None:
        value = bool(self.to_build_bonds_checkbox.get())
        self.view_model.set_to_build_bonds(value)

    def update_to_show_coordinates(self) -> None:
        value = bool(self.to_show_coordinates_checkbox.get())
        self.view_model.set_to_show_coordinates(value)

    def update_to_show_c_indexes(self) -> None:
        value = bool(self.to_show_c_indexes_checkbox.get())
        self.view_model.set_to_show_c_indexes(value)

    def update_bonds_num_of_min_distances(self) -> None:
        value = int(self.bonds_num_of_min_distances_input_field.get())
        self.view_model.set_bonds_num_of_min_distances(value)

    def update_bonds_skip_first_distances(self) -> None:
        value = int(self.bonds_skip_first_distances_input_field.get())
        self.view_model.set_bonds_skip_first_distances(value)

    def update_x_coord_limits(self) -> None:
        value_min: str = self.coord_x_limits_input_field.min_entry.get()
        value_max: str = self.coord_x_limits_input_field.max_entry.get()
        self.view_model.set_x_min(value_min)
        self.view_model.set_x_max(value_max)

    def update_y_coord_limits(self) -> None:
        value_min: str = self.coord_y_limits_input_field.min_entry.get()
        value_max: str = self.coord_y_limits_input_field.max_entry.get()
        self.view_model.set_y_min(value_min)
        self.view_model.set_y_max(value_max)

    def update_z_coord_limits(self) -> None:
        value_min: str = self.coord_z_limits_input_field.min_entry.get()
        value_max: str = self.coord_z_limits_input_field.max_entry.get()
        self.view_model.set_z_min(value_min)
        self.view_model.set_z_max(value_max)

    # def show_plot_window(self) -> None:
    #     plot_window = PlotWindow(self.input_window, title=self.structure_folder)
    #     plot_window.destroy()
