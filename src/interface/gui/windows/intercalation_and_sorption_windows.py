from pathlib import Path
import customtkinter as ctk
from tkinter import messagebox

from src.utils import Logger, FileReader

from ..viewmodels import VMIntercalationAndSorption
from ..components import Button, CheckBox, DropdownList, InputField, InputFieldCoordLimits

__all__: list[str] = [
    "UpdateAlCoordinatesTableWindow",
    "TranslateAlToOtherPlanesWindow",
    "TranslateAlToAllChannelsWindow",
    "GetAlInChannelDetailsWindow",
]


logger = Logger("IntercalationAndSorptionWindow")


class _IntercalationAndSorptionUtils:
    def __init__(self, view_model: VMIntercalationAndSorption, structure_folder: str) -> None:
        self.view_model: VMIntercalationAndSorption = view_model
        self.structure_folder: str = structure_folder

        self.file_names: list[str] = []
        self._refresh_file_name_lists()

    def _refresh_file_name_lists(self) -> None:
        path: Path = self.view_model.data_dir / self.structure_folder
        self.file_names: list[str] = FileReader.read_list_of_files(path, format=".xlsx") or ["None"]
        if (not self.view_model.file_name) or (
                self.view_model.file_name == "None") or (
                self.view_model.file_name not in self.file_names):
            self.view_model.set_file_name(self.file_names[0])


class UpdateAlCoordinatesTableWindow(_IntercalationAndSorptionUtils):
    def __init__(self, view_model: VMIntercalationAndSorption, structure_folder: str) -> None:
        super().__init__(view_model, structure_folder)

        self.create_window()

    def create_window(self) -> None:
        self.input_window = ctk.CTkToplevel()
        self.input_window.pack_propagate(True)
        self.input_window.grid_propagate(True)

        title: str = f"Update Al coordinates table ({self.structure_folder})"
        self.input_window.title(title)

        description: str = (
            "Allows to plot the selected structure, or update the Excel file with Al coordinates "
            "(to remove the empty rows and rebuild the dist matrix)."
            "If the file is not found, the program will build the Al atoms for the provided number of planes."
        )
        description_label: ctk.CTkLabel = ctk.CTkLabel(
            self.input_window, text=description, wraplength=500
        )
        description_label.pack(pady=10, padx=10)

        # Create a frame to hold the columns
        columns_frame = ctk.CTkFrame(self.input_window)
        columns_frame.pack(fill="both", expand=True, padx=10, pady=10)

        # Create frames for left and right columns inside the columns_frame
        left_frame = ctk.CTkFrame(columns_frame)
        right_frame = ctk.CTkFrame(columns_frame)

        left_frame.pack(side="left", fill="both", expand=True, padx=10, pady=10)
        right_frame.pack(side="right", fill="both", expand=True, padx=10, pady=10)

        # Left column inputs
        self.file_names_dropdown: DropdownList = DropdownList(
            left_frame,
            options=self.file_names,
            command=self.view_model.set_file_name,
        )
        self.file_names_dropdown.pack(pady=10, padx=10)

        self.number_of_planes_input_field: InputField = InputField(
            left_frame,
            text="Number of planes",
            command=self.update_number_of_planes,
            default_value=self.view_model.number_of_planes,
        )
        self.number_of_planes_input_field.pack(pady=10, padx=10)

        self.num_of_min_distances_input_field: InputField = InputField(
            left_frame,
            text="Number of min distances for bonds",
            command=self.update_num_of_min_distances,
            default_value=self.view_model.bonds_num_of_min_distances,
        )
        self.num_of_min_distances_input_field.pack(pady=10, padx=10)

        self.bonds_skip_first_distances_input_field = InputField(
            left_frame, text="Skip first distances for bonds",
            command=self.update_bonds_skip_first_distances,
            default_value=self.view_model.bonds_skip_first_distances,
        )
        self.bonds_skip_first_distances_input_field.pack(pady=10, padx=10)

        self.num_of_al_layers_input_field: InputField = InputField(
            left_frame,
            text="Number of Al layers on the plot",
            command=self.update_num_of_al_layers,
            default_value=self.view_model.num_of_al_layers,
        )
        self.num_of_al_layers_input_field.pack(pady=10, padx=10)

        # Right column inputs
        self.coord_x_limits_input_field = InputFieldCoordLimits(
            right_frame, text="X plot limits",
            command=self.update_x_coord_limits,
        )
        self.coord_x_limits_input_field.pack(pady=10, padx=10)

        self.coord_y_limits_input_field = InputFieldCoordLimits(
            right_frame, text="Y plot limits",
            command=self.update_y_coord_limits,
        )
        self.coord_y_limits_input_field.pack(pady=10, padx=10)

        self.coord_z_limits_input_field = InputFieldCoordLimits(
            right_frame, text="Z plot limits",
            command=self.update_z_coord_limits,
        )
        self.coord_z_limits_input_field.pack(pady=10, padx=10)

        self.to_show_al_indexes_checkbox = CheckBox(
            right_frame, text="Show atom's indexes on the plot",
            command=self.update_to_show_al_indexes,
            default=self.view_model.to_show_al_indexes,
        )
        self.to_show_al_indexes_checkbox.pack(pady=10, padx=10)

        # Remaining widgets below the columns
        self.plot_btn: Button = Button(
            self.input_window,
            text="Plot the model",
            command=self.plot_al_in_c_structure,
        )
        self.plot_btn.pack(pady=10, padx=10)

        self.update_tbl_btn: Button = Button(
            self.input_window,
            text="Update the Excel file for plane",
            command=self.update_al_plane_coordinates_file,
        )
        self.update_tbl_btn.pack(pady=10, padx=10)

        self.generate_tbl_btn: Button = Button(
            self.input_window,
            text="Generate the Excel file with Al coordinates for plane",
            command=self.generate_al_plane_coordinates_file,
        )
        self.generate_tbl_btn.pack(pady=(10, 25), padx=10)

    def plot_al_in_c_structure(self) -> None:
        try:
            self.view_model.plot_al_in_c_structure(self.structure_folder)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def update_al_plane_coordinates_file(self) -> None:
        try:
            path_to_file: Path = self.view_model.update_al_plane_coordinates_file(self.structure_folder)
            messagebox.showinfo("Success", f"Al plane coordinates file saved to {path_to_file}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def generate_al_plane_coordinates_file(self) -> None:
        try:
            path_to_file: Path = self.view_model.generate_al_plane_coordinates_file(self.structure_folder)
            self._refresh_file_name_lists()
            messagebox.showinfo("Success", f"Al plane coordinates file saved to {path_to_file}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def update_to_show_al_indexes(self) -> None:
        value = bool(self.to_show_al_indexes_checkbox.get())
        self.view_model.set_to_show_al_indexes(value)

    def update_number_of_planes(self) -> None:
        value = int(self.number_of_planes_input_field.get())
        self.view_model.set_number_of_planes(value)

    def update_num_of_min_distances(self) -> None:
        value = int(self.num_of_min_distances_input_field.get())
        self.view_model.set_bonds_num_of_min_distances(value)

    def update_bonds_skip_first_distances(self) -> None:
        value = int(self.bonds_skip_first_distances_input_field.get())
        self.view_model.set_bonds_skip_first_distances(value)

    def update_num_of_al_layers(self) -> None:
        value = int(self.num_of_al_layers_input_field.get())

        if value > 3:
            messagebox.showerror("Error", "The number of AL layers cannot be greater than 3.")
            return

        self.view_model.set_num_of_al_layers(value)

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


class TranslateAlToOtherPlanesWindow(_IntercalationAndSorptionUtils):
    def __init__(self, view_model: VMIntercalationAndSorption, structure_folder: str) -> None:
        super().__init__(view_model, structure_folder)

        self.create_window()

    def create_window(self) -> None:
        self.input_window = ctk.CTkToplevel()
        self.input_window.pack_propagate(True)
        self.input_window.grid_propagate(True)

        title: str = f"Translate Al to other planes ({self.structure_folder})"
        self.input_window.title(title)

        description: str = (
            "Read the selected file with coordinates for N planes "
            "and translate Al atoms to other planes."
        )
        description_label: ctk.CTkLabel = ctk.CTkLabel(
            self.input_window, text=description, wraplength=500
        )
        description_label.pack(pady=10, padx=10)

        self.file_names_dropdown: DropdownList = DropdownList(
            self.input_window,
            options=self.file_names,
            command=self.view_model.set_file_name,
            title="Al coordinates table to translate",
        )
        self.file_names_dropdown.pack(pady=10, padx=10)

        # Create a frame to hold the columns
        columns_frame = ctk.CTkFrame(self.input_window)
        columns_frame.pack(fill="both", expand=True, padx=10, pady=10)

        # Create frames for left and right columns inside the columns_frame
        left_frame = ctk.CTkFrame(columns_frame)
        right_frame = ctk.CTkFrame(columns_frame)

        left_frame.pack(side="left", fill="both", expand=True, padx=10, pady=10)
        right_frame.pack(side="right", fill="both", expand=True, padx=10, pady=10)

        # Left column inputs
        self.try_to_reflect_al_atoms_checkbox: CheckBox = CheckBox(
            left_frame,
            text="Try to reflect Al atoms to fit the plane",
            command=self.update_to_try_to_reflect_al_atoms,
        )
        self.try_to_reflect_al_atoms_checkbox.pack(pady=10, padx=10)

        self.number_of_planes_input_field: InputField = InputField(
            left_frame,
            text="Number of planes",
            command=self.update_number_of_planes,
            default_value=self.view_model.number_of_planes,
        )
        self.number_of_planes_input_field.pack(pady=10, padx=10)

        self.num_of_min_distances_input_field: InputField = InputField(
            left_frame,
            text="Number of min distances for bonds",
            command=self.update_num_of_min_distances,
            default_value=self.view_model.bonds_num_of_min_distances,
        )
        self.num_of_min_distances_input_field.pack(pady=10, padx=10)

        self.bonds_skip_first_distances_input_field = InputField(
            left_frame, text="Skip first distances for bonds",
            command=self.update_bonds_skip_first_distances,
            default_value=self.view_model.bonds_skip_first_distances,
        )
        self.bonds_skip_first_distances_input_field.pack(pady=10, padx=10)

        self.num_of_al_layers_input_field: InputField = InputField(
            left_frame,
            text="Number of Al layers on the plot",
            command=self.update_num_of_al_layers,
            default_value=self.view_model.num_of_al_layers,
        )
        self.num_of_al_layers_input_field.pack(pady=10, padx=10)

        # Right column inputs

        # Input field for coord_limits
        self.coord_x_limits_input_field = InputFieldCoordLimits(
            right_frame, text="X plot limits",
            command=self.update_x_coord_limits,
        )
        self.coord_x_limits_input_field.pack(pady=10, padx=10)

        self.coord_y_limits_input_field = InputFieldCoordLimits(
            right_frame, text="Y plot limits",
            command=self.update_y_coord_limits,
        )
        self.coord_y_limits_input_field.pack(pady=10, padx=10)

        self.coord_z_limits_input_field = InputFieldCoordLimits(
            right_frame, text="Z plot limits",
            command=self.update_z_coord_limits,
        )
        self.coord_z_limits_input_field.pack(pady=10, padx=10)

        self.to_show_al_indexes_checkbox = CheckBox(
            right_frame, text="Show atom's indexes on the plot",
            command=self.update_to_show_al_indexes,
            default=self.view_model.to_show_al_indexes,
        )
        self.to_show_al_indexes_checkbox.pack(pady=10, padx=10)

        self.translate_btn: Button = Button(
            self.input_window,
            text="Build the model",
            command=self.translate_al_to_other_planes,
        )
        self.translate_btn.pack(pady=10, padx=10)

        self.update_tbl_btn: Button = Button(
            self.input_window,
            text="Update the Excel file",
            command=self.update_file,
        )
        self.update_tbl_btn.pack(pady=(10, 25), padx=10)

    def translate_al_to_other_planes(self) -> None:
        try:
            self.view_model.translate_al_to_other_planes(self.structure_folder)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def update_file(self) -> None:
        try:
            path_to_file: Path = self.view_model.update_al_channel_coordinates(self.structure_folder)
            messagebox.showinfo("Success", f"Al coordinates table saved to {path_to_file}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def update_to_show_al_indexes(self) -> None:
        value = bool(self.to_show_al_indexes_checkbox.get())
        self.view_model.set_to_show_al_indexes(value)

    def update_number_of_planes(self) -> None:
        value = int(self.number_of_planes_input_field.get())
        self.view_model.set_number_of_planes(value)

    def update_num_of_min_distances(self) -> None:
        value = int(self.num_of_min_distances_input_field.get())
        self.view_model.set_bonds_num_of_min_distances(value)

    def update_bonds_skip_first_distances(self) -> None:
        value = int(self.bonds_skip_first_distances_input_field.get())
        self.view_model.set_bonds_skip_first_distances(value)

    def update_to_try_to_reflect_al_atoms(self) -> None:
        value = bool(self.try_to_reflect_al_atoms_checkbox.get())
        self.view_model.set_to_try_to_reflect_al_atoms(value)

    def update_num_of_al_layers(self) -> None:
        value = int(self.num_of_al_layers_input_field.get())

        if value > 3:
            self.num_of_al_layers_input_field.label.configure(text="The number of AL layers cannot be greater than 3.")
            return

        self.view_model.set_num_of_al_layers(value)

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


class GetAlInChannelDetailsWindow(_IntercalationAndSorptionUtils):
    def __init__(self, view_model: VMIntercalationAndSorption, structure_folder: str) -> None:
        super().__init__(view_model, structure_folder)

        self.create_window()

    def create_window(self) -> None:
        self.input_window = ctk.CTkToplevel()
        self.input_window.pack_propagate(True)
        self.input_window.grid_propagate(True)

        title: str = f"Get Al in channel details ({self.structure_folder})"
        self.input_window.title(title)

        description: str = (
            "Write an Excel file with Al atoms in channel details (Al atoms coordinates, "
            "distances to the plane, distances to the carbon atoms, distances to the other Al atoms)."
        )
        description_label: ctk.CTkLabel = ctk.CTkLabel(
            self.input_window, text=description, wraplength=400
        )
        description_label.pack(pady=10, padx=10)

        self.file_names_dropdown: DropdownList = DropdownList(
            self.input_window,
            options=self.file_names,
            command=self.view_model.set_file_name,
            title="Al coordinates table to analyze",
        )
        self.file_names_dropdown.pack(pady=10, padx=10)

        self.get_details_btn: Button = Button(
            self.input_window,
            text="Save Al in channel details to Excel file",
            command=self.get_al_in_channel_details,
        )
        self.get_details_btn.pack(pady=(10, 25), padx=10)

    def get_al_in_channel_details(self) -> None:
        try:
            path_to_file: Path = self.view_model.get_al_in_channel_details(self.structure_folder)
            messagebox.showinfo("Success", f"Al in channel details saved to {path_to_file}")
        except Exception as e:
            messagebox.showerror("Error", str(e))


class TranslateAlToAllChannelsWindow(_IntercalationAndSorptionUtils):
    def __init__(self, view_model: VMIntercalationAndSorption, structure_folder: str) -> None:
        super().__init__(view_model, structure_folder)

        self.create_window()

    def create_window(self) -> None:
        self.input_window = ctk.CTkToplevel()
        title: str = f"Translate Al to all channels ({self.structure_folder})"
        self.input_window.title(title)
        self.input_window.pack_propagate(True)
        self.input_window.grid_propagate(True)

        self.file_names_dropdown: DropdownList = DropdownList(
            self.input_window,
            options=self.file_names,
            command=self.view_model.set_file_name,
            title="Al coordinates table to translate",
        )
        self.file_names_dropdown.pack(pady=10, padx=10)

        # Create a frame to hold the columns
        columns_frame = ctk.CTkFrame(self.input_window)
        columns_frame.pack(fill="both", expand=True, padx=10, pady=10)

        # Create frames for left and right columns inside the columns_frame
        left_frame = ctk.CTkFrame(columns_frame)
        right_frame = ctk.CTkFrame(columns_frame)

        left_frame.pack(side="left", fill="both", expand=True, padx=10, pady=10)
        right_frame.pack(side="right", fill="both", expand=True, padx=10, pady=10)

        # Left column inputs
        self.num_of_min_distances_input_field: InputField = InputField(
            left_frame,
            text="Number of min distances for bonds",
            command=self.update_num_of_min_distances,
            default_value=self.view_model.bonds_num_of_min_distances,
        )
        self.num_of_min_distances_input_field.pack(pady=10, padx=10)

        self.bonds_skip_first_distances_input_field = InputField(
            left_frame, text="Skip first distances for bonds",
            command=self.update_bonds_skip_first_distances,
            default_value=self.view_model.bonds_skip_first_distances,
        )
        self.bonds_skip_first_distances_input_field.pack(pady=10, padx=10)

        self.num_of_al_layers_input_field: InputField = InputField(
            left_frame,
            text="Number of Al layers on the plot",
            command=self.update_num_of_al_layers,
            default_value=self.view_model.num_of_al_layers,
        )
        self.num_of_al_layers_input_field.pack(pady=10, padx=10)

        self.try_to_reflect_al_atoms_checkbox: CheckBox = CheckBox(
            left_frame,
            text="Try to reflect Al atoms to fit the plane\n(if no init file and the Al atoms will be calculated)",
            command=self.update_to_try_to_reflect_al_atoms,
        )
        self.try_to_reflect_al_atoms_checkbox.pack(pady=10, padx=10)

        # Right column inputs
        self.coord_x_limits_input_field = InputFieldCoordLimits(
            right_frame, text="X plot limits",
            command=self.update_x_coord_limits,
        )
        self.coord_x_limits_input_field.pack(pady=10, padx=10)

        self.coord_y_limits_input_field = InputFieldCoordLimits(
            right_frame, text="Y plot limits",
            command=self.update_y_coord_limits,
        )
        self.coord_y_limits_input_field.pack(pady=10, padx=10)

        self.coord_z_limits_input_field = InputFieldCoordLimits(
            right_frame, text="Z plot limits",
            command=self.update_z_coord_limits,
        )
        self.coord_z_limits_input_field.pack(pady=10, padx=10)

        self.to_show_al_indexes_checkbox = CheckBox(
            right_frame, text="Show atom's indexes on the plot",
            command=self.update_to_show_al_indexes,
            default=self.view_model.to_show_al_indexes,
        )
        self.to_show_al_indexes_checkbox.pack(pady=10, padx=10)

        self.translate_btn: Button = Button(
            self.input_window,
            text="Generate output files",
            command=self.translate_al_to_all_channels_generate_files,
        )
        self.translate_btn.pack(pady=10, padx=10)

        self.translate_btn: Button = Button(
            self.input_window,
            text="Plot structure",
            command=self.plot_al_in_c_structure,
        )
        self.translate_btn.pack(pady=(10, 25), padx=10)

    def plot_al_in_c_structure(self) -> None:
        try:
            self.view_model.translate_al_to_all_channels_plot(self.structure_folder)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def translate_al_to_all_channels_generate_files(self) -> None:
        try:
            paths: tuple[Path, Path, Path] = self.view_model.translate_al_to_all_channels_generate_files(
                self.structure_folder)
            messagebox.showinfo("Success", f"Generated files:\n{paths[0].name}\n{paths[1].name}\n{paths[2].name}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def update_file_name(self, value: str) -> None:
        self.view_model.set_file_name(value)

    def update_to_show_al_indexes(self) -> None:
        value = bool(self.to_show_al_indexes_checkbox.get())
        self.view_model.set_to_show_al_indexes(value)

    def update_num_of_min_distances(self) -> None:
        value = int(self.num_of_min_distances_input_field.get())
        self.view_model.set_bonds_num_of_min_distances(value)

    def update_bonds_skip_first_distances(self) -> None:
        value = int(self.bonds_skip_first_distances_input_field.get())
        self.view_model.set_bonds_skip_first_distances(value)

    def update_to_try_to_reflect_al_atoms(self) -> None:
        value = bool(self.try_to_reflect_al_atoms_checkbox.get())
        self.view_model.set_to_try_to_reflect_al_atoms(value)

    def update_num_of_al_layers(self) -> None:
        value = int(self.num_of_al_layers_input_field.get())
        self.view_model.set_num_of_al_layers(value)

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
