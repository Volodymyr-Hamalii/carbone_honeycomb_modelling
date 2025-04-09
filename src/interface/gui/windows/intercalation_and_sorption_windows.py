from pathlib import Path
import customtkinter as ctk
from tkinter import messagebox

from src.utils import Logger, Constants, FileReader

from ..viewmodels import VMIntercalationAndSorption
from ..components import Button, CheckBox, DropdownList, InputField

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
        if self.view_model.file_name == "None" and self.file_names:
            self.view_model.set_file_name(self.file_names[0])


class UpdateAlCoordinatesTableWindow(_IntercalationAndSorptionUtils):
    def __init__(self, view_model: VMIntercalationAndSorption, structure_folder: str) -> None:
        super().__init__(view_model, structure_folder)

        self.create_window()

    def create_window(self) -> None:
        self.input_window = ctk.CTkToplevel()
        title: str = f"Update Al coordinates table ({self.structure_folder})"
        self.input_window.title(title)
        # self.input_window.geometry("450x600")
        self.input_window.pack_propagate(True)
        self.input_window.grid_propagate(True)

        self.file_names_dropdown: DropdownList = DropdownList(
            self.input_window,
            options=self.file_names,
            command=self.view_model.set_file_name,
            # title="Al coordinates table to update",
        )
        self.file_names_dropdown.pack(pady=10, padx=10)

        self.number_of_planes_input_field: InputField = InputField(
            self.input_window,
            text="Number of planes",
            command=self.update_number_of_planes,
            # title="Number of planes",
            default_value=self.view_model.number_of_planes,
        )
        self.number_of_planes_input_field.pack(pady=10, padx=10)

        self.num_of_min_distances_input_field: InputField = InputField(
            self.input_window,
            text="Number of min distances for bonds",
            command=self.update_num_of_min_distances,
            # title="Number of min distances for bonds",
            default_value=self.view_model.bonds_num_of_min_distances,
        )
        self.num_of_min_distances_input_field.pack(pady=10, padx=10)

        # Input field for bonds_skip_first_distances
        self.bonds_skip_first_distances_input_field = InputField(
            self.input_window, text="Skip first distances for bonds",
            command=self.update_bonds_skip_first_distances,
            default_value=self.view_model.bonds_skip_first_distances,
        )
        self.bonds_skip_first_distances_input_field.pack(pady=10, padx=10)

        self.num_of_al_layers_input_field: InputField = InputField(
            self.input_window,
            text="Number of Al layers on the plot",
            command=self.update_num_of_al_layers,
            # title="Number of AL layers",
            default_value=self.view_model.num_of_al_layers,
        )
        self.num_of_al_layers_input_field.pack(pady=10, padx=10)

        # Checkbox for to_show_indexes
        self.to_show_al_indexes_checkbox = CheckBox(
            self.input_window, text="Show atom's indexes on the plot",
            command=self.update_to_show_al_indexes,
            default=self.view_model.to_show_al_indexes,
        )
        self.to_show_al_indexes_checkbox.pack(pady=10, padx=10)

        self.plot_btn: Button = Button(
            self.input_window,
            text="Build the model",
            command=self.plot_al_plane_coordinates,
        )
        self.plot_btn.pack(pady=10, padx=10)

        self.update_tbl_btn: Button = Button(
            self.input_window,
            text="Update the Excel file",
            command=self.update_al_plane_coordinates_file,
        )
        self.update_tbl_btn.pack(pady=10, padx=10)

        self.generate_tbl_btn: Button = Button(
            self.input_window,
            text="Generate the Excel file with Al coordinates for plane",
            command=self.generate_al_plane_coordinates_file,
        )
        self.generate_tbl_btn.pack(pady=10, padx=10)

    def plot_al_plane_coordinates(self) -> None:
        self.view_model.plot_al_plane_coordinates(self.structure_folder)

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

    # def update_file_name(self, value: str) -> None:
    #     self.view_model.set_file_name(value)

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

    def _refresh_file_name_lists(self) -> None:
        path: Path = self.view_model.data_dir / self.structure_folder
        self.file_names: list[str] = FileReader.read_list_of_files(path, format=".xlsx") or ["None"]
        if self.view_model.file_name == "None" and self.file_names:
            self.view_model.set_file_name(self.file_names[0])


class TranslateAlToOtherPlanesWindow(_IntercalationAndSorptionUtils):
    def __init__(self, view_model: VMIntercalationAndSorption, structure_folder: str) -> None:
        super().__init__(view_model, structure_folder)

        self.create_window()

    def create_window(self) -> None:
        self.input_window = ctk.CTkToplevel()
        title: str = f"Translate Al to other planes ({self.structure_folder})"
        self.input_window.title(title)

        description: str = (
            f"If the file {Constants.filenames.AL_FULL_CHANNEL_COORDINATES_XLSX_FILE} exists - it just plots the structure.\n"
            f"If the file above is not found, the program will build the full channel and translate Al atoms to other planes."
        )
        description_label: ctk.CTkLabel = ctk.CTkLabel(
            self.input_window, text=description, wraplength=400
        )
        description_label.pack(pady=10, padx=10)

        # self.input_window.geometry("450x600")
        self.input_window.pack_propagate(True)
        self.input_window.grid_propagate(True)

        self.file_names_dropdown: DropdownList = DropdownList(
            self.input_window,
            options=self.file_names,
            command=self.view_model.set_file_name,
            title="Al coordinates table to translate",
        )
        self.file_names_dropdown.pack(pady=10, padx=10)

        self.try_to_reflect_al_atoms_checkbox: CheckBox = CheckBox(
            self.input_window,
            text="Try to reflect Al atoms\n(if no init file and the Al atoms will be calculated)",
            command=self.update_to_try_to_reflect_al_atoms,
            # title="Try to reflect aluminum atoms",
        )
        self.try_to_reflect_al_atoms_checkbox.pack(pady=10, padx=10)

        self.number_of_planes_input_field: InputField = InputField(
            self.input_window,
            text="Number of planes",
            command=self.update_number_of_planes,
            # title="Number of planes",
            default_value=self.view_model.number_of_planes,
        )
        self.number_of_planes_input_field.pack(pady=10, padx=10)

        self.num_of_min_distances_input_field: InputField = InputField(
            self.input_window,
            text="Number of min distances for bonds",
            command=self.update_num_of_min_distances,
            # title="Number of min distances for bonds",
            default_value=self.view_model.bonds_num_of_min_distances,
        )
        self.num_of_min_distances_input_field.pack(pady=10, padx=10)

        self.bonds_skip_first_distances_input_field = InputField(
            self.input_window, text="Skip first distances for bonds",
            command=self.update_bonds_skip_first_distances,
            default_value=self.view_model.bonds_skip_first_distances,
        )
        self.bonds_skip_first_distances_input_field.pack(pady=10, padx=10)

        self.num_of_al_layers_input_field: InputField = InputField(
            self.input_window,
            text="Number of Al layers on the plot",
            command=self.update_num_of_al_layers,
            # title="Number of AL layers",
            default_value=self.view_model.num_of_al_layers,
        )
        self.num_of_al_layers_input_field.pack(pady=10, padx=10)

        # Checkbox for to_show_indexes
        self.to_show_al_indexes_checkbox = CheckBox(
            self.input_window, text="Show atom's indexes on the plot",
            command=self.update_to_show_al_indexes,
            default=self.view_model.to_show_al_indexes,
        )
        self.to_show_al_indexes_checkbox.pack(pady=10, padx=10)

        self.translate_btn: Button = Button(
            self.input_window,
            text="Translate Al to other planes",
            command=self.translate_al_to_other_planes,
        )
        self.translate_btn.pack(pady=10, padx=10)

        self.update_tbl_btn: Button = Button(
            self.input_window,
            text="Update the Excel file",
            command=self.update_file,
        )
        self.update_tbl_btn.pack(pady=10, padx=10)

    def translate_al_to_other_planes(self) -> None:
        self.view_model.translate_al_to_other_planes(self.structure_folder)

    def update_file(self) -> None:
        try:
            path_to_file: Path = self.view_model.update_al_channel_coordinates(self.structure_folder)
            messagebox.showinfo("Success", f"Al coordinates table saved to {path_to_file}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    # def update_file_name(self, value: str) -> None:
    #     self.view_model.set_file_name(value)

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


class GetAlInChannelDetailsWindow(_IntercalationAndSorptionUtils):
    def __init__(self, view_model: VMIntercalationAndSorption, structure_folder: str) -> None:
        super().__init__(view_model, structure_folder)

        self.create_window()

    def create_window(self) -> None:
        self.input_window = ctk.CTkToplevel()
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

        # self.input_window.geometry("450x300")
        self.input_window.pack_propagate(True)
        self.input_window.grid_propagate(True)

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
        self.get_details_btn.pack(pady=10, padx=10)

    # def update_file_name(self, value: str) -> None:
    #     self.view_model.set_file_name(value)

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
        # self.input_window.geometry("450x600")
        self.input_window.pack_propagate(True)
        self.input_window.grid_propagate(True)

        self.file_names_dropdown: DropdownList = DropdownList(
            self.input_window,
            options=self.file_names,
            command=self.view_model.set_file_name,
            title="Al coordinates table to translate",
        )
        self.file_names_dropdown.pack(pady=10, padx=10)

        self.try_to_reflect_al_atoms_checkbox: CheckBox = CheckBox(
            self.input_window,
            text="Try to reflect Al atoms\n(if no init file and the Al atoms will be calculated)",
            command=self.update_to_try_to_reflect_al_atoms,
            # title="Try to reflect aluminum atoms",
        )
        self.try_to_reflect_al_atoms_checkbox.pack(pady=10, padx=10)

        # self.number_of_planes_input_field: InputField = InputField(
        #     self.input_window,
        #     text="Number of planes (if no init file and the Al atoms will be calculated)",
        #     command=self.update_number_of_planes,
        #     # title="Number of planes",
        #     default_value=self.view_model.number_of_planes,
        # )
        # self.number_of_planes_input_field.pack(pady=10, padx=10)

        self.num_of_min_distances_input_field: InputField = InputField(
            self.input_window,
            text="Number of min distances for bonds",
            command=self.update_num_of_min_distances,
            # title="Number of min distances for bonds",
            default_value=self.view_model.bonds_num_of_min_distances,
        )
        self.num_of_min_distances_input_field.pack(pady=10, padx=10)

        # Input field for bonds_skip_first_distances
        self.bonds_skip_first_distances_input_field = InputField(
            self.input_window, text="Skip first distances for bonds",
            command=self.update_bonds_skip_first_distances,
            default_value=self.view_model.bonds_skip_first_distances,
        )
        self.bonds_skip_first_distances_input_field.pack(pady=10, padx=10)

        self.num_of_al_layers_input_field: InputField = InputField(
            self.input_window,
            text="Number of Al layers on the plot",
            command=self.update_num_of_al_layers,
            # title="Number of Al layers",
            default_value=self.view_model.num_of_al_layers,
        )
        self.num_of_al_layers_input_field.pack(pady=10, padx=10)

        # Checkbox for to_show_indexes
        self.to_show_al_indexes_checkbox = CheckBox(
            self.input_window, text="Show atom's indexes on the plot",
            command=self.update_to_show_al_indexes,
            default=self.view_model.to_show_al_indexes,
        )
        self.to_show_al_indexes_checkbox.pack(pady=10, padx=10)

        self.translate_btn: Button = Button(
            self.input_window,
            text="Translate Al to all channels\nand generate output files",
            command=self.translate_al_to_all_channels,
        )
        self.translate_btn.pack(pady=10, padx=10)

    def translate_al_to_all_channels(self) -> None:
        try:
            paths: tuple[Path, Path, Path] = self.view_model.translate_al_to_all_channels(self.structure_folder)
            messagebox.showinfo("Success", f"Generated files:\n{paths[0].name}\n{paths[1].name}\n{paths[2].name}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def update_file_name(self, value: str) -> None:
        self.view_model.set_file_name(value)

    # def update_number_of_planes(self) -> None:
    #     value = int(self.number_of_planes_input_field.get())
    #     self.view_model.set_number_of_planes(value)

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
