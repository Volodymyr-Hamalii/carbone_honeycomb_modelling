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
        self.input_window.geometry("450x600")

        self.file_names_dropdown: DropdownList = DropdownList(
            self.input_window,
            options=self.file_names,
            command=self.update_file_name,
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
            text="Number of min distances",
            command=self.update_num_of_min_distances,
            # title="Number of min distances for bonds",
            default_value=self.view_model.bonds_num_of_min_distances,
        )
        self.num_of_min_distances_input_field.pack(pady=10, padx=10)

        self.num_of_al_layers_input_field: InputField = InputField(
            self.input_window,
            text="Number of AL layers",
            command=self.update_num_of_al_layers,
            # title="Number of AL layers",
            default_value=self.view_model.num_of_al_layers,
        )
        self.num_of_al_layers_input_field.pack(pady=10, padx=10)

        self.next_btn: Button = Button(
            self.input_window,
            text="Build the model and update the table",
            command=self.update_al_coordinates_table,
        )
        self.next_btn.pack(pady=10, padx=10)

        self.update_tbl_btn: Button = Button(
            self.input_window,
            text="Update the Excel file",
            command=self.update_excel_file,
        )
        self.update_tbl_btn.pack(pady=10, padx=10)

    def update_al_coordinates_table(self) -> None:
        self.view_model.update_al_coordinates_tbl(self.structure_folder)

    def update_excel_file(self) -> None:
        self.view_model.update_plane_tbl_excel_file(self.structure_folder)

    def update_file_name(self) -> None:
        value: str = self.file_names_dropdown.get()
        self.view_model.set_file_name(value)

    def update_number_of_planes(self) -> None:
        value = int(self.number_of_planes_input_field.get())
        self.view_model.set_number_of_planes(value)

    def update_num_of_min_distances(self) -> None:
        value = int(self.num_of_min_distances_input_field.get())
        self.view_model.set_bonds_num_of_min_distances(value)

    def update_num_of_al_layers(self) -> None:
        value = int(self.num_of_al_layers_input_field.get())
        self.view_model.set_num_of_al_layers(value)


class TranslateAlToOtherPlanesWindow(_IntercalationAndSorptionUtils):
    def __init__(self, view_model: VMIntercalationAndSorption, structure_folder: str) -> None:
        super().__init__(view_model, structure_folder)

        self.create_window()

    def create_window(self) -> None:
        self.input_window = ctk.CTkToplevel()
        title: str = f"Translate Al to other planes ({self.structure_folder})"
        self.input_window.title(title)
        self.input_window.geometry("600x400")

        self.file_names_dropdown: DropdownList = DropdownList(
            self.input_window,
            options=self.file_names,
            command=self.update_file_name,
            title="Al coordinates table to translate",
        )
        self.file_names_dropdown.pack(pady=10, padx=10)

        self.try_to_reflect_al_atoms_checkbox: CheckBox = CheckBox(
            self.input_window,
            text="Try to reflect Al atoms (if no init file and the Al atoms will be calculated)",
            command=self.update_to_try_to_reflect_al_atoms,
            # title="Try to reflect aluminum atoms",
        )
        self.try_to_reflect_al_atoms_checkbox.pack(pady=10, padx=10)

        self.number_of_planes_input_field: InputField = InputField(
            self.input_window,
            text="Number of planes",
            command=self.update_number_of_planes,
            # title="Number of planes",
        )
        self.number_of_planes_input_field.pack(pady=10, padx=10)

        self.num_of_min_distances_input_field: InputField = InputField(
            self.input_window,
            text="Number of min distances",
            command=self.update_num_of_min_distances,
            # title="Number of min distances for bonds",
        )
        self.num_of_min_distances_input_field.pack(pady=10, padx=10)

        self.num_of_al_layers_input_field: InputField = InputField(
            self.input_window,
            text="Number of AL layers",
            command=self.update_num_of_al_layers,
            # title="Number of AL layers",
        )
        self.num_of_al_layers_input_field.pack(pady=10, padx=10)

        self.next_btn: Button = Button(
            self.input_window,
            text="Translate Al to other planes",
            command=self.translate_al_to_other_planes,
        )
        self.next_btn.pack(pady=10, padx=10)

    def translate_al_to_other_planes(self) -> None:
        try:
            path_to_file: Path = self.view_model.translate_al_to_other_planes(self.structure_folder)
            messagebox.showinfo("Success", f"Al coordinates table saved to {path_to_file}")
        except ValueError as e:
            messagebox.showerror("Error", str(e))

    def update_file_name(self) -> None:
        value: str = self.file_names_dropdown.get()
        self.view_model.set_file_name(value)

    def update_number_of_planes(self) -> None:
        value = int(self.number_of_planes_input_field.get())
        self.view_model.set_number_of_planes(value)

    def update_num_of_min_distances(self) -> None:
        value = int(self.num_of_min_distances_input_field.get())
        self.view_model.set_bonds_num_of_min_distances(value)

    def update_to_try_to_reflect_al_atoms(self) -> None:
        value = bool(self.try_to_reflect_al_atoms_checkbox.get())
        self.view_model.set_to_try_to_reflect_al_atoms(value)

    def update_num_of_al_layers(self) -> None:
        value = int(self.num_of_al_layers_input_field.get())
        self.view_model.set_num_of_al_layers(value)


class TranslateAlToAllChannelsWindow(_IntercalationAndSorptionUtils):
    def __init__(self, view_model: VMIntercalationAndSorption, structure_folder: str) -> None:
        super().__init__(view_model, structure_folder)

        self.create_window()

    def create_window(self) -> None:
        self.input_window = ctk.CTkToplevel()
        title: str = f"Translate Al to all channels ({self.structure_folder})"
        self.input_window.title(title)
        self.input_window.geometry("600x400")

        self.file_names_dropdown: DropdownList = DropdownList(
            self.input_window,
            options=self.file_names,
            command=self.update_file_name,
            title="Al coordinates table to translate",
        )
        self.file_names_dropdown.pack(pady=10, padx=10)

    def update_file_name(self) -> None:
        value: str = self.file_names_dropdown.get()
        self.view_model.set_file_name(value)

    # def update_number_of_planes(self) -> None:
    #     value = int(self.number_of_planes_input_field.get())
    #     self.view_model.set_number_of_planes(value)

    # def update_num_of_min_distances(self) -> None:
    #     value = int(self.num_of_min_distances_input_field.get())
    #     self.view_model.set_bonds_num_of_min_distances(value)

    # def update_to_try_to_reflect_al_atoms(self) -> None:
    #     value = bool(self.try_to_reflect_al_atoms_checkbox.get())
    #     self.view_model.set_to_try_to_reflect_al_atoms(value)

    # def update_num_of_al_layers(self) -> None:
    #     value = int(self.num_of_al_layers_input_field.get())
    #     self.view_model.set_num_of_al_layers(value)


class GetAlInChannelDetailsWindow(_IntercalationAndSorptionUtils):
    def __init__(self, view_model: VMIntercalationAndSorption, structure_folder: str) -> None:
        super().__init__(view_model, structure_folder)

        self.create_window()

    def create_window(self) -> None:
        self.input_window = ctk.CTkToplevel()
        title: str = f"Get Al in channel details ({self.structure_folder})"
        self.input_window.title(title)
        self.input_window.geometry("600x400")

        self.file_names_dropdown: DropdownList = DropdownList(
            self.input_window,
            options=self.file_names,
            command=self.update_file_name,
            title="Al coordinates table to analyze",
        )
        self.file_names_dropdown.pack(pady=10, padx=10)

    def update_file_name(self) -> None:
        value: str = self.file_names_dropdown.get()
        self.view_model.set_file_name(value)

    # def update_number_of_planes(self) -> None:
    #     value = int(self.number_of_planes_input_field.get())
    #     self.view_model.set_number_of_planes(value)

    # def update_num_of_min_distances(self) -> None:
    #     value = int(self.num_of_min_distances_input_field.get())
    #     self.view_model.set_bonds_num_of_min_distances(value)

    # def update_to_try_to_reflect_al_atoms(self) -> None:
    #     value = bool(self.try_to_reflect_al_atoms_checkbox.get())
    #     self.view_model.set_to_try_to_reflect_al_atoms(value)

    # def update_num_of_al_layers(self) -> None:
    #     value = int(self.num_of_al_layers_input_field.get())
    #     self.view_model.set_num_of_al_layers(value)
