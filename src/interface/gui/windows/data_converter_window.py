from pathlib import Path
import customtkinter as ctk
from tkinter import messagebox

from src.utils import Logger, Constants, FileReader

from ..viewmodels import VMDataConverter
from ..components import Button, CheckBox, DropdownList


logger = Logger("DataOperationsWindow")


class DataConverterWindow:
    def __init__(self, view_model: VMDataConverter, structure_folder: str) -> None:
        self.view_model: VMDataConverter = view_model
        self.structure_folder: str = structure_folder

        # self.excel_file_names: list[str] = ["None"]
        # self.dat_file_names: list[str] = ["None"]
        # self.pdb_file_names: list[str] = ["None"]
        self.file_names: list[str] = ["None"]
        self._refresh_file_name_lists()

        self.create_window()

    def create_window(self) -> None:
        self.input_window = ctk.CTkToplevel()
        title: str = f"Data operations ({self.structure_folder})"
        self.input_window.title(title)

        # self.input_window.pack_propagate(False)
        # self.input_window.grid_propagate(False)
        self.input_window.geometry("600x400")

        # Create dropdown lists
        # self.excel_file_names_dropdown: DropdownList = DropdownList(
        #     self.input_window,
        #     options=self.excel_file_names,
        #     command=self.update_excel_file_name,
        # )
        # self.excel_file_names_dropdown.pack(pady=10, padx=10)

        # self.dat_file_names_dropdown: DropdownList = DropdownList(
        #     self.input_window,
        #     options=self.dat_file_names,
        #     command=self.update_dat_file_name,
        # )
        # self.dat_file_names_dropdown.pack(pady=10, padx=10)

        # self.pdb_file_names_dropdown: DropdownList = DropdownList(
        #     self.input_window,
        #     options=self.pdb_file_names,
        #     command=self.update_pdb_file_name,
        # )
        # self.pdb_file_names_dropdown.pack(pady=10, padx=10)

        self.data_dir_dropdown: DropdownList = DropdownList(
            self.input_window,
            # options=[Constants.filenames.RESULT_DATA_DIR, Constants.filenames.INIT_DATA_DIR],
            options=[Constants.filenames.RESULT_DATA_DIR],
            command=self.update_data_dir,
            title="Data directory",
            is_disabled=True,
        )
        self.data_dir_dropdown.pack(pady=10, padx=10)

        self.file_names_dropdown: DropdownList = DropdownList(
            self.input_window,
            options=self.file_names,
            command=self.update_file_name,
            title="File to convert",
        )
        self.file_names_dropdown.pack(pady=10, padx=10)

        self.formats_dropdown: DropdownList = DropdownList(
            self.input_window,
            options=self.view_model.available_formats,
            command=self.update_formats,
            title="Target file format",
        )
        self.formats_dropdown.pack(pady=10, padx=10)

        # Button to proceed to the next step
        self.next_btn = Button(
            self.input_window, text="Convert", command=self.convert_file
        )
        self.next_btn.pack(pady=10, padx=10)

    # def update_excel_file_name(self, value: str) -> None:
    #     self.view_model.set_excel_file_name(value)

    # def update_dat_file_name(self, value: str) -> None:
    #     self.view_model.set_dat_file_name(value)

    # def update_pdb_file_name(self, value: str) -> None:
    #     self.view_model.set_pdb_file_name(value)

    def update_data_dir(self, value: str) -> None:
        self.view_model.set_data_dir(value)
        self._refresh_file_name_lists()

    def update_file_name(self, value: str) -> None:
        self.view_model.set_file_name(value)

    def update_formats(self, value: str) -> None:
        self.view_model.set_file_format(value)

    def _refresh_file_name_lists(self) -> None:
        path: Path = self.view_model.data_dir / self.structure_folder
        # self.excel_file_names: list[str] = FileReader.read_list_of_files(path, ".xlsx") or ["None"]
        # self.dat_file_names: list[str] = FileReader.read_list_of_files(path, ".dat") or ["None"]
        # self.pdb_file_names: list[str] = FileReader.read_list_of_files(path, ".pdb") or ["None"]
        self.file_names: list[str] = FileReader.read_list_of_files(path) or ["None"]

        if self.view_model.file_name == "None":
            self.view_model.set_file_name(self.file_names[0])

        if self.view_model.file_format == "None":
            self.view_model.set_file_format(self.view_model.available_formats[0])

    def convert_file(self) -> None:
        try:
            if self.view_model.file_name.endswith(self.view_model.file_format):
                messagebox.showinfo("Info", "File already in the correct format.")
                return

            converted_file_path: Path = self.view_model.convert_file(
                init_file_path=self.view_model.data_dir / self.structure_folder / self.view_model.file_name,
                target_format=self.view_model.file_format
            )
            messagebox.showinfo("Success", f"File converted successfully and saved to {converted_file_path}.")

        except Exception as e:
            logger.error(f"Error converting file: {e}")
            messagebox.showerror("Error", str(e))
