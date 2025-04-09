import customtkinter as ctk

from ..viewmodels import VMShowInitData
from ..components import Button, CheckBox, PlotWindow


class ChannelDetailsWindow:
    def __init__(self, view_model: VMShowInitData, structure_folder: str) -> None:
        self.structure_folder: str = structure_folder
        self.view_model: VMShowInitData = view_model
        self.create_window()

    def create_window(self) -> None:
        self.input_window = ctk.CTkToplevel()
        self.input_window.pack_propagate(True)
        self.input_window.grid_propagate(True)

        self.input_window.title(f"Show channel parameters ({self.structure_folder})")

        # Checkbox for to_show_channel_angles
        self.to_show_channel_angles_checkbox = CheckBox(
            self.input_window, text="Show channel angles",
            command=self.update_to_show_channel_angles,
            default=self.view_model.to_show_channel_angles,
        )
        self.to_show_channel_angles_checkbox.pack(pady=10, padx=10)

        # Checkbox for to_show_dists_to_plane
        self.to_show_dists_to_plane_checkbox = CheckBox(
            self.input_window, text="Show distances to plane",
            command=self.update_to_show_dists_to_plane,
            default=self.view_model.to_show_dists_to_plane,
        )
        self.to_show_dists_to_plane_checkbox.pack(pady=10, padx=10)

        # Checkbox for to_show_dists_to_edges
        self.to_show_dists_to_edges_checkbox = CheckBox(
            self.input_window, text="Show distances to edges",
            command=self.update_to_show_dists_to_edges,
            default=self.view_model.to_show_dists_to_edges,
        )
        self.to_show_dists_to_edges_checkbox.pack(pady=10, padx=10)

        # Checkbox for to_show_coordinates
        self.to_show_coordinates_checkbox = CheckBox(
            self.input_window, text="Show coordinates",
            command=self.update_to_show_coordinates,
            default=self.view_model.to_show_coordinates,
        )
        self.to_show_coordinates_checkbox.pack(pady=10, padx=10)

        # Button to proceed to the next step
        self.next_btn = Button(
            self.input_window, text="Show",
            command=self.get_channel_details,
        )
        self.next_btn.pack(pady=(10, 25))

    def get_channel_details(self) -> None:
        self.view_model.get_channel_details(self.structure_folder)
        self.show_plot_window()

    def update_to_show_coordinates(self) -> None:
        value = bool(self.to_show_coordinates_checkbox.get())
        self.view_model.set_to_show_coordinates(value)

    def update_to_show_dists_to_plane(self) -> None:
        value = bool(self.to_show_dists_to_plane_checkbox.get())
        self.view_model.set_to_show_dists_to_plane(value)

    def update_to_show_dists_to_edges(self) -> None:
        value = bool(self.to_show_dists_to_edges_checkbox.get())
        self.view_model.set_to_show_dists_to_edges(value)

    def update_to_show_channel_angles(self) -> None:
        value = bool(self.to_show_channel_angles_checkbox.get())
        self.view_model.set_to_show_channel_angles(value)

    def show_plot_window(self) -> None:
        plot_window = PlotWindow(self.input_window, title="Channel Details")
        plot_window.destroy()
