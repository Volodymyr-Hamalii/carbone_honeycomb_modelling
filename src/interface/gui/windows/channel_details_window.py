from ..viewmodels import VMShowInitData
from ..components import Button, CheckBox, PlotWindow
from .windows_template import WindowsTemplate


class ChannelDetailsWindow(WindowsTemplate):
    def __init__(self, view_model: VMShowInitData, structure_folder: str) -> None:
        super().__init__()
        self.structure_folder: str = structure_folder
        self.view_model: VMShowInitData = view_model
        self.create_window(
            title=f"Show channel parameters ({self.structure_folder})",
        )
        self.create_ui()

    def create_ui(self) -> None:
        # Checkbox for to_show_channel_angles
        self.to_show_channel_angles_checkbox: CheckBox = self.pack_check_box(
            self.window, text="Show channel angles",
            command=self.update_to_show_channel_angles,
            default=self.view_model.to_show_channel_angles,
        )

        # Checkbox for to_show_dists_to_plane
        self.to_show_dists_to_plane_checkbox: CheckBox = self.pack_check_box(
            self.window, text="Show distances to plane",
            command=self.update_to_show_dists_to_plane,
            default=self.view_model.to_show_dists_to_plane,
        )

        # Checkbox for to_show_dists_to_edges
        self.to_show_dists_to_edges_checkbox: CheckBox = self.pack_check_box(
            self.window, text="Show distances to edges",
            command=self.update_to_show_dists_to_edges,
            default=self.view_model.to_show_dists_to_edges,
        )

        # Checkbox for to_show_coordinates
        self.to_show_coordinates_checkbox: CheckBox = self.pack_check_box(
            self.window, text="Show coordinates",
            command=self.update_to_show_coordinates,
            default=self.view_model.to_show_coordinates,
        )

        # Button to proceed to the next step
        self.btn: Button = self.pack_button(
            self.window, text="Show",
            command=self.get_channel_details,
            pady=(10, 25),
        )

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
        plot_window: PlotWindow = PlotWindow(self.window, title="Channel Details")
        plot_window.destroy()
