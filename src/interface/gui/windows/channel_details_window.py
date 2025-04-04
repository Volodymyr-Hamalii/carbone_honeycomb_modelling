import customtkinter as ctk

from ..viewmodels.show_init_data import ViewModelShowInitData
from ..components import Button, CheckBox, PlotWindow


class ChannelDetailsWindow:
    def __init__(self, view_model: ViewModelShowInitData) -> None:
        self.view_model: ViewModelShowInitData = view_model
        self.create_window()

    def create_window(self) -> None:
        self.input_window = ctk.CTkToplevel()
        self.input_window.title("Show channel parameters")
        self.input_window.pack_propagate(False)
        self.input_window.grid_propagate(False)

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
        self.next_btn.pack(pady=10)

    def get_channel_details(self) -> None:
        self.view_model.get_channel_details()
        self.show_plot_window()

    def update_to_show_coordinates(self) -> None:
        value = bool(self.to_show_coordinates_checkbox.get())
        self.view_model.set_to_show_coordinates(value)

    def show_plot_window(self) -> None:
        plot_window = PlotWindow(self.input_window, title="Channel Details")
        plot_window.destroy()
