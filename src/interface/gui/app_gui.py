import customtkinter as ctk
from src.interface.gui.viewmodels.viewmodels import AppViewModel
from src.interface.gui.components import Button, CheckBox, PlotWindow


class AppGui(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Carbon Honeycomb Viewer")
        self.geometry("600x400")

        # Initialize ViewModel
        self.view_model = AppViewModel()

        # Create GUI components
        self.create_widgets()

    def create_widgets(self):
        # Button to show init structure
        self.show_init_structure_btn = Button(
            self, text="Show Init Structure", command=self.open_show_init_structure_window)

        self.show_init_structure_btn.pack(pady=10)

    def open_show_init_structure_window(self):
        # Create a new window for input parameters
        self.input_window = ctk.CTkToplevel(self)
        self.input_window.title("Input Parameters")
        self.input_window.geometry("300x200")

        # Checkbox for to_build_bonds
        self.to_build_bonds_checkbox = CheckBox(
            self.input_window, text="Build Bonds", command=self.update_to_build_bonds)
        self.to_build_bonds_checkbox.pack(pady=10)

        # Button to proceed to the next step
        self.next_btn = Button(
            self.input_window, text="Next", command=self.run_show_structure)
        self.next_btn.pack(pady=10)

    def update_to_build_bonds(self):
        value = bool(self.to_build_bonds_checkbox.get())
        print(value)
        self.view_model.set_show_bonds(value)

    def run_show_structure(self):
        # Close the input window
        self.input_window.destroy()

        # Run the action and show the plot
        self.view_model.run_action(2.5)
        self.show_plot_window()

    def show_plot_window(self):
        # Create a new window to display the plot
        plot_window = PlotWindow(self, title="Structure Visualization")
        # Here you would integrate the actual plotting logic
        # For example, using matplotlib to draw the plot in the plot_window


app = AppGui()
app.mainloop()
