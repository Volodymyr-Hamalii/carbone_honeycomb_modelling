import customtkinter as ctk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


class MyApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Carbon Honeycomb Viewer")
        self.geometry("600x400")

        # Checkbox
        self.show_bonds = ctk.CTkCheckBox(self, text="Show Bonds")
        self.show_bonds.pack()

        # Input box
        self.param_label = ctk.CTkLabel(self, text="Distance Threshold:")
        self.param_label.pack()
        self.param_input = ctk.CTkEntry(self)
        self.param_input.insert(0, "2.5")
        self.param_input.pack()

        # Button
        self.run_btn = ctk.CTkButton(self, text="Visualize", command=self.visualize)
        self.run_btn.pack(pady=10)

    def visualize(self):
        value = float(self.param_input.get())
        show_bonds = self.show_bonds.get()

        # Placeholder for your visualization logic
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter([0, 1, 2], [0, 1, 0], [0, 1, 2])

        canvas = FigureCanvasTkAgg(fig, master=self)
        canvas.draw()
        canvas.get_tk_widget().pack()


app = MyApp()
app.mainloop()
