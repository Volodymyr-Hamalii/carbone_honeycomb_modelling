import tkinter as tk
from tkinter import ttk


class AppGUI:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("App Actions")

        # Create a frame for inputs
        self.frame = ttk.Frame(root, padding="10")
        self.frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Boolean input example
        self.to_set_var = tk.BooleanVar(value=False)
        self.to_set_checkbox = ttk.Checkbutton(self.frame, text="To Set", variable=self.to_set_var)
        self.to_set_checkbox.grid(row=0, column=0, sticky=tk.W)

        # Text input example
        self.text_var = tk.StringVar(value="default")
        self.text_entry = ttk.Entry(self.frame, textvariable=self.text_var)
        self.text_entry.grid(row=1, column=0, sticky=(tk.W, tk.E))

        # Dropdown example
        self.dropdown_var = tk.StringVar(value="Option1")
        self.dropdown = ttk.Combobox(self.frame, textvariable=self.dropdown_var)
        self.dropdown['values'] = ("Option1", "Option2", "Option3")
        self.dropdown.grid(row=2, column=0, sticky=(tk.W, tk.E))

        # Run button
        self.run_button = ttk.Button(self.frame, text="Run", command=self.run_action)
        self.run_button.grid(row=3, column=0, sticky=tk.W)

    def run_action(self):
        # Retrieve values from GUI elements
        to_set = self.to_set_var.get()
        text_input = self.text_var.get()
        dropdown_selection = self.dropdown_var.get()

        # Call your existing logic here
        print(f"To Set: {to_set}, Text Input: {text_input}, Dropdown: {dropdown_selection}")


if __name__ == "__main__":
    root = tk.Tk()
    app = AppGUI(root)
    root.mainloop()
