import customtkinter as ctk
from typing import Callable


class InputField(ctk.CTkEntry):
    def __init__(
            self,
            master,
            text: str,
            command: Callable,
            state: str = "normal",
            default_value: str | int | float | None = None,
            **kwargs,
    ) -> None:
        # Create a frame to hold the entry and button
        self.frame = ctk.CTkFrame(master)
        self.frame.pack(fill="x", padx=5, pady=5)

        # Initialize the CTkEntry within the frame
        super().__init__(self.frame, **kwargs)
        self.configure(state=state)

        self.command = command

        if default_value is not None:
            self.insert(0, default_value)

        self.bind("<Return>", lambda event: self.command())

        # Create and pack the label above the frame
        self.label = ctk.CTkLabel(master, text=text)
        self.label.pack(in_=self.frame, side="top", fill="x")

        # Pack the entry in the frame
        self.pack(side="left", fill="x", expand=True)

        # Create and pack the "Apply" button in the frame
        self.apply_button = ctk.CTkButton(self.frame, text="Apply", command=self.command)
        self.apply_button.pack(side="right")
