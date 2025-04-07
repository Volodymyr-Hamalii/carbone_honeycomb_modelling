import customtkinter as ctk
from typing import Callable


class Button(ctk.CTkButton):
    def __init__(self, master, text: str, command: Callable, state: str = "active", **kwargs) -> None:
        super().__init__(master, text=text, command=command, **kwargs)
        self.configure(state=state)
