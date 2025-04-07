import customtkinter as ctk
from typing import Callable


class CheckBox(ctk.CTkCheckBox):
    def __init__(self, master, text: str, command: Callable, default: bool = False, **kwargs) -> None:
        self.var = ctk.BooleanVar(value=default)
        super().__init__(master, text=text, command=command, variable=self.var, **kwargs)
