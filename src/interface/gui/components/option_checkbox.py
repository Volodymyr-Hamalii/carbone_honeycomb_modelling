import customtkinter as ctk
from typing import Callable, List


class CheckBox(ctk.CTkCheckBox):
    def __init__(self, master, text: str, command: Callable, **kwargs) -> None:
        super().__init__(master, text=text, command=command, **kwargs)


class DropdownList(ctk.CTkOptionMenu):
    def __init__(self, master, options: List[str], command: Callable, **kwargs) -> None:
        super().__init__(master, values=options, command=command, **kwargs)
