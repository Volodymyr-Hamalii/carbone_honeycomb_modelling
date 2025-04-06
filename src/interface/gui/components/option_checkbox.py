import customtkinter as ctk
from typing import Callable, List


class CheckBox(ctk.CTkCheckBox):
    def __init__(self, master, text: str, command: Callable, default: bool = False, **kwargs) -> None:
        self.var = ctk.BooleanVar(value=default)
        super().__init__(master, text=text, command=command, variable=self.var, **kwargs)


class DropdownList(ctk.CTkOptionMenu):
    def __init__(
            self,
            master,
            options: List[str],
            command: Callable,
            title: str = "",
            is_disabled: bool = False,
            **kwargs,
    ) -> None:
        kwargs['fg_color'] = 'white'
        kwargs['text_color'] = 'black'

        super().__init__(master, values=options, command=command, **kwargs)

        if title:
            self.title: ctk.CTkLabel = ctk.CTkLabel(master, text=title)
            self.title.pack(pady=(10, 0), padx=10)

        if is_disabled:
            self.configure(state=ctk.DISABLED)
