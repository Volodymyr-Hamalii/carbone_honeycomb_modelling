import customtkinter as ctk
from typing import Callable


class DropdownList(ctk.CTkOptionMenu):
    def __init__(
            self,
            master,
            options: list[str],
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

    def set_options(self, options: list[str], default_value: str = "") -> None:
        self.configure(values=options)
        if default_value:
            self.set(default_value)
