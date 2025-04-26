import customtkinter as ctk
import pandas as pd


class Table(ctk.CTkFrame):
    def __init__(
            self,
            data: pd.DataFrame,
            title: str = "",
            **kwargs,
    ) -> None:
        super().__init__(**kwargs)

        # Create a label for the title
        title_label = ctk.CTkLabel(self, text=title)
        title_label.pack()

        # Create a table using the DataFrame
        for i, column in enumerate(data.columns):
            header = ctk.CTkLabel(self, text=column)
            header.grid(row=0, column=i)

        for i, row in data.iterrows():
            for j, value in enumerate(row):
                cell = ctk.CTkLabel(self, text=value)
                cell.grid(row=i+1, column=j)  # type: ignore
