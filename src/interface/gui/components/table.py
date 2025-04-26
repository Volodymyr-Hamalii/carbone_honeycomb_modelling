import customtkinter as ctk
import pandas as pd
import tkinter as tk


class Table(ctk.CTkFrame):
    def __init__(
            self,
            data: pd.DataFrame,
            master: ctk.CTk | ctk.CTkToplevel,
            title: str = "",
            **kwargs,
    ) -> None:
        super().__init__(master, **kwargs)

        # Set the background color for the entire table
        self.configure(bg_color="white")

        # Create a label for the title using grid
        title_label = ctk.CTkLabel(self, text=title, bg_color="white")
        title_label.grid(row=0, column=0, columnspan=len(data.columns), sticky="ew")

        # Create a canvas to hold the table and scrollbars
        canvas = ctk.CTkCanvas(self, bg="white")
        canvas.grid(row=1, column=0, sticky="nsew")

        # Add scrollbars
        v_scrollbar = ctk.CTkScrollbar(self, orientation='vertical', command=canvas.yview)
        v_scrollbar.grid(row=1, column=1, sticky='ns')
        h_scrollbar = ctk.CTkScrollbar(self, orientation='horizontal', command=canvas.xview)
        h_scrollbar.grid(row=2, column=0, sticky='ew')

        canvas.configure(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)

        # Create a frame inside the canvas
        table_frame = ctk.CTkFrame(canvas, bg_color="white")
        canvas.create_window((0, 0), window=table_frame, anchor='nw')

        # Calculate column widths based on content
        col_widths = []
        for col in data.columns:
            if isinstance(col, tuple):
                # Calculate the maximum width for each level of the MultiIndex
                max_content_width: int = max(
                    max(data[col].astype(str).apply(len).max() for col in data.columns),
                    max(len(str(subcol)) for subcol in col)
                )
            else:
                # Calculate the maximum width for a single-level column
                max_content_width = max(data[col].astype(str).apply(len).max(), len(str(col)))

            col_widths.append(max_content_width + 2)  # Add small padding

        # Handle MultiIndex columns
        if isinstance(data.columns, pd.MultiIndex):
            # Create headers for each level of the MultiIndex
            for level in range(data.columns.nlevels):
                col_start: int = 0
                while col_start < len(data.columns):
                    col_end: int = col_start
                    # Find the range of columns with the same top-level header
                    while (col_end + 1 < len(data.columns) and
                           data.columns[col_start][0] == data.columns[col_end + 1][0]):
                        col_end += 1

                    # Create a header for the top-level
                    if level == 0:
                        header_frame = ctk.CTkFrame(table_frame, bg_color="black")
                        header_frame.grid(
                            row=level,
                            column=col_start + 1,
                            columnspan=(col_end - col_start + 1),
                            sticky="nsew",
                            padx=1,
                            pady=1,
                        )
                        header = tk.Text(
                            header_frame,
                            height=1,
                            width=sum(col_widths[col_start:col_end + 1]),
                            font=("Arial", 10, "bold"),  # Bold font for headers
                            bg="lightgray",  # Background color for headers
                            bd=0,  # No border
                            highlightthickness=0,  # No highlight border
                            wrap="none"  # No text wrapping
                        )
                        header.insert("1.0", data.columns[col_start][0])
                        header.tag_configure("center", justify='center')
                        header.tag_add("center", "1.0", "end")
                        header.config(state="disabled")  # Make the text read-only
                        header.pack(fill="both", expand=True)

                    # Create headers for the second level
                    for col in range(col_start, col_end + 1):
                        header_frame = ctk.CTkFrame(table_frame, bg_color="black")
                        header_frame.grid(row=level + 1, column=col, sticky="nsew", padx=1, pady=1)
                        header = tk.Text(
                            header_frame,
                            height=1,
                            width=col_widths[col],
                            font=("Arial", 10, "bold"),  # Bold font for headers
                            bg="lightgray",  # Background color for headers
                            bd=0,  # No border
                            highlightthickness=0,  # No highlight border
                            wrap="none"  # No text wrapping
                        )
                        header.insert("1.0", data.columns[col][1])
                        header.tag_configure("center", justify='center')
                        header.tag_add("center", "1.0", "end")
                        header.config(state="disabled")  # Make the text read-only
                        header.pack(fill="both", expand=True)

                    col_start: int = col_end + 1

        else:
            # Create a single row of headers
            for i, column in enumerate(data.columns):
                header_frame = ctk.CTkFrame(table_frame, bg_color="black")
                header_frame.grid(row=0, column=i, sticky="nsew", padx=1, pady=1)
                header = tk.Text(
                    header_frame,
                    height=1,
                    width=col_widths[i],
                    font=("Arial", 10, "bold"),  # Bold font for headers
                    bg="lightgray",  # Background color for headers
                    bd=0,  # No border
                    highlightthickness=0,  # No highlight border
                    wrap="none"  # No text wrapping
                )
                header.insert("1.0", column)
                header.tag_configure("center", justify='center')
                header.tag_add("center", "1.0", "end")
                header.config(state="disabled")  # Make the text read-only
                header.pack(fill="both", expand=True)

        # Add index header
        index_header_frame = ctk.CTkFrame(table_frame, bg_color="black")
        index_header_frame.grid(row=0, column=0, sticky="nsew", padx=1, pady=1)
        index_header = tk.Text(
            index_header_frame,
            height=1,
            width=index_width,
            font=("Arial", 10, "bold"),  # Bold font for headers
            bg="lightgray",  # Background color for headers
            bd=0,  # No border
            highlightthickness=0,  # No highlight border
            wrap="none"  # No text wrapping
        )
        index_header.insert("1.0", "Index")
        index_header.tag_configure("center", justify='center')
        index_header.tag_add("center", "1.0", "end")
        index_header.config(state="disabled")  # Make the text read-only
        index_header.pack(fill="both", expand=True)

        # Create the table cells
        for i, (index, row) in enumerate(data.iterrows()):
            # Add index cell
            index_frame = ctk.CTkFrame(table_frame, bg_color="black")
            index_frame.grid(row=i + data.columns.nlevels, column=0, sticky="nsew", padx=1, pady=1)
            index_cell = tk.Text(
                index_frame,
                height=1,
                width=index_width,
                font=("Arial", 9),  # Reduced font size for cells
                bg="lightgray",  # Background color for index cells
                bd=0,  # No border
                highlightthickness=0,  # No highlight border
                wrap="none"  # No text wrapping
            )
            index_cell.insert("1.0", str(index))
            index_cell.tag_configure("center", justify='center')
            index_cell.tag_add("center", "1.0", "end")
            index_cell.config(state="disabled")  # Make the text read-only
            index_cell.pack(fill="both", expand=True)

            for j, value in enumerate(row):
                cell_frame = ctk.CTkFrame(table_frame, bg_color="black")
                cell_frame.grid(
                    row=i + data.columns.nlevels,
                    column=j + 1,
                    sticky="nsew",
                    padx=1,
                    pady=1,
                )
                cell = tk.Text(
                    cell_frame,
                    height=1,
                    width=col_widths[j],
                    font=("Arial", 9),  # Reduced font size for cells
                    bg="white",  # Background color for cells
                    bd=0,  # No border
                    highlightthickness=0,  # No highlight border
                    wrap="none"  # No text wrapping
                )
                cell.insert("1.0", str(value))
                cell.tag_configure("center", justify='center')
                cell.tag_add("center", "1.0", "end")
                cell.config(state="disabled")  # Make the text read-only
                cell.pack(fill="both", expand=True)

        # Make the table adaptive
        self.grid_rowconfigure(1, weight=1)
        self.grid_columnconfigure(0, weight=1)

        # Update scroll region
        table_frame.update_idletasks()
        canvas.config(scrollregion=canvas.bbox("all"))
