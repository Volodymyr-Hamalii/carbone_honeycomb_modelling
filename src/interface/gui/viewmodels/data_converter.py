from pathlib import Path
import numpy as np
import pandas as pd

from src.utils import Constants, Logger, Inputs, FileReader, FileWriter, DataConverter
from .view_model_params_setter import VMParamsSetter

logger = Logger("Actions")


class VMDataConverter(VMParamsSetter):
    @classmethod
    def convert_file(
            cls,
            init_file_path: Path,
            target_format: str,
            is_init_data_dir: bool = False,
    ) -> Path:
        """
        Convert a file from one format to another. Return the path to the converted file.
        """
        # Determine the source format
        source_format: str = init_file_path.suffix[1:]  # Remove the dot from the suffix

        # Read the data based on the source format
        if source_format == "xlsx":
            df: pd.DataFrame | None = FileReader.read_excel_file(
                structure_folder=init_file_path.parent.name,
                file_name=init_file_path.name,
                # folder_path=init_file_path.parent,
                is_init_data_dir=is_init_data_dir,
            )
            if df is None:
                raise ValueError(f"Failed to read Excel file: {init_file_path}")

            # If more than 3 columns, find columns with "X", "Y", "Z"
            if len(df.columns) > 3:
                x_col: str | None = next((col for col in df.columns if "x" in col.lower()), None)
                y_col: str | None = next((col for col in df.columns if "y" in col.lower()), None)
                z_col: str | None = next((col for col in df.columns if "z" in col.lower()), None)
                if x_col and y_col and z_col:
                    df = df[[x_col, y_col, z_col]]
                else:
                    raise ValueError("Could not find X, Y, Z columns in Excel file.")

        elif source_format == "dat":
            data: np.ndarray = FileReader.read_dat_file(
                structure_folder=init_file_path.parent.name,
                file_name=init_file_path.name,
                # folder_path=init_file_path.parent,
                is_init_data_dir=is_init_data_dir,
            )
            df = pd.DataFrame(data, columns=["X", "Y", "Z"])

        elif source_format == "pdb":
            # Assuming a method to read PDB files into a DataFrame
            data: np.ndarray = FileReader.read_pdb_file(
                structure_folder=init_file_path.parent.name,
                # folder_path=init_file_path.parent,
                file_name=init_file_path.name,
                is_init_data_dir=False
            )
            df = pd.DataFrame(data, columns=["X", "Y", "Z"])

        else:
            raise ValueError(f"Unsupported source format: {source_format}")

        # Write the data based on the target format
        if target_format == "xlsx":
            FileWriter.write_excel_file(
                df=df,
                structure_folder=init_file_path.parent.name,
                file_name=init_file_path.stem + ".xlsx",
                sheet_name="Sheet1",
                # folder_path=init_file_path.parent,
                is_init_data_dir=is_init_data_dir,
            )

        elif target_format == "dat":
            dat_lines: list[str] = DataConverter.convert_df_to_dat(df)
            FileWriter.write_dat_file(
                data_lines=dat_lines,
                path_to_file=init_file_path.with_suffix(".dat"),
                structure_folder=init_file_path.parent.name
            )

        elif target_format == "pdb":
            pdb_lines: list[str] = DataConverter.convert_df_to_pdb(df)
            FileWriter.write_pdb_file(
                data_lines=pdb_lines,
                path_to_file=init_file_path.with_suffix(".pdb"),
                structure_folder=init_file_path.parent.name
            )

        else:
            raise ValueError(f"Unsupported target format: {target_format}")

        return init_file_path.with_suffix(f".{target_format}")

    # @staticmethod
    # def convert_init_dat_to_pdb(structure_folder: str) -> None:
    #     """
    #     Convert init_data/{structure_folder}/ljout.dat into result_data/{structure_folder}/ljout-from-init-dat.pdb
    #     Also create result_data/{structure_folder}/structure_settings.json template if it didn't exist.
    #     """

    #     FileConverter.dat_to_pdb(structure_folder=structure_folder)

    # @staticmethod
    # def convert_excel_to_dat(structure_folder: str) -> None:
    #     """
    #     Convert result_data/{structure_folder}/{file_name}.xlsx into result_data/{structure_folder}/{file_name}.dat.
    #     """

    #     file_name: str = Inputs.text_input(
    #         to_set,
    #         default_value=Constants.filenames.AL_FULL_CHANNEL_COORDINATES_XLSX_FILE,
    #         text="File name to convert",
    #         env_id="file_name_to_convert")

    #     folder_path: Path = Constants.path.RESULT_DATA_PATH

    #     df: pd.DataFrame | None = FileReader.read_excel_file(
    #         structure_folder=structure_folder,
    #         file_name=file_name,
    #         folder_path=folder_path,
    #     )

    #     if df is None:
    #         logger.error(f"File {file_name} not found at {folder_path}")
    #         return

    #     file_name_dat: str = file_name.replace(".xlsx", ".dat")

    #     FileWriter.write_dat_file(
    #         data_lines=df.to_numpy(),
    #         structure_folder=structure_folder,
    #         path_to_file=folder_path / structure_folder / file_name_dat,
    #     )
