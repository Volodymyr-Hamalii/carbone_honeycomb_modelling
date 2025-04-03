from pathlib import Path
import pandas as pd

from src.utils import Constants, Logger, Inputs, FileReader, FileConverter, FileWriter


logger = Logger("Actions")


class AppActionsInitDataParsing:
    @staticmethod
    def convert_init_dat_to_pdb(structure_folder: str, to_set: bool) -> None:
        """
        Convert init_data/{structure_folder}/ljout.dat into result_data/{structure_folder}/ljout-from-init-dat.pdb
        Also create result_data/{structure_folder}/structure_settings.json template if it didn't exist.
        """

        FileConverter.dat_to_pdb(structure_folder=structure_folder)

    @staticmethod
    def convert_excel_to_dat(structure_folder: str, to_set: bool) -> None:
        """
        Convert result_data/{structure_folder}/{file_name}.xlsx into result_data/{structure_folder}/{file_name}.dat.
        """

        file_name: str = Inputs.text_input(
            to_set,
            default_value=Constants.filenames.AL_FULL_CHANNEL_COORDINATES_XLSX_FILE,
            text="File name to convert",
            env_id="file_name_to_convert")

        folder_path: Path = Constants.path.RESULT_DATA_PATH

        df: pd.DataFrame | None = FileReader.read_excel_file(
            structure_folder=structure_folder,
            file_name=file_name,
            folder_path=folder_path,
        )

        if df is None:
            logger.error(f"File {file_name} not found at {folder_path}")
            return

        file_name_dat: str = file_name.replace(".xlsx", ".dat")

        FileWriter.write_dat_file(
            data_lines=df.to_numpy(),
            structure_folder=structure_folder,
            path_to_file=folder_path / structure_folder / file_name_dat,
        )
