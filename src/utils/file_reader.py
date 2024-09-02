import os
import json
from .path_builder import PathBuilder
from .logger import Logger

logger = Logger(__name__)


class FileReader:
    @staticmethod
    def read_json_file(
            structure_folder: str, folder_path: str | None = None, file_name: str = "structure_settings.json"
    ) -> dict | list | None:
        """
        Read JSON file (by default 'structure_settings.json').
        If folder_path=None -- uses path to 'relust_data' folder.
        """

        if folder_path is None:
            path_to_file: str = PathBuilder.build_path_to_result_data_file(
                structure_folder=structure_folder, file=file_name)
        else:
            path_to_file: str = os.path.join(folder_path, file_name)

        if not os.path.exists(path_to_file):
            logger.warning(f"File {path_to_file} not exists.")
            return None

        with open(path_to_file, "r") as f:
            data_json: str = f.read()
            return json.loads(data_json)
