import json
from pathlib import Path

from ..logger import Logger
from .path_builder import PathBuilder

logger = Logger(__name__)


class FileReader:
    @staticmethod
    def read_json_file(
            structure_folder: str, folder_path: Path | str | None = None, file_name: str = "structure_settings.json"
    ) -> dict | list | None:
        """
        Read JSON file (by default 'structure_settings.json').
        If folder_path=None -- uses path to 'result_data' folder.
        """

        if folder_path is None:
            path_to_file: Path = PathBuilder.build_path_to_result_data_file(
                structure_folder=structure_folder, file=file_name)
        else:
            path_to_file: Path = Path(structure_folder) / file_name

        if not path_to_file.exists():
            logger.warning(f"File {path_to_file} not exists.")
            return None

        data_json: str = Path(path_to_file).read_text(encoding="utf-8")
        return json.loads(data_json)
