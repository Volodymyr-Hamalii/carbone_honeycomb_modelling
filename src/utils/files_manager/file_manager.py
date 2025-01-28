from pathlib import Path
from .path_builder import PathBuilder


class FileManager:
    @staticmethod
    def _get_path_to_file(
            structure_folder: str,
            file_name: str,
            folder_path: Path | str | None = None,
            is_init_data_dir: bool = True,
    ) -> Path:
        if folder_path is None:
            if is_init_data_dir:
                return PathBuilder.build_path_to_init_data_file(
                    structure_folder=structure_folder, file=file_name)
            else:
                return PathBuilder.build_path_to_result_data_file(
                    structure_folder=structure_folder, file=file_name)
        else:
            return Path(folder_path) / structure_folder / file_name
