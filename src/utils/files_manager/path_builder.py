from pathlib import Path
from ..constants import Constants


class PathBuilder:

    ### INIT DATA ###

    @staticmethod
    def build_path_to_init_data_dir(
            path_to_script: Path | str = Constants.path.path_to_root_script_dir,
            init_dir: Path | str = Constants.filenames.INIT_DATA_DIR) -> Path:

        return Path(path_to_script) / init_dir

    @classmethod
    def build_path_to_init_data_file(
            cls,
            structure_folder: str | None = None,
            path_to_script: Path = Constants.path.path_to_root_script_dir,
            init_dir: Path | str = Constants.filenames.INIT_DATA_DIR,
            file: Path | str = Constants.filenames.INIT_DAT_FILE) -> Path:

        path_dir: Path = cls.build_path_to_init_data_dir(path_to_script, init_dir)

        if structure_folder is not None:
            return path_dir / structure_folder / file

        return path_dir / file

    ### RESULT DATA ###

    @staticmethod
    def build_path_to_result_data_dir(
            path_to_script: Path | str = Constants.path.path_to_root_script_dir,
            result_dir: Path | str = Constants.filenames.RESULT_DATA_DIR) -> Path:

        path_to_result_dir: Path = Path(path_to_script) / result_dir

        # Create directory if it doesn't exist
        path_to_result_dir.mkdir(parents=True, exist_ok=True)

        return path_to_result_dir

    @classmethod
    def build_path_to_result_data_file(
            cls,
            structure_folder: str,
            path_to_script: Path | str = Constants.path.path_to_root_script_dir,
            result_dir: Path | str = Constants.filenames.RESULT_DATA_DIR,
            file: Path | str = Constants.filenames.INIT_PDB_FILE) -> Path:

        path_to_result_dir: Path = cls.build_path_to_result_data_dir(path_to_script, result_dir)

        path_to_data_in_result_dir: Path = path_to_result_dir / structure_folder

        # Create directory if it doesn't exist
        path_to_data_in_result_dir.mkdir(parents=True, exist_ok=True)

        return path_to_data_in_result_dir / file
