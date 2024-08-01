import os
from .constants import Constants


class PathBuilder:

    ### INIT DATA ###

    @staticmethod
    def build_path_to_init_data_dir(
            path_to_script: str = Constants.path.path_to_root_script_dir,
            init_dir: str = Constants.path.INIT_DATA_DIR) -> str:

        return os.path.join(path_to_script, init_dir)

    @classmethod
    def build_path_to_init_data_file(
            cls,
            structure_folder: str | None = None,
            path_to_script: str = Constants.path.path_to_root_script_dir,
            init_dir: str = Constants.path.INIT_DATA_DIR,
            file: str = "ljout.dat") -> str:

        path_dir: str = cls.build_path_to_init_data_dir(path_to_script, init_dir)

        if structure_folder is not None:
            return os.path.join(path_dir, structure_folder, file)
        return os.path.join(path_dir, file)

    ### RESULT DATA ###

    @staticmethod
    def build_path_to_result_data_dir(
            path_to_script: str = Constants.path.path_to_root_script_dir,
            result_dir: str = Constants.path.RESULT_DATA_DIR) -> str:

        path_to_result_dir: str = os.path.join(path_to_script, result_dir)

        # Create directory if it doesn't exist
        os.makedirs(path_to_result_dir, exist_ok=True)

        return path_to_result_dir

    @classmethod
    def build_path_to_result_data_file(
            cls,
            structure_folder: str,
            path_to_script: str = Constants.path.path_to_root_script_dir,
            result_dir: str = Constants.path.RESULT_DATA_DIR,
            file: str = Constants.path.INIT_PDB_FILE) -> str:

        path_to_result_dir: str = cls.build_path_to_result_data_dir(path_to_script, result_dir)

        path_to_data_in_result_dir: str = os.path.join(path_to_result_dir, structure_folder)

        # Create directory if it doesn't exist
        os.makedirs(path_to_data_in_result_dir, exist_ok=True)

        return os.path.join(path_to_data_in_result_dir, file)
