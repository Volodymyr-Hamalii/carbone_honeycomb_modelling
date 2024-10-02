import logging
import os
from pathlib import Path
import dotenv

dotenv.load_dotenv()


class PathConstants:
    # path_to_root_script = os.path.join()

    path_to_utils_dir: Path = Path(__file__).resolve().parent
    path_to_root_script_dir: Path = path_to_utils_dir.parent.parent


class FilenamesConstants:
    # filenames
    INIT_DATA_DIR: str = "init_data"
    RESULT_DATA_DIR: str = "result_data"

    INIT_DAT_FILE: str = "ljout.dat"
    INIT_PDB_FILE: str = "ljout-from-init-dat.pdb"
    PDB_FILE_ONE_CHANNEL: str = "ljout-from-init-dat-one-channel.pdb"

    STRUCTURE_SETTINGS_FILE: str = "structure_settings.json"

    AL_FILE: str = "al.pdb"


class LoggerConstants:
    LEVELS: dict[str, int] = {
        "debug": logging.DEBUG,  # 10
        "info": logging.INFO,  # 20
        "performance": 25,
        "metrics": 27,
        "warning": logging.WARNING,  # 30
        "error": logging.ERROR,  # 40
    }
    DEFAULT_LEVEL: str = "info"  # set 'warning' as default
    LEVEL: int = int(os.environ.get("LEVEL", 0)) or LEVELS[DEFAULT_LEVEL]


class PhysicsConstants:
    AL_LATTICE_PARAM = 4.0414


class Constants:
    path = PathConstants
    logger = LoggerConstants
    filenames = FilenamesConstants
    physics = PhysicsConstants

    MAX_NUMBER_OF_THREADS = 1
    DEFAULT_ACTION: str = os.environ.get("DEFAULT_ACTION") or "full_flow"
    DEFAULT_STRUCTURE_FOLDER: str = os.environ.get("DEFAULT_STRUCTURE_FOLDER") or "A1-7_h3"

    DEV_MODE: bool = os.environ.get("DEV_MODE", "false") == "true"
