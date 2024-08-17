import logging
import os
import dotenv

dotenv.load_dotenv()


class Path:
    INIT_DATA_DIR: str = "init_data"
    RESULT_DATA_DIR: str = "result_data"

    INIT_DAT_FILE: str = "ljout.dat"
    INIT_PDB_FILE: str = "ljout-result.pdb"

    # path_to_root_script = os.path.join()

    path_to_utils_dir = os.path.dirname(os.path.abspath(__file__))
    path_to_root_script_dir = os.path.abspath(os.path.join(path_to_utils_dir, "../.."))


class Logger:
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


class Constants:
    path = Path
    logger = Logger

    MAX_NUMBER_OF_THREADS = 1
    DEV_MODE: bool = os.environ.get("DEV_MODE", "false") == "true"
