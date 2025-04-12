import os
from pathlib import Path
from math import sqrt
import logging
import dotenv
import sys


def _clear_env_cache() -> None:
    env_vars_to_clear: list[str] = [
        "DEFAULT_ACTION",
        "DEFAULT_STRUCTURE_FOLDER",
        "open_calculated_atoms",
        "set_equidistant",
        "al_lattice_type",
        "number_of_planes",
        "number_of_layers",
        "number_of_min_distances",
        "try_to_reflect_al_atoms",
        "show_al_layers",
        "interactive_mode",
    ]

    for var in env_vars_to_clear:
        if var in os.environ:
            del os.environ[var]


# Remove cached env vars and load .env file
_clear_env_cache()
dotenv.load_dotenv()


class _ConstantsSettings:
    MULTI_THREADED: bool = False
    DEV_MODE: bool = os.environ.get("DEV_MODE", "false") == "true"  # False by default


class _ConstantsFilenames:
    # file_names
    INIT_DATA_DIR: str = "init_data"
    RESULT_DATA_DIR: str = "result_data"

    INIT_DAT_FILE: str = "ljout.dat"
    INIT_PDB_FILE: str = "ljout-from-init-dat.pdb"
    PDB_FILE_ONE_CHANNEL: str = "ljout-from-init-dat-one-channel.pdb"

    AL_PLANE_COORDINATES_XLSX_FILE: str = "al-plane-coordinates.xlsx"
    AL_CHANNEL_COORDINATES_XLSX_FILE: str = "al-channel-coordinates.xlsx"
    AL_FULL_CHANNEL_COORDINATES_XLSX_FILE: str = "al-full-channel-coordinates.xlsx"
    AL_ALL_CHANNELS_COORDINATES_XLSX_FILE: str = "al-all-channels-coordinates.xlsx"
    AL_ALL_CHANNELS_COORDINATES_DAT_FILE: str = "Al.dat"
    C_ALL_CHANNELS_COORDINATES_DAT_FILE: str = "C.dat"
    AL_CHANNEL_DETAILS_XLSX_FILE: str = "built-al-structure-details.xlsx"

    STRUCTURE_SETTINGS_FILE: str = "structure_settings.json"

    AL_FILE: str = "al.pdb"


class _ConstantsPath:
    if getattr(sys, 'frozen', False):
        # If the application is frozen (e.g., packaged with PyInstaller)
        UTILS_DIR_PATH: Path = Path(sys.executable).resolve().parent
        ROOT_DIR_PATH: Path = UTILS_DIR_PATH
    else:
        # If running in a normal Python environment
        UTILS_DIR_PATH: Path = Path(__file__).resolve().parent
        ROOT_DIR_PATH: Path = UTILS_DIR_PATH.parent.parent

    INIT_DATA_PATH: Path = ROOT_DIR_PATH / _ConstantsFilenames.INIT_DATA_DIR
    RESULT_DATA_PATH: Path = ROOT_DIR_PATH / _ConstantsFilenames.RESULT_DATA_DIR


class _ConstantsLogger:
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


class _ConstantsMath:
    COORDINATE_INDEX_MAP: dict[str, int] = {
        "x": 0,
        "y": 1,
        "z": 2,
    }

    INDEX_COORDINATE_MAP: dict[int, str] = {
        i: c for c, i in COORDINATE_INDEX_MAP.items()
    }


class _ConstantsAlParams:
    # Aluminium
    LATTICE_PARAM: float = 4.049  # A
    DIST_BETWEEN_ATOMS: float = LATTICE_PARAM / sqrt(2)

    MIN_ALLOWED_DIST_BETWEEN_ATOMS: float = DIST_BETWEEN_ATOMS * 0.7
    MIN_RECOMENDED_DIST_BETWEEN_ATOMS: float = DIST_BETWEEN_ATOMS * 0.92


class _ConstantsPhys:
    al = _ConstantsAlParams

    MIN_ALLOWED_DIST_BETWEEN_AL_C: float = 2.15


class Constants:
    path = _ConstantsPath
    logger = _ConstantsLogger
    file_names = _ConstantsFilenames
    phys = _ConstantsPhys
    math = _ConstantsMath
    settings = _ConstantsSettings

    DEFAULT_ACTION: str = os.environ.get("DEFAULT_ACTION") or "full_flow"
    DEFAULT_STRUCTURE_FOLDER: str = os.environ.get("DEFAULT_STRUCTURE_FOLDER") or "A1-7_h3"
