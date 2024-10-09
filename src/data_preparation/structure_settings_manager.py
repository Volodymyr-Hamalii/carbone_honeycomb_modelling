import json
from pathlib import Path

from ..utils import PathBuilder, FileReader, Constants, Logger
from .structure_settings_map import StructureSettings, ChannelLimits, ChannelPoints


logger = Logger("StructureSettingsManager")


class StructureSettingsManager:
    file_name: str = Constants.filenames.STRUCTURE_SETTINGS_FILE

    @classmethod
    def read_file(cls, structure_folder: str) -> None | StructureSettings:
        structure_settings: dict | None = FileReader.read_json_file(
            structure_folder=structure_folder, file_name=cls.file_name)  # type: ignore

        if structure_settings is None:
            logger.warning(f"Settings file for {structure_folder} not found")
            return

        # TODO: to automate searching for limits instead of providing them in the file
        channel_limits: dict = structure_settings["channel_limits"]
        return StructureSettings(
            channel_limits=ChannelLimits(
                x_min=channel_limits["x_min"],
                x_max=channel_limits["x_max"],
                y_min=channel_limits["y_min"],
                y_max=channel_limits["y_max"],
                z_min=channel_limits.get("z_min"),
                z_max=channel_limits.get("z_max"),
            ),

            points_to_set_channel_planes=[
                ChannelPoints(
                    points=points_data["points"],
                    direction=points_data["direction"],
                ) for points_data in structure_settings.get("points_to_set_channel_planes", [])
            ],

            distance_from_plane=structure_settings.get("distance_from_plane", 0),
            max_distance_to_carbon_atoms=structure_settings.get("max_distance_to_carbon_atoms", 0),
            al_lattice_parameter=structure_settings.get("al_lattice_parameter", 0),
        )

    @classmethod
    def create_structure_settings_template(cls, structure_folder: str) -> None:
        structure_settings_template: dict = {
            "channel_limits": {
                "x_min": None,
                "x_max": None,
                "y_min": None,
                "y_max": None,
                "z_min": None,
                "z_max": None,
            },
            "points_to_set_channel_planes": [
                {
                    "points": [
                        [None, None, None],
                        [None, None, None],
                        [None, None, None]
                    ],
                    "direction": 1
                },
            ],
            "max_distance_to_carbon_atoms": None,
            "al_lattice_parameter": Constants.physics.AL_LATTICE_PARAM,
        }

        path_to_file: Path = PathBuilder.build_path_to_result_data_file(
            structure_folder=structure_folder, file=cls.file_name)

        if path_to_file.exists():
            logger.info(f"File {path_to_file} already exists.")
            return

        # Create a structure_settings.json template
        Path(path_to_file).write_text(json.dumps(structure_settings_template), encoding="utf-8")
        logger.info(f"Created structure_settings template in {path_to_file}. Please, fill it out.")
