import json
from pathlib import Path

from ..utils import PathBuilder, FileReader, Constants, Logger
from .structure_settings import StructureSettings, ChannelLimits, ChannelPoints


logger = Logger("StructureSettingsManager")


class StructureSettingsManager:
    file_name: str = Constants.filenames.STRUCTURE_SETTINGS_FILE

    @classmethod
    def get_structure_settings(cls, structure_folder: str) -> None | StructureSettings:
        """ Read structure_settings file, parse it and calculate some fields. """

        structure_settings: dict | None = FileReader.read_json_file(
            structure_folder=structure_folder, file_name=cls.file_name)  # type: ignore

        if structure_settings is None:
            logger.warning(f"Settings file for {structure_folder} not found")
            return

        points_to_set_channel_planes: list[ChannelPoints] = [
            ChannelPoints(
                points=points_data["points"],
                direction=points_data["direction"],
            ) for points_data in structure_settings.get("points_to_set_channel_planes", [])
        ]

        # Get channel_limits
        channel_limits_from_channel_points: ChannelLimits = cls.build_channel_limits_from_channel_points(
            points_to_set_channel_planes)

        channel_limits_data: dict = structure_settings.get("channel_limits", {})

        channel_limits = ChannelLimits(
            x_min=channel_limits_data.get("x_min") or channel_limits_from_channel_points.x_min,
            x_max=channel_limits_data.get("x_max") or channel_limits_from_channel_points.x_max,
            y_min=channel_limits_data.get("y_min") or channel_limits_from_channel_points.y_min,
            y_max=channel_limits_data.get("y_max") or channel_limits_from_channel_points.y_max,
            z_min=channel_limits_data.get("z_min") or channel_limits_from_channel_points.z_min,
            z_max=channel_limits_data.get("z_max") or channel_limits_from_channel_points.z_max,
        )

        return StructureSettings(
            points_to_set_channel_planes=points_to_set_channel_planes,
            channel_limits=channel_limits,
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
                    "direction": 1 if i % 2 else -1
                } for i in range(1, 7)
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
        Path(path_to_file).write_text(json.dumps(structure_settings_template, indent=4), encoding="utf-8")
        logger.info(f"Created structure_settings template in {path_to_file}. Please, fill it out.")

    @staticmethod
    def build_channel_limits_from_channel_points(channel_points: list[ChannelPoints]) -> ChannelLimits:
        all_points: list[list[float]] = [
            point
            for channel_points_item in channel_points
            for point in channel_points_item.points
        ]

        return ChannelLimits(
            x_min=min(point[0] for point in all_points),
            x_max=max(point[0] for point in all_points),
            y_min=min(point[1] for point in all_points),
            y_max=max(point[1] for point in all_points),
            z_min=min(point[2] for point in all_points),
            z_max=max(point[2] for point in all_points),
        )
