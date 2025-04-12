from ..utils import FileReader, Constants, Logger
from .structure_settings import StructureSettings


logger = Logger("StructureSettingsManager")


class StructureSettingsManager:
    file_name: str = Constants.file_names.STRUCTURE_SETTINGS_FILE

    @classmethod
    def get_structure_settings(cls, structure_folder: str) -> StructureSettings:
        """ Read structure_settings file, parse it and calculate some fields. """

        structure_settings: dict | None = FileReader.read_json_file(
            folder_path=Constants.path.RESULT_DATA_PATH / structure_folder, file_name=cls.file_name)

        if structure_settings is None:
            # Read default structure settings
            structure_settings: dict | None = FileReader.read_json_file(
                folder_path=Constants.path.RESULT_DATA_PATH, file_name=cls.file_name)

        if structure_settings:
            distance_from_plane: float = structure_settings.get("distance_from_plane", 0)
            max_distance_to_carbon_atoms: float = structure_settings.get("max_distance_to_carbon_atoms", 0)
        else:
            logger.warning(f"structure_settings.json file not found for {structure_folder} structure.")
            distance_from_plane: float = 0
            max_distance_to_carbon_atoms: float = 0

        return StructureSettings(
            distance_from_plane=distance_from_plane,
            max_distance_to_carbon_atoms=max_distance_to_carbon_atoms,
        )
