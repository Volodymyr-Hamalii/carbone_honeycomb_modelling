from pathlib import Path
from src.utils import Constants


class VMParamsSetter:
    def __init__(self) -> None:
        # Plot details
        self.to_build_bonds: bool = True
        self.to_show_coordinates: bool = False
        self.to_show_indexes: bool = False

        self.x_min: float | None = None
        self.x_max: float | None = None
        self.y_min: float | None = None
        self.y_max: float | None = None
        self.z_min: float | None = None
        self.z_max: float | None = None

        self.bonds_num_of_min_distances: int = 2
        self.bonds_skip_first_distances: int = 0

        # Channel details
        self.to_show_dists_to_plane: bool = True
        self.to_show_dists_to_edges: bool = True
        self.to_show_channel_angles: bool = True

        # Files and paths
        self.data_dir: Path = Constants.path.RESULT_DATA_PATH
        self.file_name: str = "None"
        self.file_format: str = "None"  # "xlsx", "dat", "pdb"
        self.available_formats: list[str] = ["xlsx", "dat", "pdb"]
        self.excel_file_name: str = Constants.filenames.AL_FULL_CHANNEL_COORDINATES_XLSX_FILE
        self.dat_file_name: str = Constants.filenames.AL_ALL_CHANNELS_COORDINATES_DAT_FILE
        self.pdb_file_name: str = Constants.filenames.PDB_FILE_ONE_CHANNEL

        # Intercalation and sorption
        self.number_of_planes: int = 1
        self.num_of_al_layers: int = 1
        self.to_translate_al: bool = True
        self.to_try_to_reflect_al_atoms: bool = True
        self.to_equidistant_al_points: bool = True
        self.to_filter_al_atoms: bool = True
        self.al_lattice_type: str = "FCC"

    ######### Plot details #########

    def set_to_build_bonds(self, value: bool) -> None:
        self.to_build_bonds: bool = value

    def set_to_show_coordinates(self, value: bool) -> None:
        self.to_show_coordinates: bool = value

    def set_to_show_indexes(self, value: bool) -> None:
        self.to_show_indexes: bool = value

    def set_x_min(self, value: float) -> None:
        self.x_min = value

    def set_x_max(self, value: float) -> None:
        self.x_max = value

    def set_y_min(self, value: float) -> None:
        self.y_min = value

    def set_y_max(self, value: float) -> None:
        self.y_max = value

    def set_z_min(self, value: float) -> None:
        self.z_min = value

    def set_z_max(self, value: float) -> None:
        self.z_max = value

    def set_bonds_num_of_min_distances(self, value: int) -> None:
        self.bonds_num_of_min_distances: int = value

    def set_bonds_skip_first_distances(self, value: int) -> None:
        self.bonds_skip_first_distances: int = value

    ######### Channel details #########

    def set_to_show_dists_to_plane(self, value: bool) -> None:
        self.to_show_dists_to_plane: bool = value

    def set_to_show_dists_to_edges(self, value: bool) -> None:
        self.to_show_dists_to_edges: bool = value

    def set_to_show_channel_angles(self, value: bool) -> None:
        self.to_show_channel_angles: bool = value

    ######### Files and paths #########

    def set_data_dir(self, value: str) -> None:
        if value == Constants.filenames.INIT_DATA_DIR:
            self.data_dir: Path = Constants.path.INIT_DATA_PATH
        elif value == Constants.filenames.RESULT_DATA_DIR:
            self.data_dir: Path = Constants.path.RESULT_DATA_PATH
        else:
            raise ValueError(f"Invalid data directory: {value}")

    def set_file_name(self, value: str) -> None:
        self.file_name: str = value

    def set_file_format(self, value: str) -> None:
        self.file_format: str = value

    def set_excel_file_name(self, value: str) -> None:
        self.excel_file_name: str = value

    def set_dat_file_name(self, value: str) -> None:
        self.dat_file_name: str = value

    def set_pdb_file_name(self, value: str) -> None:
        self.pdb_file_name: str = value

    ######### Intercalation and sorption #########

    def set_number_of_planes(self, value: int) -> None:
        self.number_of_planes: int = value

    def set_num_of_al_layers(self, value: int) -> None:
        self.num_of_al_layers: int = value

    def set_to_try_to_reflect_al_atoms(self, value: bool) -> None:
        self.to_try_to_reflect_al_atoms: bool = value

    def set_to_equidistant_al_points(self, value: bool) -> None:
        self.to_equidistant_al_points: bool = value

    def set_to_filter_al_atoms(self, value: bool) -> None:
        self.to_filter_al_atoms: bool = value

    def set_al_lattice_type(self, value: str) -> None:
        self.al_lattice_type: str = value
