from src.utils import Constants


class VMParamsSetter:
    def __init__(self) -> None:
        # self.structure_folder: str = structure_folder

        self.to_build_bonds: bool = False
        self.to_show_coordinates: bool = False
        self.to_show_indexes: bool = False

        # Channel details
        self.to_show_dists_to_plane: bool = True
        self.to_show_dists_to_edges: bool = False
        self.to_show_channel_angles: bool = True

        self.excel_file_name: str = Constants.filenames.AL_FULL_CHANNEL_COORDINATES_XLSX_FILE
        self.dat_file_name: str = Constants.filenames.INIT_DAT_FILE

    # def set_structure_folder(self, value: str) -> None:
    #     self.structure_folder: str = value

    def set_to_build_bonds(self, value: bool) -> None:
        self.to_build_bonds: bool = value

    def set_to_show_coordinates(self, value: bool) -> None:
        self.to_show_coordinates: bool = value

    def set_to_show_indexes(self, value: bool) -> None:
        self.to_show_indexes: bool = value

    def set_excel_file_name(self, value: str) -> None:
        self.excel_file_name: str = value

    def set_dat_file_name(self, value: str) -> None:
        self.dat_file_name: str = value

    def set_to_show_dists_to_plane(self, value: bool) -> None:
        self.to_show_dists_to_plane: bool = value

    def set_to_show_dists_to_edges(self, value: bool) -> None:
        self.to_show_dists_to_edges: bool = value

    def set_to_show_channel_angles(self, value: bool) -> None:
        self.to_show_channel_angles: bool = value
