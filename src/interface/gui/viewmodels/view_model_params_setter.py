class ViewModelParamsSetter:
    def __init__(self) -> None:
        # self.structure_folder: str = structure_folder

        self.to_build_bonds: bool = False
        self.to_show_coordinates: bool = False
        self.to_show_indexes: bool = False

    # def set_structure_folder(self, value: str) -> None:
    #     self.structure_folder: str = value

    def set_to_build_bonds(self, value: bool) -> None:
        self.to_build_bonds: bool = value

    def set_to_show_coordinates(self, value: bool) -> None:
        self.to_show_coordinates: bool = value

    def set_to_show_indexes(self, value: bool) -> None:
        self.to_show_indexes: bool = value
