from src.interface.cli.show_init_data import AppActionsShowInitData


class AppViewModel:
    def __init__(self) -> None:
        self.show_bonds = False

    def set_show_bonds(self, value: bool):
        self.show_bonds: bool = value

    def run_action(self, distance_threshold: float):
        AppActionsShowInitData.show_init_structure("example_structure", self.show_bonds)
