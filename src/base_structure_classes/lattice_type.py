from abc import ABC, abstractmethod

class LatticeType(ABC):
    def __init__(self, structure_type: str) -> None:
        """ Base class for different lattice structures """
        self.is_cell: bool = structure_type == "cell"  # Cubic cell
        self.is_fcc: bool = structure_type == "FCC"  # Face-centered cubic
        self.is_hcp: bool = structure_type == "HCP"  # Hexagonal close-packed

        # Check the valid structure types in derived classes
        if not self._is_valid_structure(structure_type):
            raise ValueError(f"Invalid structure type: {structure_type}. Available options are: {self._get_available_types()}")


    def _is_valid_structure(self, structure_type: str) -> bool:
        return structure_type in self._get_available_types()


    @staticmethod
    @abstractmethod
    def get_info() -> str:
        """ Abstract method to provide information about the lattice type """
        ...


    @staticmethod
    @abstractmethod
    def _get_available_types() -> list[str]:
        """ Abstract method to provide a list of available structure types """
        ...