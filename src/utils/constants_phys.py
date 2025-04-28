from abc import ABC
from math import sqrt


__all__: list[str] = [
    "ConstantsAtomParams",
    "ConstantsAlParams",
    "ConstantsArParams",
    "ATOM_PARAMS_MAP",
]


class ConstantsAtomParams(ABC):
    """ Abstract class for atom params """
    ATOMS_NAME: str
    ATOM_SYMBOL: str

    LATTICE_PARAM: float
    DIST_BETWEEN_ATOMS: float
    DIST_BETWEEN_LAYERS: float

    MIN_RECOMENDED_DIST_BETWEEN_ATOMS: float
    MIN_ALLOWED_DIST_BETWEEN_ATOMS: float


class ConstantsAlParams(ConstantsAtomParams):
    """ Aluminium params """
    ATOMS_NAME: str = "Aluminium"
    ATOM_SYMBOL: str = "Al"

    LATTICE_PARAM: float = 4.049  # A
    DIST_BETWEEN_ATOMS: float = LATTICE_PARAM / sqrt(2)
    DIST_BETWEEN_LAYERS: float = LATTICE_PARAM / sqrt(3)

    MIN_RECOMENDED_DIST_BETWEEN_ATOMS: float = DIST_BETWEEN_ATOMS * 0.92
    MIN_ALLOWED_DIST_BETWEEN_ATOMS: float = DIST_BETWEEN_ATOMS * 0.7


class ConstantsArParams(ConstantsAtomParams):
    """ Argon params """
    ATOMS_NAME: str = "Argon"
    ATOM_SYMBOL: str = "Ar"

    LATTICE_PARAM: float = 5.310  # A
    DIST_BETWEEN_ATOMS: float = LATTICE_PARAM / sqrt(2)
    DIST_BETWEEN_LAYERS: float = LATTICE_PARAM / sqrt(3)

    MIN_RECOMENDED_DIST_BETWEEN_ATOMS: float = DIST_BETWEEN_ATOMS * 0.92
    MIN_ALLOWED_DIST_BETWEEN_ATOMS: float = DIST_BETWEEN_ATOMS * 0.7


ATOM_PARAMS_MAP: dict = {
    "al": ConstantsAlParams,  # Aluminium
    "ar": ConstantsArParams,  # Argon
}
