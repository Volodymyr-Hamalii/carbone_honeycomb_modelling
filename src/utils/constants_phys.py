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
    LATTICE_PARAM: float
    DIST_BETWEEN_ATOMS: float

    MIN_ALLOWED_DIST_BETWEEN_ATOMS: float
    MIN_RECOMENDED_DIST_BETWEEN_ATOMS: float

    MIN_ALLOWED_DIST_TO_C: float


class ConstantsAlParams(ConstantsAtomParams):
    """ Aluminium params """
    LATTICE_PARAM: float = 4.049  # A
    DIST_BETWEEN_ATOMS: float = LATTICE_PARAM / sqrt(2)

    MIN_ALLOWED_DIST_BETWEEN_ATOMS: float = DIST_BETWEEN_ATOMS * 0.7
    MIN_RECOMENDED_DIST_BETWEEN_ATOMS: float = DIST_BETWEEN_ATOMS * 0.92

    MIN_ALLOWED_DIST_TO_C: float = 2.15


class ConstantsArParams(ConstantsAtomParams):
    """ Argon params """
    LATTICE_PARAM: float = 3.755  # A
    DIST_BETWEEN_ATOMS: float = LATTICE_PARAM / sqrt(2)

    MIN_ALLOWED_DIST_BETWEEN_ATOMS: float = DIST_BETWEEN_ATOMS * 0.7
    MIN_RECOMENDED_DIST_BETWEEN_ATOMS: float = DIST_BETWEEN_ATOMS * 0.92

    MIN_ALLOWED_DIST_TO_C: float = 2.0


ATOM_PARAMS_MAP: dict = {
    "al": ConstantsAlParams,  # Aluminium
    "ar": ConstantsArParams,  # Argon
}
