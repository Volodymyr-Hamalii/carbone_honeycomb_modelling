from dataclasses import dataclass
import numpy as np

@dataclass
class Points:
    """ Template for any class with points array as a property. """
    points: np.ndarray

    def __len__(self) -> int:
        return len(self.points)
