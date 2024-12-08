from typing import Protocol
import numpy as np


class CarbonHoneycombPlaneProtocol(Protocol):
    points: np.ndarray
    direction: bool
