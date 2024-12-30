from dataclasses import dataclass


@dataclass
class CoordinateLimits:
    x_min: float = -float("inf")
    x_max: float = float("inf")

    y_min: float = -float("inf")
    y_max: float = float("inf")

    z_min: float = -float("inf")
    z_max: float = float("inf")
