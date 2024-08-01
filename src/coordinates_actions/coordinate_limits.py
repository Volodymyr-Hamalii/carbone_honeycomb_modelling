class CoordinateLimits:
    """ Coordinates to limit structure (e.g. only one honeycomb channel) """

    def __init__(
            self,
            x_min: float,
            x_max: float,
            y_min: float,
            y_max: float,
            z_min: float | None = None,
            z_max: float | None = None,
    ) -> None:
        self.x_min: float = x_min
        self.x_max: float = x_max
        self.y_min: float = y_min
        self.y_max: float = y_max

        if z_min is not None:
            self.z_min: float = z_min

        if z_max is not None:
            self.z_max: float = z_max
