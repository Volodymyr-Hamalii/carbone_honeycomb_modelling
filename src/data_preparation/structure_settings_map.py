class ChannelLimits:
    # TODO: to automate searching for limits instead of providing them in the file

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


class ChannelPoints:
    def __init__(self, points: list[list[float]], direction: int) -> None:
        self.points: list[list[float]] = points
        self.direction: int = direction


class StructureSettings:
    def __init__(self, structure_setting: dict) -> None:
        self.channel_limits = ChannelLimits(
            x_min=structure_setting["channel_limits"]["x_min"],
            x_max=structure_setting["channel_limits"]["x_max"],
            y_min=structure_setting["channel_limits"]["y_min"],
            y_max=structure_setting["channel_limits"]["y_max"],
            z_min=structure_setting["channel_limits"].get("z_min"),
            z_max=structure_setting["channel_limits"].get("z_max"),
        )

        self.points_to_set_channel_planes: list[ChannelPoints] = [
            ChannelPoints(
                points=points_data["points"],
                direction=points_data["direction"],
            )
            for points_data
            in structure_setting.get("points_to_set_channel_planes", [])
        ]

        self.distance_from_plane: float = structure_setting.get("distance_from_plane", 0)
        self.max_distance_to_carbone_atoms: float = structure_setting.get("max_distance_to_carbone_atoms", 0)
