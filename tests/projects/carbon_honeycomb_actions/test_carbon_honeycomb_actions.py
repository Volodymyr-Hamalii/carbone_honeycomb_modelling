from itertools import combinations
import numpy as np
from collections import defaultdict

# from src.projects import CarbonHoneycombActions
from src.coordinate_operations import PointsOrganizer


groups_by_xy1 = {
    (np.float32(4), np.float32(1)): np.array([...]),
    (np.float32(7), np.float32(14)): np.array([...]),
    (np.float32(2), np.float32(9)): np.array([...]),
    (np.float32(0), np.float32(5)): np.array([...]),
    (np.float32(6), np.float32(-7)): np.array([...]),
    (np.float32(-3), np.float32(-1)): np.array([...]),
    (np.float32(1), np.float32(13)): np.array([...]),
}

groups_by_the_xy_lines_fact = PointsOrganizer._group_by_the_xy_lines(groups_by_xy1)


groups_by_the_xy_lines_expected = [
    {
        (np.float32(2), np.float32(9)): np.array([...]),
        (np.float32(-3), np.float32(-1)): np.array([...]),
        (np.float32(0), np.float32(5)): np.array([...]),
    },
    {
        (np.float32(6), np.float32(-7)): np.array([...]),
        (np.float32(1), np.float32(13)): np.array([...]),
        (np.float32(4), np.float32(1)): np.array([...]),
        (np.float32(2), np.float32(9)): np.array([...]),
    },
    # {
    #     (np.float32(7), np.float32(14)): np.array([...]),
    # },
]
