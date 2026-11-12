"""Temperature grids and stable case-directory names."""

import numpy as np


def temperature_points(values):
    """Expand inclusive [minimum, maximum, step] temperatures robustly."""
    minimum, maximum, step = (float(value) for value in values)
    count = int(np.floor((maximum - minimum) / step + 1.0e-10))
    return [minimum + index * step for index in range(count + 1)]


def temperature_directory_name(temperature):
    """Return a sortable, filesystem-safe temperature directory name."""
    if float(temperature).is_integer():
        return f"T{int(temperature):04d}K"
    value = f"{temperature:.6f}".rstrip("0").rstrip(".").replace(".", "p")
    return f"T{value}K"
