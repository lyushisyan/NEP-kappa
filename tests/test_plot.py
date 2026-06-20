from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nepkappa.plot import (
    ANGSTROM3_TO_M3,
    EV_TO_J,
    effective_geometry_correction,
    volume_heat_capacity,
)


def test_volume_heat_capacity_uses_supplied_primitive_volume():
    heat_capacity = np.array([[[2.0], [4.0]]])
    weight = np.array([1.0, 3.0])
    mesh = np.array([2, 1, 1])
    volume_angstrom3 = 10.0

    cv = volume_heat_capacity(heat_capacity, weight, mesh, volume_angstrom3)
    expected = ((2.0 * 1.0 + 4.0 * 3.0) / 2) * EV_TO_J / (
        volume_angstrom3 * ANGSTROM3_TO_M3
    )

    assert np.allclose(cv, [expected])


def test_film_effective_geometry_correction_uses_effective_thickness():
    config = SimpleNamespace(
        dimensionality=2,
        effective_thickness=10.0,
        vacuum_axis="z",
    )
    cell = np.diag([5.0, 6.0, 40.0])
    correction = effective_geometry_correction(config, cell, 5.0 * 6.0 * 40.0)

    assert correction["dimensionality"] == 2
    assert np.isclose(correction["factor"], 4.0)


def test_nanowire_effective_geometry_correction_uses_effective_area():
    config = SimpleNamespace(
        dimensionality=1,
        effective_area=20.0,
        periodic_axis="z",
    )
    cell = np.diag([5.0, 6.0, 40.0])
    correction = effective_geometry_correction(config, cell, 5.0 * 6.0 * 40.0)

    assert correction["dimensionality"] == 1
    assert np.isclose(correction["factor"], 1.5)
