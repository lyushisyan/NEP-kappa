from types import SimpleNamespace

import h5py
import numpy as np


from nepkappa.transport_data import read_transport_data


def test_transport_hdf5_numerical_regression(tmp_path):
    """Lock the HDF5 interpretation and SI unit conversion to reference values."""
    path = tmp_path / "kappa-m211.hdf5"
    with h5py.File(path, "w") as handle:
        handle["temperature"] = np.array([100.0, 300.0])
        handle["kappa"] = np.array([[10.0, 11.0, 12.0], [20.0, 21.0, 22.0]])
        handle["heat_capacity"] = np.array([[[2.0], [4.0]], [[3.0], [5.0]]])
        handle["group_velocity"] = np.zeros((2, 1, 3))
        handle["frequency"] = np.array([[1.0], [2.0]])
        handle["weight"] = np.array([1.0, 3.0])
        handle["mesh"] = np.array([2, 1, 1])
        handle["gamma"] = np.ones((2, 2, 1))

    config = SimpleNamespace(
        dimensionality=3,
        plot_temperature=300.0,
        plot_tau="total",
        plot_kappa="all",
    )
    data = read_transport_data(path, 10.0, np.diag([1.0, 2.0, 5.0]), config)

    np.testing.assert_allclose(
        data["volume_heat_capacity"],
        [1.1215236438e11, 1.4419589706e11],
        rtol=1.0e-10,
    )
    np.testing.assert_allclose(
        data["kappa"],
        [[10.0, 11.0, 12.0], [20.0, 21.0, 22.0]],
    )
    assert data["tau_temperature_index"] == 1
    assert data["geometry_correction"]["factor"] == 1.0
