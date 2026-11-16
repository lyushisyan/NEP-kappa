"""Spatial plotting reads saved occupations and keeps the initial probe cells fixed."""
import json
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

from nepkappa.tdbte import H_THz, KB
from nepkappa.tdbte_observables import Observables
from nepkappa.tdbte_plot import plot_tdbte_results


@pytest.fixture
def spatial_plot_data(tmp_path):
    folder = tmp_path/'tdbte'
    folder.mkdir()
    frequency = np.array([0., 1., 2., 3.])
    ids = np.array([2, 2, 7, 7])
    # The hottest and coldest cells swap later, deliberately testing fixed probes.
    branch_T = np.array([
        [[310, 310], [302, 302], [290, 290], [300, 300]],
        [[305, 306], [302, 303], [299, 300], [301, 302]],
        [[295, 298], [303, 305], [310, 313], [300, 304]],
    ], dtype=float)
    n = np.zeros((3, 4, 4))
    for mode in range(1, 4):
        temperature = branch_T[:, :, 0 if ids[mode] == 2 else 1]
        n[:, :, mode] = 1/np.expm1(H_THz*frequency[mode]/(KB*temperature))
    kernel = SimpleNamespace(frequency=frequency, zero=frequency == 0)
    obs = Observables(kernel, ids, .5, 1., np.zeros(4))
    temperature = obs.temperature(n)
    with h5py.File(folder/'spatial.hdf5', 'w') as f:
        f.attrs.update(complete=True, uniform_mode_weight=.5)
        for name, value in dict(time_ps=[0., 1., 2.], x_nm=[1., 3., 5., 7.],
                                temperature_K=temperature, frequency_THz=frequency,
                                mode_branch_indices=ids, branch_indices=[2, 7], occupation=n).items():
            f[name] = value
    (folder/'audit.json').write_text(json.dumps({'numerical_checks_passed': True}))
    return tmp_path, branch_T, temperature


@pytest.mark.parametrize('layout,count', [('combined', 2), ('separate', 4), ('both', 6)])
def test_four_temperature_panels_and_physical_branch_values(spatial_plot_data, monkeypatch, layout, count):
    from matplotlib.figure import Figure
    root, branch_T, temperature = spatial_plot_data
    original = Figure.savefig
    combined = []
    def inspect(fig, path, **kwargs):
        if Path(path).name == 'tdbte_spatial.png':
            assert len(fig.axes) == 4
            assert all(not ax.texts for ax in fig.axes)  # no a/b/c/d annotations
            assert len(fig.axes[0].lines) == 3
            np.testing.assert_allclose(fig.axes[0].lines[0].get_ydata(), temperature[0])
            for line, cell in zip(fig.axes[1].lines, (0, 1, 2)):
                np.testing.assert_allclose(line.get_ydata(), temperature[:, cell])
            for ax, cell in zip(fig.axes[2:], (0, 2)):
                assert len(ax.lines) == 2
                for line, expected in zip(ax.lines, branch_T[:, cell].T):
                    np.testing.assert_allclose(line.get_ydata(), expected, atol=1e-7)
                    assert line.get_linestyle() == '-'
                    assert line.get_marker() == 'None'
            combined.append(True)
        return original(fig, path, **kwargs)
    monkeypatch.setattr(Figure, 'savefig', inspect)
    config = SimpleNamespace(result_dir=root, plot_layout=layout, plot_dpi=40)
    paths = plot_tdbte_results(config)
    assert len(paths) == count and all(Path(p).is_file() for p in paths)
    assert bool(combined) == (layout != 'separate')
    output = root/'plots'
    meta = json.loads((output/'tdbte_spatial_plot_metadata.json').read_text())
    assert meta['probe_positions_nm'] == {'hot': 1., 'intermediate': 3., 'cold': 5.}
    for name, cell in (('hot', 0), ('cold', 2)):
        table = np.loadtxt(output/f'tdbte_spatial_{name}_branch_temperature.csv', delimiter=',', skiprows=1)
        np.testing.assert_allclose(table[:, 1:], branch_T[:, cell], atol=1e-7)


def test_spatial_plot_honours_snapshot_times_and_requires_complete_audit(spatial_plot_data):
    root, _, _ = spatial_plot_data
    config = SimpleNamespace(result_dir=root, plot_layout='combined', plot_dpi=40,
                             plot_tdbte_times=[0., 1., 2.])
    plot_tdbte_results(config)
    meta = json.loads((root/'plots/tdbte_spatial_plot_metadata.json').read_text())
    assert meta['snapshot_times_ps'] == [0., 1., 2.]
    with h5py.File(root/'tdbte/spatial.hdf5', 'r+') as f:
        f.attrs['complete'] = False
    with pytest.raises(ValueError, match='incomplete'):
        plot_tdbte_results(config)
    (root/'tdbte/audit.json').write_text(json.dumps({'numerical_checks_passed': False}))
    with pytest.raises(ValueError, match='audit'):
        plot_tdbte_results(config)
