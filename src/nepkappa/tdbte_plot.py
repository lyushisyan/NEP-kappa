"""Plots for completed time-dependent BTE trajectories."""

from __future__ import annotations

import os
import json
from pathlib import Path
import tempfile
import warnings

import numpy as np


def _read_branch_energy(path: Path):
    if not path.is_file():
        raise FileNotFoundError(f"TD-BTE branch-energy data not found: {path}")
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().strip().split(",")
    if len(header) < 2 or header[0] != "time_ps":
        raise ValueError("TD-BTE branch-energy.csv has an invalid header.")
    data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    if data.shape[0] < 2 or data.shape[1] != len(header):
        raise ValueError("TD-BTE branch-energy.csv has an invalid data shape.")
    if not np.isfinite(data).all() or not np.all(np.diff(data[:, 0]) > 0):
        raise ValueError("TD-BTE branch-energy.csv needs finite data and increasing time.")
    labels = []
    for name in header[1:]:
        if not name.startswith("branch_") or not name.endswith("_excess_eV_per_cell"):
            raise ValueError(f"Unrecognized TD-BTE branch-energy column: {name}")
        labels.append(name[len("branch_"):-len("_excess_eV_per_cell")])
    return data[:, 0], data[:, 1:], labels


def _plot_energy_only(config):
    """Draw branch-resolved energy transfer and the total-energy check."""
    result_dir = Path(config.result_dir).resolve()
    times, energies, labels = _read_branch_energy(
        result_dir / "tdbte" / "branch-energy.csv"
    )
    layout = config.plot_layout
    dpi = int(config.plot_dpi)
    if dpi <= 0:
        raise ValueError("plot.dpi must be a positive integer.")

    cache_dir = Path(tempfile.gettempdir()) / "nepkappa-matplotlib"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_dir = result_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    paths = []

    def draw_branches(ax):
        for index, label in enumerate(labels):
            ax.plot(times, energies[:, index], label=f"Branch {label}")
        ax.set(xlabel="Time (ps)", ylabel="Excess energy (eV / primitive cell)",
               title="Energy redistribution")
        ax.legend(fontsize="small", ncol=2)
        ax.grid(alpha=0.3)

    def draw_total(ax):
        total = energies.sum(axis=1)
        ax.plot(times, total, color="tab:purple")
        ax.set(xlabel="Time (ps)", ylabel="Excess energy (eV / primitive cell)",
               title="Total excess energy")
        ax.grid(alpha=0.3)

    if layout in {"separate", "both"}:
        for name, draw in (("tdbte_branch_energy", draw_branches),
                           ("tdbte_total_energy", draw_total)):
            figure, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
            try:
                draw(ax)
                path = plot_dir / f"{name}.png"
                figure.savefig(path, dpi=dpi)
                paths.append(path)
            finally:
                plt.close(figure)
    if layout in {"combined", "both"}:
        figure, axes = plt.subplots(1, 2, figsize=(12, 4.5),
                                   constrained_layout=True)
        try:
            draw_branches(axes[0])
            draw_total(axes[1])
            path = plot_dir / "tdbte_combined.png"
            figure.savefig(path, dpi=dpi)
            paths.append(path)
        finally:
            plt.close(figure)
    print("  - TD-BTE plots:")
    for path in paths:
        print(f"    {path}")
    return paths


def branch_temperatures(frequency, occupations, branch_ids, weight):
    """Invert each branch's thermal energy; fixed zero modes carry no energy."""
    from scipy.optimize import brentq
    from .tdbte import H_THz, KB

    labels = np.unique(branch_ids)
    temperatures = np.empty((len(labels), occupations.shape[1]))
    for index, branch in enumerate(labels):
        mask = (branch_ids == branch) & (frequency > 0)
        quantum = H_THz * frequency[mask]
        if not len(quantum):
            raise ValueError(f'Branch {branch} has no positive-frequency modes for temperature inversion')
        targets = weight * (quantum @ occupations[mask])

        def equilibrium_energy(temperature):
            if temperature == 0:
                return 0.
            x = quantum / (KB * temperature)
            return weight * np.sum(quantum * np.exp(-x) / (-np.expm1(-x)))

        for step, target in enumerate(targets):
            if not np.isfinite(target) or target < 0:
                raise ValueError('Invalid branch energy for equivalent temperature')
            if target == 0:
                temperatures[index, step] = 0.
                continue
            upper = 300.
            while equilibrium_energy(upper) < target:
                upper *= 2
                if not np.isfinite(upper):
                    raise ValueError('Unable to bracket branch temperature')
            temperatures[index, step] = brentq(
                lambda temperature: equilibrium_energy(temperature) - target,
                0., upper, xtol=1e-9,
            )
    return temperatures


def interpolate_mesh_energy(addresses, mesh, values, qpoints):
    """Periodic trilinear interpolation of frequency-sorted branch energies."""
    mesh = np.asarray(mesh, dtype=int)
    addresses = np.asarray(addresses)
    if mesh.shape != (3,) or np.any(mesh < 2):
        raise ValueError('TD-BTE dispersion needs a full diagonal mesh')
    nq = int(np.prod(mesh))
    if (addresses.shape != (nq, 3) or not np.issubdtype(addresses.dtype, np.integer)
            or values.shape[0] != nq):
        raise ValueError('TD-BTE grid addresses do not match trajectory modes')
    addresses = addresses % mesh
    if len(np.unique(addresses, axis=0)) != nq:
        raise ValueError('TD-BTE grid addresses contain duplicates')
    grid = np.empty(tuple(mesh) + values.shape[1:])
    grid[tuple(addresses.T)] = values
    scaled = (np.asarray(qpoints) % 1.) * mesh
    base = np.floor(scaled).astype(int)
    fraction = scaled - base
    result = np.zeros((len(scaled),) + values.shape[1:])
    for offset in np.ndindex(2, 2, 2):
        offset = np.asarray(offset)
        indices = (base + offset) % mesh
        weights = np.prod(np.where(offset, fraction, 1. - fraction), axis=1)
        result += grid[tuple(indices.T)] * weights.reshape((-1,) + (1,) * (values.ndim - 1))
    return result


def _snapshot_indices(times, requested):
    if requested is None:
        requested = [times[0], (times[0] + times[-1]) / 2, times[-1]]
    requested = np.asarray(requested, dtype=float)
    if (requested.shape != (3,) or not np.isfinite(requested).all()
            or np.any(np.diff(requested) <= 0)
            or requested[0] < times[0] or requested[-1] > times[-1]):
        raise ValueError('plot.tdbte-times needs three increasing times within the saved trajectory')
    indices = np.array([np.argmin(abs(times - value)) for value in requested])
    if len(np.unique(indices)) != 3:
        raise ValueError('TD-BTE snapshots resolve to the same saved time; save more samples')
    if not np.allclose(times[indices], requested, rtol=0, atol=1e-10):
        warnings.warn(f'TD-BTE snapshots use nearest saved times: {times[indices].tolist()} ps',
                      RuntimeWarning)
    return indices


def _load_dynamics(config, times, energies, csv_labels):
    from .tdbte import H_THz, KB
    from .tdbte_storage import sha256
    from .plot import load_phono3py_yaml, make_phonopy, make_band_paths

    root = Path(config.result_dir).resolve()
    audit = json.loads((root / 'tdbte' / 'audit.json').read_text())
    if audit.get('numerical_checks_passed') is not True:
        raise ValueError('Six-panel TD-BTE plots require a passed numerical audit')
    with np.load(root / 'tdbte' / 'trajectories.npz', allow_pickle=False) as data:
        frequency = data['frequency_THz'].copy()
        occupations = data['excited'].copy()
        branch_ids = data['mode_branch_indices'].copy()
        saved_times = data['time_ps'].copy()
        weight = float(data['uniform_mode_weight'])
    if (frequency.ndim != 1 or branch_ids.shape != frequency.shape
            or not np.issubdtype(branch_ids.dtype, np.integer)
            or occupations.shape != (len(frequency), len(times))
            or not np.array_equal(saved_times, times)
            or not np.isfinite(frequency).all() or np.any(frequency < 0)
            or not np.isfinite(occupations).all() or np.any(occupations < 0)
            or not np.isfinite(weight) or weight <= 0):
        raise ValueError('Malformed TD-BTE trajectory arrays')
    labels = np.unique(branch_ids)
    nb = len(labels)
    if nb < 4 or not np.array_equal(labels, np.arange(nb)):
        raise ValueError('Acoustic/optical panels require three acoustic branches and optical branches')
    if csv_labels != [str(b) for b in labels]:
        raise ValueError('Branch-energy labels differ from trajectory branches')
    temperature = float(audit['temperature_K'])
    if not np.isfinite(temperature) or temperature <= 0:
        raise ValueError('Invalid reference temperature in TD-BTE audit')
    reference = np.zeros_like(frequency)
    positive = frequency > 0
    x = H_THz * frequency[positive] / (KB * temperature)
    reference[positive] = np.exp(-x) / (-np.expm1(-x))
    excess = H_THz * frequency[:, None] * (occupations - reference[:, None])
    reconstructed = np.array([weight * excess[branch_ids == b].sum(axis=0) for b in labels]).T
    if not np.allclose(energies, reconstructed, rtol=1e-8, atol=1e-12):
        raise ValueError('Branch-energy CSV differs from saved occupations/reference temperature')

    # Kernels already record the exact q ordering and input identities. No new
    # collision events or dynamics are constructed while plotting.
    candidates = [root / 'tdbte-kernel' / 'manifest.json']
    if getattr(config, 'tdbte_kernel', None):
        candidates.append(Path(config.tdbte_kernel).expanduser())
    candidates.append(Path(audit['kernel_path']))
    metadata = None
    for path in candidates:
        if path.suffix == '.npz':
            path = path.with_suffix('.json')
        if path.is_file() and sha256(path) == audit.get('metadata_sha256'):
            metadata = json.loads(path.read_text())
            break
    if metadata is None:
        raise ValueError('Matching TD-BTE kernel metadata is required for dispersion colouring')
    mesh = np.asarray(metadata.get('mesh', []), dtype=int)
    addresses = np.asarray(metadata.get('grid_addresses', []))
    if mesh.shape != (3,) or len(frequency) != int(np.prod(mesh)) * nb:
        raise ValueError('Kernel mesh does not match the trajectory')
    if (not np.array_equal(branch_ids, np.tile(np.arange(nb), int(np.prod(mesh))))
            or metadata.get('mode_branch_indices') != branch_ids.tolist()
            or not np.isclose(weight, 1. / np.prod(mesh), rtol=1e-12)):
        raise ValueError('TD-BTE q/branch ordering or normalization is inconsistent')

    sources = [root]
    for source in (getattr(config, 'tdbte_force_constants', None), metadata.get('source')):
        if source:
            sources.append(Path(source).expanduser())
    expected = metadata.get('input_sha256', {})
    required = ('phono3py_disp.yaml', 'fc2.hdf5')
    source = next((p for p in sources if all(
        (p / name).is_file() and expected.get(name) == sha256(p / name) for name in required
    )), None)
    if source is None:
        raise FileNotFoundError('TD-BTE dispersion needs the original matching fc2.hdf5 and '
                                'phono3py_disp.yaml; keep them or copy them to result-dir')
    phonon = make_phonopy(load_phono3py_yaml(source / required[0]), source / required[1])
    # Check frequency ordering against the original grid before mapping colours.
    interpolate_mesh_energy(addresses, mesh, excess.reshape(-1, nb, len(times))[:, :, :1],
                            [[0., 0., 0.]])
    phonon.run_qpoints(addresses / mesh)
    check = np.asarray(phonon.get_qpoints_dict()['frequencies'])
    check[abs(check) <= 1e-4] = 0.
    if not np.allclose(check.ravel(), frequency, rtol=1e-7, atol=1e-5):
        raise ValueError('Dispersion FC2 frequencies do not match TD-BTE modes')
    paths, path_labels = make_band_paths(phonon.primitive, config)
    phonon.run_band_structure(paths, is_band_connection=False)
    band = phonon.get_band_structure_dict()
    band['labels'] = path_labels
    indices = _snapshot_indices(times, getattr(config, 'plot_tdbte_times', None))
    snapshots = excess[:, indices].reshape(-1, nb, 3)
    colours = [interpolate_mesh_energy(addresses, mesh, snapshots, q) for q in paths]
    return dict(band=band, colours=colours, indices=indices, labels=labels,
                temperatures=branch_temperatures(frequency, occupations, branch_ids, weight),
                reference_temperature=temperature, source=source)


def _draw_energy_dispersion(ax, band, colours, snapshot, norm, cmap):
    from matplotlib.collections import LineCollection
    ticks, tick_labels, previous = [], [], None
    for distances, frequency, labels, colour in zip(
        band['distances'], band['frequencies'], band['labels'], colours
    ):
        distances, frequency = np.asarray(distances), np.asarray(frequency)
        start, end = labels
        if not ticks:
            ticks.append(distances[0]); tick_labels.append(start)
        elif previous != start:
            if np.isclose(ticks[-1], distances[0]):
                tick_labels[-1] += '|' + start
            else:
                ticks.append(distances[0]); tick_labels.append(start)
        ticks.append(distances[-1]); tick_labels.append(end)
        previous = end
        for branch in range(frequency.shape[1]):
            points = np.column_stack((distances, frequency[:, branch]))
            segments = np.stack((points[:-1], points[1:]), axis=1)
            collection = LineCollection(segments, cmap=cmap, norm=norm, linewidth=2.)
            collection.set_array((colour[:-1, branch, snapshot] + colour[1:, branch, snapshot]) / 2)
            ax.add_collection(collection)
        ax.axvline(distances[-1], color='0.85', linewidth=.8, zorder=0)
    all_frequency = np.concatenate(band['frequencies'])
    span = max(float(np.ptp(all_frequency)), 1.)
    ax.set(xlim=(ticks[0], ticks[-1]),
           ylim=(float(all_frequency.min()) - .03 * span, float(all_frequency.max()) + .03 * span),
           ylabel='Frequency (THz)', xlabel='Wave vector')
    ax.set_xticks(ticks, tick_labels)
    ax.grid(axis='y', alpha=.25)


def plot_tdbte_results(config):
    """Plot three energy-coloured dispersions and branch energy/temperatures."""
    root = Path(config.result_dir).resolve()
    times, energies, csv_labels = _read_branch_energy(root / 'tdbte' / 'branch-energy.csv')
    if not (root / 'tdbte' / 'trajectories.npz').is_file():
        if getattr(config, 'plot_tdbte_times', None) is not None:
            raise FileNotFoundError('Six-panel TD-BTE plots require trajectories.npz and kernel metadata')
        warnings.warn('Only branch-energy.csv is available: writing legacy two-panel energy plots',
                      RuntimeWarning)
        return _plot_energy_only(config)
    dpi = int(config.plot_dpi)
    if dpi <= 0:
        raise ValueError('plot.dpi must be a positive integer')
    dynamics = _load_dynamics(config, times, energies, csv_labels)
    cache = Path(tempfile.gettempdir()) / 'nepkappa-matplotlib'
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR', str(cache))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable

    low = min(float(c.min()) for c in dynamics['colours'])
    high = max(float(c.max()) for c in dynamics['colours'])
    extent = max(abs(low), abs(high), np.finfo(float).tiny)
    norm = Normalize(-extent, extent) if low < 0 else Normalize(0, extent)
    cmap = 'coolwarm' if low < 0 else 'viridis'
    colourbar = ScalarMappable(norm=norm, cmap=cmap)
    reference = dynamics['reference_temperature']
    branch_colours = plt.get_cmap('tab10')

    def draw_branches(ax):
        for index, label in enumerate(dynamics['labels']):
            ax.plot(times, energies[:, index], color=branch_colours(index), label=f'Branch {label}')
        ax.set(xlabel='Time (ps)', ylabel='Excess energy (eV / primitive cell)',
               title='Branch energy redistribution')
        ax.legend(fontsize='small', ncol=2)
        ax.grid(alpha=.3)

    def draw_temperature(ax, branches, title):
        for branch in branches:
            ax.plot(times, dynamics['temperatures'][branch], color=branch_colours(branch),
                    label=f'Branch {branch}')
        ax.axhline(reference, color='0.5', linestyle='--', linewidth=1., label=f'{reference:g} K reference')
        ax.set(xlabel='Time (ps)', ylabel='Equivalent temperature (K)', title=title)
        ax.legend(fontsize='small')
        ax.grid(alpha=.3)

    def draw_panel(index, ax):
        if index < 3:
            _draw_energy_dispersion(ax, dynamics['band'], dynamics['colours'], index, norm, cmap)
            ax.set_title(f't = {times[dynamics["indices"][index]]:g} ps')
        elif index == 3:
            draw_branches(ax)
        elif index == 4:
            draw_temperature(ax, range(3), 'Acoustic branches')
        else:
            draw_temperature(ax, range(3, len(csv_labels)), 'Optical branches')

    output = root / 'plots'
    output.mkdir(parents=True, exist_ok=True)
    paths = []
    names = [f'tdbte_dispersion_{i + 1}' for i in range(3)] + [
        'tdbte_branch_energy', 'tdbte_acoustic_temperature', 'tdbte_optical_temperature']
    if config.plot_layout in {'separate', 'both'}:
        for index, name in enumerate(names):
            fig, ax = plt.subplots(figsize=(6.2, 4.6), constrained_layout=True)
            try:
                draw_panel(index, ax)
                if index < 3:
                    fig.colorbar(colourbar, ax=ax, label='Excess energy (eV / mode)')
                path = output / (name + '.png')
                fig.savefig(path, dpi=dpi)
                paths.append(path)
            finally:
                plt.close(fig)
    if config.plot_layout in {'combined', 'both'}:
        fig, axes = plt.subplots(2, 3, figsize=(17.4, 9.0), constrained_layout=True)
        try:
            for index, ax in enumerate(axes.ravel()):
                draw_panel(index, ax)
            fig.colorbar(colourbar, ax=axes[0, :].tolist(), location='top',
                         fraction=.055, pad=.06, label='Excess energy (eV / mode; interpolated on q path)')
            path = output / 'tdbte_combined.png'
            fig.savefig(path, dpi=dpi)
            paths.append(path)
        finally:
            plt.close(fig)
    np.savetxt(output / 'tdbte_branch_temperature.csv',
               np.column_stack((times, dynamics['temperatures'].T)), delimiter=',', comments='',
               header='time_ps,' + ','.join(f'branch_{b}_temperature_K' for b in dynamics['labels']))
    (output / 'tdbte_plot_metadata.json').write_text(json.dumps({
        'snapshot_times_ps': times[dynamics['indices']].tolist(),
        'colour_quantity': 'H*f*(n(t)-n_Bose(T0)) in eV/mode; no zero-point contribution',
        'colour_interpolation': 'periodic trilinear, frequency-sorted branches; visual interpolation only',
        'colour_limits_eV_per_mode': [norm.vmin, norm.vmax],
        'reference_temperature_K': reference,
        'branch_temperature': 'Bose energy-equivalent temperature on the same q mesh; fixed zero modes excluded',
        'force_constants_source': str(dynamics['source']),
    }, indent=2) + '\n')
    print('  - TD-BTE snapshots (ps): ' + ', '.join(f'{t:g}' for t in times[dynamics['indices']]))
    print('  - Dispersion colours interpolate mode excess energy on the q mesh; frequencies are fixed.')
    print('  - Branch temperatures are energy-equivalent Bose temperatures, not fitted distributions.')
    print('  - TD-BTE plots:')
    for path in paths:
        print(f'    {path}')
    return paths
