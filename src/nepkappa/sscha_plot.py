"""Compare completed SSCHA auxiliary phonon models across temperatures."""
from __future__ import annotations

import copy
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import yaml

from nepkappa.plot import (
    build_plot_data, build_harmonic_reference, build_harmonic_properties, find_phonon_metadata,
    load_phono3py_yaml, make_phonopy, draw_dispersion, journal_style,
    THZ_ANGSTROM_TO_KM_PER_S,
)
from nepkappa.sscha import temperature_points, temperature_directory_name

FIGURES = ('dispersion', 'dos', 'heat_capacity', 'group_velocity')


def load_sscha_series(config):
    """Require every requested temperature, and evaluate each FC2 only at its T."""
    root = Path(config.result_dir).expanduser().resolve()
    cases = []
    for temperature in temperature_points(config.scph_temps):
        label = temperature_directory_name(temperature)
        base = root / 'qha-sscha' / label if config.qha_enabled else root
        directory = base / config.scph_workdir / label
        summary_path = directory / 'summary.yaml'
        if not summary_path.is_file():
            raise FileNotFoundError(f'Missing completed SSCHA result: {summary_path}')
        summary = yaml.safe_load(summary_path.read_text()) or {}
        if summary.get('status') != 'complete' or not np.isclose(
            float(summary.get('temperature', float('nan'))), temperature
        ):
            raise ValueError(f'Incomplete or mismatched SSCHA temperature: {summary_path}')
        metadata = find_phonon_metadata(directory)
        fc2 = directory / 'fc2.hdf5'
        if not fc2.is_file():
            fc2 = directory / 'force_constants.hdf5'
        for path in (metadata, fc2):
            if not path.is_file():
                raise FileNotFoundError(f'Missing SSCHA plotting data: {path}')
        selected = copy.copy(config)
        selected.temps = [float(temperature)]
        selected._plot_harmonic_reference = str(
            base / 'fc2.hdf5' if config.qha_enabled
            else Path(getattr(config, 'scph_initial_fc2', None) or base / 'fc2.hdf5')
        )
        data = load_phono3py_yaml(metadata)
        phonon = make_phonopy(data, fc2)
        if not cases:
            anchor_phonon = phonon
        current = build_plot_data(phonon, data.unitcell, None, selected)
        reference = build_harmonic_reference(phonon, selected)
        cases.append({'temperature': float(temperature), 'data': current,
                      'reference': reference, 'primitive_numbers': phonon.primitive.numbers})
    # Volumes can change in QHA. Use the same reduced-q path and the first
    # temperature's distance axis so equivalent high-symmetry points align.
    anchor = cases[0]
    anchor['harmonic_heat_capacity'] = load_harmonic_heat_capacity(config, anchor_phonon)
    for case in cases:
        if not np.array_equal(case['primitive_numbers'], anchor['primitive_numbers']):
            raise ValueError('Cannot overlay SSCHA results with different primitive species/order')
        for data in (case['data'], case['reference']):
            if data is not None:
                align_band_axis(data['band'], anchor['data']['band'])
    return cases


def load_harmonic_heat_capacity(config, anchor):
    """Use one static-volume FC2 for the entire reference Cv temperature grid."""
    root = Path(config.result_dir).expanduser().resolve()
    expected_volume = None
    if config.qha_enabled:
        summary_path = root / 'qha-summary.yaml'
        summary = yaml.safe_load(summary_path.read_text()) if summary_path.is_file() else {}
        points = [p for p in (summary or {}).get('points', [])
                  if np.isclose(p['volume_ratio'], 1., rtol=0, atol=1e-8)]
        if (summary or {}).get('status') != 'complete' or len(points) != 1:
            print('  - Static-volume heat-capacity reference unavailable: completed QHA ratio=1 point required')
            return None
        point = points[0]
        fc2 = root / 'qha-points' / Path(point['directory']).name / 'fc2.hdf5'
        if not fc2.is_file():
            fc2 = fc2.with_name('force_constants.hdf5')
        expected_volume = float(point['primitive_volume_A3'])
    else:
        fc2 = Path(getattr(config, 'scph_initial_fc2', None) or root / 'fc2.hdf5').expanduser().resolve()
    metadata = find_phonon_metadata(fc2.parent)
    if not fc2.is_file() or not metadata.is_file():
        print('  - Static-volume heat-capacity reference unavailable: FC2 and matching metadata required')
        return None
    data = load_phono3py_yaml(metadata)
    phonon = make_phonopy(data, fc2)
    original, current = phonon.primitive, anchor.primitive
    if not np.array_equal(original.numbers, current.numbers):
        raise ValueError('Static-volume heat-capacity reference has different primitive species/order')
    if expected_volume is not None:
        if not np.isclose(original.volume, expected_volume, rtol=1e-6, atol=1e-6):
            raise ValueError('Static-volume heat-capacity reference does not match the QHA ratio=1 volume')
    else:
        delta = current.scaled_positions - original.scaled_positions
        if (not np.allclose(current.cell, original.cell, rtol=1e-6, atol=1e-6)
                or not np.allclose(delta - np.rint(delta), 0, atol=1e-5)):
            print('  - Static-volume heat-capacity reference unavailable: primitive cells or mappings differ')
            return None
    print(f'  - Static-volume heat-capacity reference: {fc2}')
    # No DOS, bands, transport results, or force calculations are needed here.
    return build_harmonic_properties(phonon, config)


def align_band_axis(band, anchor):
    """Reject different q paths instead of overlaying unrelated segments."""
    if band['labels'] != anchor['labels'] or len(band['qpoints']) != len(anchor['qpoints']):
        raise ValueError('SSCHA temperatures have different band paths; choose one custom path')
    for points, expected in zip(band['qpoints'], anchor['qpoints']):
        if np.shape(points) != np.shape(expected) or not np.allclose(points, expected):
            raise ValueError('SSCHA temperatures have different q points; choose one custom path')
    band['distances'] = anchor['distances']


def draw_series(name, ax, cases, coupled=False):
    colors = plt.get_cmap('viridis')(np.linspace(0.12, 0.88, len(cases)))
    if name == 'dispersion':
        # Reuse standard labels/grid, removing only the temporary band curves.
        draw_dispersion(ax, cases[0]['data']['band'])
        # Vertical symmetry lines have two x values; identify band artists by
        # their multi-point x arrays, without relying on insertion order.
        for line in list(ax.lines):
            if len(line.get_xdata()) > 2:
                line.remove()
    for index, (case, color) in enumerate(zip(cases, colors)):
        data = case['data']
        reference = case['reference']
        if name == 'heat_capacity':
            continue
        if reference is not None and (coupled or index == 0):
            draw_model(name, ax, reference, color if coupled else '0.45', True)
        draw_model(name, ax, data, color, False)
    if name == 'heat_capacity':
        temps = [case['temperature'] for case in cases]
        cv = [float(case['data']['transport']['volume_heat_capacity'][0]) for case in cases]
        ax.scatter(temps, cv, c=colors, s=40, zorder=3)
        reference_cv = cases[0].get('harmonic_heat_capacity')
        if reference_cv is not None:
            ax.plot(reference_cv['temperature'], reference_cv['volume_heat_capacity'],
                    '--', color='0.45', zorder=1)
        ax.set_xlabel('Temperature (K)')
        ax.set_ylabel(r'Volumetric heat capacity (J m$^{-3}$ K$^{-1}$)')
        ax.set_title('SSCHA points: FC2 at each temperature', fontsize=10)
    elif name == 'dos':
        ax.set_xlabel('DOS')
        ax.set_ylabel('Frequency (THz)')
    elif name == 'group_velocity':
        ax.set_xlabel('Frequency (THz)')
        ax.set_ylabel('Group velocity (km/s)')
    ax.grid(color='0.9', linewidth=.8)
    handles = [Line2D([], [], color=color, marker='o' if name in ('heat_capacity', 'group_velocity') else None,
                      linestyle='None' if name in ('heat_capacity', 'group_velocity') else '-',
                      label=f'{case["temperature"]:g} K') for case, color in zip(cases, colors)]
    has_reference = (cases[0].get('harmonic_heat_capacity') is not None
                     if name == 'heat_capacity' else any(c['reference'] is not None for c in cases))
    if has_reference:
        handles.append(Line2D([], [], color='0.45', linestyle='None' if name == 'group_velocity' else '--',
                              marker='o' if name == 'group_velocity' else None, markerfacecolor='none',
                              label=('Harmonic at fixed volume' if name == 'heat_capacity'
                                     else 'Harmonic at QHA volumes' if coupled else 'Harmonic reference')))
    ax.legend(handles=handles, fontsize=8, loc='best')


def draw_model(name, ax, data, color, reference):
    style = '--' if reference else '-'
    if name == 'dispersion':
        for distances, frequencies in zip(data['band']['distances'], data['band']['frequencies']):
            ax.plot(distances, frequencies, color=color, linestyle=style,
                    linewidth=1. if reference else 1.5, alpha=.65 if reference else 1.,
                    zorder=1 if reference else 2)
    elif name == 'dos':
        ax.plot(data['dos']['total_dos'], data['dos']['frequency_points'], color=color,
                linestyle=style, linewidth=1. if reference else 1.5)
    elif name == 'group_velocity':
        transport = data['transport']
        frequency = transport['frequency']
        speed = np.linalg.norm(transport['group_velocity'] * THZ_ANGSTROM_TO_KM_PER_S, axis=2)
        mask = np.isfinite(frequency) & np.isfinite(speed) & (frequency > 0)
        ax.scatter(frequency[mask], speed[mask], s=7 if reference else 9,
                   facecolors='none' if reference else color, edgecolors=color,
                   linewidths=.4, alpha=.22 if reference else .32,
                   zorder=1 if reference else 2)


def plot_sscha_results(config):
    cases = load_sscha_series(config)
    directory = Path(config.result_dir).expanduser().resolve() / 'plots'
    directory.mkdir(parents=True, exist_ok=True)
    saved = []
    with journal_style():
        if config.plot_layout in ('both', 'separate'):
            for name in FIGURES:
                fig, ax = plt.subplots(figsize=(7.4, 4.8), constrained_layout=True)
                draw_series(name, ax, cases, config.qha_enabled)
                path = directory / f'{name}.png'
                fig.savefig(path, dpi=config.plot_dpi)
                plt.close(fig)
                saved.append(path)
        if config.plot_layout in ('both', 'combined'):
            fig, axes = plt.subplots(2, 2, figsize=(11.6, 9), constrained_layout=True)
            for name, ax in zip(FIGURES, axes.ravel()):
                draw_series(name, ax, cases, config.qha_enabled)
            path = directory / 'combined.png'
            fig.savefig(path, dpi=config.plot_dpi)
            plt.close(fig)
            saved.append(path)
    print('  - SSCHA temperatures: ' + ', '.join(f'{c["temperature"]:g} K' for c in cases))
    print('  - Heat capacity: each temperature uses its own auxiliary FC2 and volume')
    print(f'  - Generated {len(saved)} SSCHA plot files: {directory}')
    return saved
